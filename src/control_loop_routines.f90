!> \file
!> \author Chris Bradley
!> \brief This module handles all control loop routines.
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

!> This module handles all control loop routines.
MODULE CONTROL_LOOP_ROUTINES

  USE BaseRoutines
  USE ControlLoopAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE KINDS
  USE PROBLEM_CONSTANTS
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines  
  USE STRINGS
  USE TYPES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

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
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE ControlLoop_CreateStart
    MODULE PROCEDURE CONTROL_LOOP_CREATE_START
  END INTERFACE ControlLoop_CreateStart

  INTERFACE ControlLoop_CreateFinish
    MODULE PROCEDURE CONTROL_LOOP_CREATE_FINISH
  END INTERFACE ControlLoop_CreateFinish

  INTERFACE CONTROL_LOOP_CURRENT_TIMES_GET
    MODULE PROCEDURE ControlLoop_CurrentTimesGet
  END INTERFACE CONTROL_LOOP_CURRENT_TIMES_GET

  INTERFACE ControlLoop_IterationsSet
    MODULE PROCEDURE CONTROL_LOOP_ITERATIONS_SET
  END INTERFACE ControlLoop_IterationsSet

  INTERFACE CONTROL_LOOP_LABEL_GET
    MODULE PROCEDURE CONTROL_LOOP_LABEL_GET_C
    MODULE PROCEDURE CONTROL_LOOP_LABEL_GET_VS
  END INTERFACE CONTROL_LOOP_LABEL_GET
  
  INTERFACE ControlLoop_LabelGet
    MODULE PROCEDURE CONTROL_LOOP_LABEL_GET_C
    MODULE PROCEDURE CONTROL_LOOP_LABEL_GET_VS
  END INTERFACE ControlLoop_LabelGet
  
  INTERFACE CONTROL_LOOP_LABEL_SET
    MODULE PROCEDURE CONTROL_LOOP_LABEL_SET_C
    MODULE PROCEDURE CONTROL_LOOP_LABEL_SET_VS
  END INTERFACE CONTROL_LOOP_LABEL_SET
  
  INTERFACE ControlLoop_LabelSet
    MODULE PROCEDURE CONTROL_LOOP_LABEL_SET_C
    MODULE PROCEDURE CONTROL_LOOP_LABEL_SET_VS
  END INTERFACE ControlLoop_LabelSet
  
  INTERFACE ControlLoop_MaximumIterationsSet
    MODULE PROCEDURE CONTROL_LOOP_MAXIMUM_ITERATIONS_SET
  END INTERFACE ControlLoop_MaximumIterationsSet

  INTERFACE ControlLoop_LoadOutputSet
    MODULE PROCEDURE CONTROL_LOOP_LOAD_OUTPUT_SET
  END INTERFACE ControlLoop_LoadOutputSet
  
  INTERFACE ControlLoop_NumberOfSubLoopsGet
    MODULE PROCEDURE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET
  END INTERFACE ControlLoop_NumberOfSubLoopsGet

  INTERFACE ControlLoop_NumberOfSubLoopsSet
    MODULE PROCEDURE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET
  END INTERFACE ControlLoop_NumberOfSubLoopsSet

  INTERFACE ControlLoop_OutputTypeGet
    MODULE PROCEDURE CONTROL_LOOP_OUTPUT_TYPE_GET
  END INTERFACE ControlLoop_OutputTypeGet

  INTERFACE ControlLoop_OutputTypeSet
    MODULE PROCEDURE CONTROL_LOOP_OUTPUT_TYPE_SET
  END INTERFACE ControlLoop_OutputTypeSet

  INTERFACE ControlLoop_SubLoopGet
    MODULE PROCEDURE CONTROL_LOOP_SUB_LOOP_GET
  END INTERFACE ControlLoop_SubLoopGet

  INTERFACE ControlLoop_SolversDestroy
    MODULE PROCEDURE CONTROL_LOOP_SOLVERS_DESTROY
  END INTERFACE ControlLoop_SolversDestroy

  INTERFACE ControlLoop_SolverEquationsDestroy
    MODULE PROCEDURE CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY
  END INTERFACE ControlLoop_SolverEquationsDestroy

  INTERFACE ControlLoop_TimesGet
    MODULE PROCEDURE CONTROL_LOOP_TIMES_GET
  END INTERFACE ControlLoop_TimesGet

  INTERFACE ControlLoop_TimesSet
    MODULE PROCEDURE CONTROL_LOOP_TIMES_SET
  END INTERFACE ControlLoop_TimesSet

  INTERFACE ControlLoop_TypeSet
    MODULE PROCEDURE CONTROL_LOOP_TYPE_SET
  END INTERFACE ControlLoop_TypeSet

  INTERFACE ControlLoop_TimeInputSet
    MODULE PROCEDURE CONTROL_LOOP_TIME_INPUT_SET
  END INTERFACE ControlLoop_TimeInputSet
  
  INTERFACE ControlLoop_TimeOutputSet
    MODULE PROCEDURE CONTROL_LOOP_TIME_OUTPUT_SET
  END INTERFACE ControlLoop_TimeOutputSet

  PUBLIC CONTROL_LOOP_NO_OUTPUT,CONTROL_LOOP_PROGRESS_OUTPUT,CONTROL_LOOP_TIMING_OUTPUT
  
  PUBLIC CONTROL_LOOP_CREATE_FINISH,CONTROL_LOOP_CREATE_START

  PUBLIC ControlLoop_CreateFinish,ControlLoop_CreateStart

  PUBLIC CONTROL_LOOP_CURRENT_TIMES_GET

  PUBLIC ControlLoop_CurrentTimesGet

  PUBLIC ControlLoop_CurrentTimeInformationGet
  
  PUBLIC ControlLoop_Destroy

  PUBLIC ControlLoop_FieldVariablesCalculate
  
  PUBLIC CONTROL_LOOP_ITERATIONS_SET

  PUBLIC ControlLoop_IterationsSet

  PUBLIC CONTROL_LOOP_LABEL_GET,CONTROL_LOOP_LABEL_SET

  PUBLIC ControlLoop_LabelGet,ControlLoop_LabelSet

  PUBLIC CONTROL_LOOP_MAXIMUM_ITERATIONS_SET

  PUBLIC ControlLoop_MaximumIterationsSet

  PUBLIC CONTROL_LOOP_LOAD_OUTPUT_SET

  PUBLIC ControlLoop_LoadOutputSet

  PUBLIC ControlLoop_AbsoluteToleranceSet
  
  PUBLIC ControlLoop_RelativeToleranceSet

  PUBLIC CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET,CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET
  
  PUBLIC CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET,CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET

  PUBLIC ControlLoop_NumberOfSubLoopsGet,ControlLoop_NumberOfSubLoopsSet

  PUBLIC CONTROL_LOOP_OUTPUT_TYPE_GET,CONTROL_LOOP_OUTPUT_TYPE_SET

  PUBLIC ControlLoop_OutputTypeGet,ControlLoop_OutputTypeSet

  PUBLIC ControlLoop_PreviousValuesUpdate

  PUBLIC CONTROL_LOOP_SUB_LOOP_GET

  PUBLIC ControlLoop_SubLoopGet

  PUBLIC CONTROL_LOOP_SOLVERS_DESTROY

  PUBLIC ControlLoop_SolversDestroy

  PUBLIC CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY

  PUBLIC ControlLoop_SolverEquationsDestroy

  PUBLIC CONTROL_LOOP_TIMES_GET,CONTROL_LOOP_TIMES_SET

  PUBLIC ControlLoop_TimesGet,ControlLoop_TimesSet

  PUBLIC CONTROL_LOOP_TYPE_SET

  PUBLIC ControlLoop_TypeSet
  
  PUBLIC CONTROL_LOOP_TIME_INPUT_SET

  PUBLIC ControlLoop_TimeInputSet

  PUBLIC CONTROL_LOOP_TIME_OUTPUT_SET

  PUBLIC ControlLoop_TimeOutputSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the process of creating a control loop
  RECURSIVE SUBROUTINE CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to finish.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    
    ENTERS("CONTROL_LOOP_CREATE_FINISH",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        !Finish the sub-loops first
        IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
          DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
            CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP2,ERR,ERROR,*999)
          ENDDO !loop_idx
        ENDIF
        !Finish this control loop
        CONTROL_LOOP%CONTROL_LOOP_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_CREATE_FINISH",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the process of creating a control loop for a problem.
  SUBROUTINE CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER, INTENT(INOUT) :: PROBLEM !<A pointer to the problem to initialise the control for.
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<On exit, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("CONTROL_LOOP_CREATE_START",ERR,ERROR,*998)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ASSOCIATED(CONTROL_LOOP)) THEN
        CALL FlagError("Control loop is already associated.",ERR,ERROR,*998)
      ELSE
        CALL ControlLoop_Initialise(problem%CONTROL_LOOP,err,error,*999)        
        PROBLEM%CONTROL_LOOP%PROBLEM=>PROBLEM
        PROBLEM%CONTROL_LOOP%CONTROL_LOOP_LEVEL=1
        CONTROL_LOOP=>PROBLEM%CONTROL_LOOP
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*998)
    ENDIF
              
    EXITS("CONTROL_LOOP_CREATE_START")
    RETURN
999 NULLIFY(CONTROL_LOOP)
998 ERRORSEXITS("CONTROL_LOOP_CREATE_START",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_CREATE_START

  !
  !================================================================================================================================
  !

  !>Gets the current time parameters for a time control loop. If the specified loop is not a time loop the next time loop up the chain will be used. \see OpenCMISS::cmfe_ControlLoop_CurrentTimesGet
  SUBROUTINE ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<The control loop to get the current times for
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
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<The control loop to get the time information for
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
    TYPE(CONTROL_LOOP_TYPE), POINTER :: parentLoop
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: timeLoop

    ENTERS("ControlLoop_CurrentTimeInformationGet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    IF(.NOT.controlLoop%CONTROL_LOOP_FINISHED) CALL FlagError("Control loop has not been finished.",err,error,*999)

    !Find a time loop from either the specified control loop or the next time loop up the chain.
    controlLoopLevel=controlLoop%CONTROL_LOOP_LEVEL
    parentLoop=>controlLoop
    DO levelIdx=controlLoopLevel,1,-1
      IF(controlLoopLevel==0) THEN
        CALL FlagError("Could not find a time loop for the specified control loop.",err,error,*999)
      ELSE
        IF(parentLoop%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          timeLoop=>parentLoop%TIME_LOOP
          IF(ASSOCIATED(timeLoop)) THEN
            currentTime=timeLoop%CURRENT_TIME
            timeIncrement=timeLoop%TIME_INCREMENT
            startTime=timeLoop%START_TIME
            stopTime=timeLoop%STOP_TIME
            currentIteration=timeLoop%ITERATION_NUMBER
            outputIteration=timeLoop%OUTPUT_NUMBER
          ELSE
            CALL FlagError("Control loop time loop is not associated.",err,error,*999)
          ENDIF
          EXIT
        ELSE
          parentLoop=>parentLoop%PARENT_LOOP
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

  !>Destroy a control loop \see OpenCMISS::cmfe_ControlLoop_Destroy
  SUBROUTINE ControlLoop_Destroy(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    CALL ControlLoop_Finalise(controlLoop,err,error,*999)
       
    EXITS("ControlLoop_Destroy")
    RETURN
999 ERRORSEXITS("ControlLoop_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a ControlLoopFieldVariableType and deallocates all memory
  SUBROUTINE ControlLoop_FieldVariableFinalise(controlLoopFieldVariable,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariableType), INTENT(INOUT) :: controlLoopFieldVariable !<The control loop field variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_FieldVariableFinalise",err,error,*999)

    NULLIFY(controlLoopFieldVariable%fieldVariable)
    controlLoopFieldVariable%timeDependence=0
    controlLoopFieldVariable%linearity=0
      
    EXITS("ControlLoop_FieldVariableFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariableFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a ControlLoopFieldVariableType
  SUBROUTINE ControlLoop_FieldVariableInitialise(controlLoopFieldVariable,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariableType), INTENT(INOUT) :: controlLoopFieldVariable !<The control loop field variable to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_FieldVariableInitialise",err,error,*999)

    NULLIFY(controlLoopFieldVariable%fieldVariable)
    controlLoopFieldVariable%timeDependence=0
    controlLoopFieldVariable%linearity=0
      
    EXITS("ControlLoop_FieldVariableInitialise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariableInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableInitialise

  !
  !================================================================================================================================
  !

  !>Adds a field variable to the list of control loop field variables
  SUBROUTINE ControlLoop_FieldVariableAdd(controlLoopFieldVariables,variableLinearity,variableTimeDependence,fieldVariable, &
    & err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariablesType), POINTER :: controlLoopFieldVariables !<The control loop field variables to add the field variable to
    INTEGER(INTG), INTENT(IN) :: variableLinearity !<The linearity of the field variable
    INTEGER(INTG), INTENT(IN) :: variableTimeDependence !<The time dependence of the field variable
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable !<A pointer to the field variable to add.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    LOGICAL :: found
    TYPE(ControlLoopFieldVariableType), ALLOCATABLE :: newFieldVariables(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: controlLoopVariable
   
    ENTERS("ControlLoop_FieldVariableAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoopFieldVariables)) CALL FlagError("Control loop field variables is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    !See if we already have this field variable in the list for the control loop
    found=.FALSE.
    DO variableIdx=1,controlLoopFieldVariables%numberOfFieldVariables
      controlLoopVariable=>controlLoopFieldVariables%fieldVariables(variableIdx)%fieldVariable
      IF(ASSOCIATED(fieldVariable,controlLoopVariable)) THEN
        found=.TRUE.
        EXIT
      ENDIF
    ENDDO !variableIdx
    IF(found) THEN
      !We have found the variable. Check if the time dependence and nonlinearity needs to be updated.
      IF(variableTimeDependence>controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence) &
        & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence=variableTimeDependence
      IF(variableLinearity>controlLoopFieldVariables%fieldVariables(variableIdx)%linearity) &
        & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence=variableLinearity
    ELSE
      !We have not found the field variable in the control loop list so add it to the list.
      !Reallocate the list to take the new variable
      ALLOCATE(newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new field variables.",err,error,*999)
      DO variableIdx=1,controlLoopFieldVariables%numberOfFieldVariables
        newFieldVariables(variableIdx)%fieldVariable=> &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%fieldVariable
        newFieldVariables(variableIdx)%timeDependence= &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence
        newFieldVariables(variableIdx)%linearity= &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%linearity
      ENDDO !variableIdx
      !Add in the field variable and it's information
      CALL ControlLoop_FieldVariableInitialise(newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1),err,error,*999)
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%fieldVariable=>fieldVariable
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%timeDependence=variableTimeDependence
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%linearity=variableLinearity
      !Move alloc the new list
      CALL MOVE_ALLOC(newFieldVariables,controlLoopFieldVariables%fieldVariables)
      controlLoopFieldVariables%numberOfFieldVariables=controlLoopFieldVariables%numberOfFieldVariables+1
    ENDIF
     
    EXITS("ControlLoop_FieldVariableAdd")
    RETURN
999 IF(ALLOCATED(newFieldVariables)) DEALLOCATE(newFieldVariables)
    ERRORSEXITS("ControlLoop_FieldVariableAdd",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableAdd

  !
  !================================================================================================================================
  !

  !>Calculate the dependent field variables involved in a control loop
  RECURSIVE SUBROUTINE ControlLoop_FieldVariablesCalculate(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to calculate the field variables for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx,matrixIdx,solverIdx,variableIdx,variableLinearity,variableTimeDependence
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop2
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver
    TYPE(FIELD_TYPE), POINTER :: field
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("ControlLoop_FieldVariablesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%NUMBER_OF_SUB_LOOPS>0) THEN
      !We have sub loops so recursively calculate the field variables for the underlying control loops
      DO loopIdx=1,controlLoop%NUMBER_OF_SUB_LOOPS
        controlLoop2=>controlLoop%SUB_LOOPS(loopIdx)%PTR
        CALL ControlLoop_FieldVariablesCalculate(controlLoop2,err,error,*999)
      ENDDO !loopIdx
      !Add all the variables from the sub-loops to this loop
      !Initialise this control loop field variables
      CALL ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*999)
      DO loopIdx=1,controlLoop%NUMBER_OF_SUB_LOOPS
        controlLoop2=>controlLoop%SUB_LOOPS(loopIdx)%ptr
        IF(.NOT.ASSOCIATED(controlLoop2%fieldVariables)) THEN
          localError="Field variables is not associated for sub loop index "// &
            & TRIM(NumberToVString(loopIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO variableIdx=1,controlLoop2%fieldVariables%numberOfFieldVariables
          fieldVariable=>controlLoop2%fieldVariables%fieldVariables(variableIdx)%fieldVariable
          variableLinearity=controlLoop2%fieldVariables%fieldVariables(variableIdx)%linearity
          variableTimeDependence=controlLoop2%fieldVariables%fieldVariables(variableIdx)%timeDependence
          CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
            & fieldVariable,err,error,*999)
        ENDDO !variableIdx
      ENDDO !loopIdx
    ELSE
      !There are no sub loops so process the solvers for this loop to calculate the field variables
      !Initialise this control loop field variables
      CALL ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      !Loop over the solvers
      DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
        !Get the solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
        solverEquations=>solver%SOLVER_EQUATIONS
!!TODO: Need to think about solvers that do not have solver equations.
        IF(ASSOCIATED(solverEquations)) THEN
          !If we have solver equations then find the variables.
!!TODO: could flag the linearity and time dependence of variable in equations.
          SELECT CASE(solverEquations%linearity)
          CASE(SOLVER_EQUATIONS_LINEAR)
            variableLinearity=CONTROL_LOOP_FIELD_VARIABLE_LINEAR
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            variableLinearity=CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR
          CASE DEFAULT
            localError="The solver equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          SELECT CASE(solverEquations%timeDependence)
          CASE(SOLVER_EQUATIONS_STATIC)
            variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_STATIC
          CASE(SOLVER_EQUATIONS_QUASISTATIC)
            variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC
          CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
            dynamicSolver=>solver%DYNAMIC_SOLVER
            IF(ASSOCIATED(dynamicSolver)) THEN
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE)  THEN
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE)  THEN
                  variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC
                ELSE
                  variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC
                ENDIF
              ELSE
                variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC
              ENDIF
            ELSE
              CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The solver equations time dependence type of "// &
              & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Get the solver mapping
          NULLIFY(solverMatrices)
          CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
          !Loop over the solver matrices
          DO matrixIdx=1,solverMatrices%NUMBER_OF_MATRICES
            !Loop over the field variables associated with the solver mapping
            DO variableIdx=1,solverMapping%VARIABLES_LIST(matrixIdx)%NUMBER_OF_VARIABLES
              fieldVariable=>solverMapping%VARIABLES_LIST(matrixIdx)%variables(variableIdx)%variable
              IF(ASSOCIATED(fieldVariable)) THEN
                CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
                  & fieldVariable,err,error,*999)
              ELSE
                localError="The field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
                  & " of matrix index "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ENDDO !variableIdx
          ENDDO !matrixIdx
          !Add in the RHS
          DO variableIdx=1,solverMapping%rhsVariablesList%NUMBER_OF_VARIABLES
            fieldVariable=>solverMapping%rhsVariablesList%variables(variableIdx)%variable
            IF(ASSOCIATED(fieldVariable)) THEN
              CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
                & fieldVariable,err,error,*999)
            ELSE
              localError="The field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
                & " of matrix index "//TRIM(NumberToVString(matrixIdx,"*",err,error))//" is not associated."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !equationSetIdx
        ENDIF
      ENDDO !solverIdx
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Control loop field variables:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Loop level = ",controlLoop%CONTROL_LOOP_LEVEL,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sub loop index = ",controlLoop%SUB_LOOP_INDEX,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",controlLoop%label,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Loop type = ",controlLoop%LOOP_TYPE,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Field Variables:",err,error,*999)
      IF(ASSOCIATED(controlLoop%fieldVariables)) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of field variables = ",controlLoop%fieldVariables% &
          & numberOfFieldVariables,err,error,*999)
        IF(ALLOCATED(controlLoop%fieldVariables%fieldVariables)) THEN
          DO variableIdx=1,controlLoop%fieldVariables%numberOfFieldVariables          
            fieldVariable=>controlLoop%fieldVariables%fieldVariables(variableIdx)%fieldVariable
            IF(ASSOCIATED(fieldVariable)) THEN
              field=>fieldVariable%field
              IF(ASSOCIATED(field)) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Variable index : ",variableIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable field user number = ",field%USER_NUMBER, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type = ",fieldVariable%VARIABLE_TYPE, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable linearity = ",controlLoop%fieldVariables% &
                  & fieldVariables(variableIdx)%linearity,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable time dependence = ",controlLoop%fieldVariables% &
                  & fieldVariables(variableIdx)%timeDependence,err,error,*999)
              ELSE
                localError="Field variable field is not associated for variable index "// &
                  & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              localError="Field variable is not associated for variable index "// &
                & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)              
            ENDIF
          ENDDO !variableIdx
        ELSE
          CALL FlagError("Control loop field variables is not allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of field variables = ",0,err,error,*999)
      ENDIF      
    ENDIF
    
    EXITS("ControlLoop_FieldVariablesCalculate")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariablesCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesCalculate

  !
  !================================================================================================================================
  !

  !>Finalises a ControlLoopFieldVariablesType and deallocates all memory
  SUBROUTINE ControlLoop_FieldVariablesFinalise(controlLoopFieldVariables,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariablesType), POINTER :: controlLoopFieldVariables !<The control loop field variables to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldVariableIdx
   
    ENTERS("ControlLoop_FieldVariablesFinalise",err,error,*999)

    IF(ASSOCIATED(controlLoopFieldVariables)) THEN
      DO fieldVariableIdx=1,SIZE(controlLoopFieldVariables%fieldVariables,1)
        CALL ControlLoop_FieldVariableFinalise(controlLoopFieldVariables%fieldVariables(fieldVariableIdx),err,error,*998)
      ENDDO !fieldVariableIdx
998   DEALLOCATE(controlLoopFieldVariables%fieldVariables)
      controlLoopFieldVariables%numberOfFieldVariables=0
    ENDIF
      
    EXITS("ControlLoop_FieldVariablesFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariablesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a ControlLoopFieldVariableType for a controlLoop
  SUBROUTINE ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<The control loop to initialsie the field variables for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
   
    ENTERS("ControlLoop_FieldVariablesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    
    IF(ASSOCIATED(controlLoop%fieldVariables)) CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,err,error,*999)
    ALLOCATE(controlLoop%fieldVariables,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate control loop field variables.",err,error,*999)
    controlLoop%fieldVariables%numberOfFieldVariables=0
      
    EXITS("ControlLoop_FieldVariablesInitialise")
    RETURN
999 CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_FieldVariablesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesInitialise

  !
  !=================================================================================================================================
  !

  !>Finalise a control loop and deallocate all memory
  RECURSIVE SUBROUTINE ControlLoop_Finalise(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop2
 
    ENTERS("ControlLoop_Finalise",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      !Finalise any sub control loops first
      IF(controlLoop%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loopIdx=1,controlLoop%NUMBER_OF_SUB_LOOPS
          controlLoop2=>controlLoop%SUB_LOOPS(loopIdx)%ptr
          CALL ControlLoop_Finalise(controlLoop2,err,error,*999)
        ENDDO !loop_idx
        DEALLOCATE(controlLoop%SUB_LOOPS)
      ENDIF
      !Finalise any solvers
      IF(ASSOCIATED(controlLoop%solvers)) CALL SOLVERS_DESTROY(controlLoop%solvers,err,error,*999)
      !Now finalise this control loop
      controlLoop%label=""
      CALL CONTROL_LOOP_SIMPLE_FINALISE(controlLoop%SIMPLE_LOOP,err,error,*999)
      CALL CONTROL_LOOP_FIXED_FINALISE(controlLoop%FIXED_LOOP,err,error,*999)
      CALL CONTROL_LOOP_LOAD_INCREMENT_FINALISE(controlLoop%LOAD_INCREMENT_LOOP,err,error,*999)
      CALL CONTROL_LOOP_TIME_FINALISE(controlLoop%TIME_LOOP,err,error,*999)
      CALL CONTROL_LOOP_WHILE_FINALISE(controlLoop%WHILE_LOOP,err,error,*999)
      CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,err,error,*999)
      DEALLOCATE(controlLoop)
    ENDIF
       
    EXITS("ControlLoop_Finalise")
    RETURN
999 ERRORSEXITS("ControlLoop_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Finalise

  !
  !=================================================================================================================================
  !

  !>Initialise a control loop.
  SUBROUTINE ControlLoop_Initialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("ControlLoop_Initialise",err,error,*998)

    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    
    ALLOCATE(controlLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate control loop.",err,error,*999)
    NULLIFY(controlLoop%problem)
    NULLIFY(controlLoop%PARENT_LOOP)
    controlLoop%CONTROL_LOOP_FINISHED=.FALSE.
    controlLoop%label=" "
    controlLoop%CONTROL_LOOP_LEVEL=0
    controlLoop%SUB_LOOP_INDEX=0
    controlLoop%outputType=CONTROL_LOOP_NO_OUTPUT
    NULLIFY(controlLoop%SIMPLE_LOOP)
    NULLIFY(controlLoop%FIXED_LOOP)
    NULLIFY(controlLoop%TIME_LOOP)
    NULLIFY(controlLoop%WHILE_LOOP)
    NULLIFY(controlLoop%LOAD_INCREMENT_LOOP)
    controlLoop%NUMBER_OF_SUB_LOOPS=0
    NULLIFY(controlLoop%fieldVariables)
    NULLIFY(controlLoop%solvers)
    controlLoop%LOOP_TYPE=PROBLEM_CONTROL_SIMPLE_TYPE
    CALL CONTROL_LOOP_SIMPLE_INITIALISE(controlLoop,err,error,*999)
              
    EXITS("ControlLoop_Initialise")
    RETURN
999 CALL ControlLoop_Finalise(controlLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a fixed control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_FIXED_FINALISE(FIXED_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER, INTENT(INOUT) :: FIXED_LOOP !<A pointer to the fixed control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_FIXED_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(FIXED_LOOP)) THEN
      DEALLOCATE(FIXED_LOOP)
    ENDIF
       
    EXITS("CONTROL_LOOP_FIXED_FINALISE")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_FIXED_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_FIXED_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a fixed loop for a control loop.
  SUBROUTINE CONTROL_LOOP_FIXED_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to initialise the fixed loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("CONTROL_LOOP_FIXED_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%FIXED_LOOP)) THEN
        CALL FlagError("The fixed loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%FIXED_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate fixed loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%FIXED_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%FIXED_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%FIXED_LOOP%START_ITERATION=1
        CONTROL_LOOP%FIXED_LOOP%STOP_ITERATION=100
        CONTROL_LOOP%FIXED_LOOP%ITERATION_INCREMENT=1
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("CONTROL_LOOP_FIXED_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_FIXED_FINALISE(CONTROL_LOOP%FIXED_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("CONTROL_LOOP_FIXED_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_FIXED_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the iteration parameters for a fixed control loop. \see OpenCMISS::cmfe_ControlLoop_IterationsSet
  SUBROUTINE CONTROL_LOOP_ITERATIONS_SET(CONTROL_LOOP,START_ITERATION,STOP_ITERATION,ITERATION_INCREMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to fixed control loop to set the iterations for
    INTEGER(INTG), INTENT(IN) :: START_ITERATION !<The start iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: STOP_ITERATION !<The stop iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: ITERATION_INCREMENT !<The iteration increment for the fixed control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_FIXED_TYPE), POINTER :: FIXED_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("CONTROL_LOOP_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_FIXED_LOOP_TYPE) THEN
          FIXED_LOOP=>CONTROL_LOOP%FIXED_LOOP
          IF(ASSOCIATED(FIXED_LOOP)) THEN
            IF(ITERATION_INCREMENT==0) THEN
              LOCAL_ERROR="The specified time increment of "//TRIM(NUMBER_TO_VSTRING(ITERATION_INCREMENT,"*",ERR,ERROR))// &
                & " is invalid. The iteration increment must not be zero."          
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            ELSE
              IF(ITERATION_INCREMENT>0) THEN
                IF(STOP_ITERATION<=START_ITERATION) THEN
                  LOCAL_ERROR="The specified stop iteration of "//TRIM(NUMBER_TO_VSTRING(STOP_ITERATION,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified start increment of "// &
                    & TRIM(NUMBER_TO_VSTRING(START_ITERATION,"*",ERR,ERROR))// &
                    & ". For a positive iteration increment the stop iteration must be > than the start iteration."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ELSE
                IF(START_ITERATION<=STOP_ITERATION) THEN
                  LOCAL_ERROR="The specified start iteration of "//TRIM(NUMBER_TO_VSTRING(START_ITERATION,"*",ERR,ERROR))// &
                    & " is incompatiable with a specified stop iteration of "// &
                    & TRIM(NUMBER_TO_VSTRING(STOP_ITERATION,"*",ERR,ERROR))// &
                    & ". For a negative iteration increment the stop iteration must be < than the start iteration."
                  CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                ENDIF
              ENDIF
            ENDIF
            FIXED_LOOP%START_ITERATION=START_ITERATION
            FIXED_LOOP%STOP_ITERATION=STOP_ITERATION
            FIXED_LOOP%ITERATION_INCREMENT=ITERATION_INCREMENT
          ELSE
            CALL FlagError("Control loop fixed loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The specified control loop is not a fixed control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_ITERATIONS_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_ITERATIONS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_ITERATIONS_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the label of a control loop. \see OpenCMISS::cmfe_ControlLoop_LabelGet
  SUBROUTINE CONTROL_LOOP_LABEL_GET_C(CONTROL_LOOP,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: LABEL !<On return, the control loop label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: C_LENGTH,VS_LENGTH

    ENTERS("CONTROL_LOOP_LABEL_GET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      C_LENGTH=LEN(LABEL)
      VS_LENGTH=LEN_TRIM(CONTROL_LOOP%LABEL)
      IF(C_LENGTH>VS_LENGTH) THEN
        LABEL=CHAR(CONTROL_LOOP%LABEL,VS_LENGTH)
      ELSE
        LABEL=CHAR(CONTROL_LOOP%LABEL,C_LENGTH)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("CONTROL_LOOP_LABEL_GET_C")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LABEL_GET_C",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_LABEL_GET_C

   !
  !================================================================================================================================
  !

  !>Returns the label of a control loop. \see OpenCMISS::cmfe_ControlLoop_LabelGet
  SUBROUTINE CONTROL_LOOP_LABEL_GET_VS(CONTROL_LOOP,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: LABEL !<On return, the control loop label.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("CONTROL_LOOP_LABEL_GET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      LABEL=VAR_STR(CHAR(CONTROL_LOOP%LABEL))
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("CONTROL_LOOP_LABEL_GET_VS")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LABEL_GET_VS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_LABEL_GET_VS

  !
  !================================================================================================================================
  !

  !>Sets the label of a control loop. \see OpenCMISS::cmfe_ControlLoop_LabelSet
  SUBROUTINE CONTROL_LOOP_LABEL_SET_C(CONTROL_LOOP,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("CONTROL_LOOP_LABEL_SET_C",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        CONTROL_LOOP%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("CONTROL_LOOP_LABEL_SET_C")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LABEL_SET_C",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_LABEL_SET_C

  !
  !================================================================================================================================
  !

  !>Sets the label of a control loop. \see OpenCMISS::cmfe_ControlLoop_LabelSet
  SUBROUTINE CONTROL_LOOP_LABEL_SET_VS(CONTROL_LOOP,LABEL,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: LABEL !<The label to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables

    ENTERS("CONTROL_LOOP_LABEL_SET_VS",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        CONTROL_LOOP%LABEL=LABEL
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("CONTROL_LOOP_LABEL_SET_VS")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LABEL_SET_VS",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_LABEL_SET_VS

  !
  !================================================================================================================================
  !

  !>Sets the maximum number of iterations for a while or load increment control loop. \see OpenCMISS_cmfe_ControlLoop_MaximumIterationsSet
  SUBROUTINE CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(CONTROL_LOOP,MAXIMUM_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to the while or load incremented control loop to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_ITERATIONS !<The maximum number of iterations for the while or load incremented control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: WHILE_LOOP
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER :: LOAD_INCREMENT_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        !allow to update the maximum number of iterations at a later time for the load increment loop type
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
          LOAD_INCREMENT_LOOP=>CONTROL_LOOP%LOAD_INCREMENT_LOOP
          IF(ASSOCIATED(LOAD_INCREMENT_LOOP)) THEN
            IF(MAXIMUM_ITERATIONS<=0) THEN
              LOCAL_ERROR="The specified maximum number of iterations of "// &
                & TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                & " is invalid. The maximum number of iterations must be greater than zero."          
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)            
            ENDIF
            LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
          ELSE
            CALL FlagError("Control loop load incremented loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
        ENDIF
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
          WHILE_LOOP=>CONTROL_LOOP%WHILE_LOOP
          IF(ASSOCIATED(WHILE_LOOP)) THEN
            IF(MAXIMUM_ITERATIONS<=0) THEN
              LOCAL_ERROR="The specified maximum number of iterations of "// &
                & TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                & " is invalid. The maximum number of iterations must be greater than zero."          
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)            
            ENDIF
            WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
          ELSE
            CALL FlagError("Control loop while loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSEIF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
          LOAD_INCREMENT_LOOP=>CONTROL_LOOP%LOAD_INCREMENT_LOOP
          IF(ASSOCIATED(LOAD_INCREMENT_LOOP)) THEN
            IF(MAXIMUM_ITERATIONS<=0) THEN
              LOCAL_ERROR="The specified maximum number of iterations of "// &
                & TRIM(NUMBER_TO_VSTRING(MAXIMUM_ITERATIONS,"*",ERR,ERROR))// &
                & " is invalid. The maximum number of iterations must be greater than zero."          
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)            
            ENDIF
            LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=MAXIMUM_ITERATIONS
          ELSE
            CALL FlagError("Control loop load increment loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The specified control loop is not a while or load increment control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_MAXIMUM_ITERATIONS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_MAXIMUM_ITERATIONS_SET

  !
  !================================================================================================================================
  !

  !>Sets/changes the output for a load incremented control loop identified by an object. \see OpenCMISS_cmfe_ControlLoop_LoadOutputSet
  SUBROUTINE CONTROL_LOOP_LOAD_OUTPUT_SET(CONTROL_LOOP,OUTPUT_FREQUENCY,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to the load incremented control loop to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_FREQUENCY !<The output frequency modulo to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER :: LOAD_INCREMENT_LOOP
 
    ENTERS("CONTROL_LOOP_LOAD_OUTPUT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
          LOAD_INCREMENT_LOOP=>CONTROL_LOOP%LOAD_INCREMENT_LOOP
          IF(ASSOCIATED(LOAD_INCREMENT_LOOP)) THEN
            LOAD_INCREMENT_LOOP%OUTPUT_NUMBER=OUTPUT_FREQUENCY
          ELSE
            CALL FlagError("Control loop load increment loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The specified control loop is not a load increment control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_LOAD_OUTPUT_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LOAD_OUTPUT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_LOAD_OUTPUT_SET

  !
  !================================================================================================================================
  !

  !>Sets the absolute tolerance (convergence condition tolerance) for a while control loop. \see OpenCMISS::cmfe_ControlLoop_AbsoluteToleranceSet
  SUBROUTINE ControlLoop_AbsoluteToleranceSet(controlLoop,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<A pointer to while control loop to set the maximum iterations for
    REAL(DP), INTENT(IN) :: absoluteTolerance !<The absolute tolerance value for a while control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: whileLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_AbsoluteToleranceSet",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(controlLoop%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",err,error,*999)
      ELSE
        IF(controlLoop%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
          whileLoop=>controlLoop%WHILE_LOOP
          IF(ASSOCIATED(whileLoop)) THEN
            IF(absoluteTolerance<=0) THEN
              localError="The specified absolute tolerance of "// &
                & TRIM(NUMBER_TO_VSTRING(absoluteTolerance,"*",err,error))// &
                & " is invalid for a while loop. The tolerance must be greater than zero."          
              CALL FlagError(localError,err,error,*999)            
            ENDIF
            whileLoop%ABSOLUTE_TOLERANCE=absoluteTolerance
          ELSE
            CALL FlagError("Control loop while loop is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ControlLoop_AbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("ControlLoop_AbsoluteToleranceSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AbsoluteToleranceSet

  !
  !================================================================================================================================
  !

  !>Sets the relative tolerance (convergence condition tolerance) for a while control loop. \see OpenCMISS_cmfe_ControlLoopRelativeToleranceSet
  SUBROUTINE ControlLoop_RelativeToleranceSet(controlLoop,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<A pointer to while control loop to set the maximum iterations for
    REAL(DP), INTENT(IN) :: relativeTolerance !<The relative tolerance value for a while control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: whileLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_RelativeToleranceSet",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      IF(controlLoop%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",err,error,*999)
      ELSE
        IF(controlLoop%LOOP_TYPE==PROBLEM_CONTROL_WHILE_LOOP_TYPE) THEN
          whileLoop=>controlLoop%WHILE_LOOP
          IF(ASSOCIATED(whileLoop)) THEN
            IF(relativeTolerance<=0) THEN
              localError="The specified relative tolerance of "// &
                & TRIM(NUMBER_TO_VSTRING(relativeTolerance,"*",err,error))// &
                & " is invalid for a while loop. The tolerance must be greater than zero."          
              CALL FlagError(localError,err,error,*999)            
            ENDIF
            whileLoop%RELATIVE_TOLERANCE=relativeTolerance
          ELSE
            CALL FlagError("Control loop while loop is not associated.",err,error,*999)
          ENDIF
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
       
    EXITS("ControlLoop_RelativeToleranceSet")
    RETURN
999 ERRORSEXITS("ControlLoop_RelativeToleranceSet",err,error)
    RETURN 1
  END SUBROUTINE ControlLoop_RelativeToleranceSet

  !
  !================================================================================================================================
  !

  !>Sets the number of iterations for a time type control loop. If set to 0 (default), it will be computed by start and stop time and time increment. \see OPENCMISS_CMISSControlLoopNumberOfIterationsSet
  SUBROUTINE CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET(CONTROL_LOOP,NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to time control loop to set the number of iterations for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_ITERATIONS !<The number of iterations for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
        IF(ASSOCIATED(TIME_LOOP)) THEN
          IF(NUMBER_OF_ITERATIONS<0) THEN
            LOCAL_ERROR="The specified number of iterations of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_ITERATIONS,"*",ERR,ERROR))// &
              & " is invalid. The number must be non-negative."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ENDIF
          TIME_LOOP%NUMBER_OF_ITERATIONS=NUMBER_OF_ITERATIONS
          
          !Update time increment if number of iterations differs from time stepping settings
          IF (CEILING((TIME_LOOP%STOP_TIME-TIME_LOOP%START_TIME)/TIME_LOOP%TIME_INCREMENT) &
            & /= TIME_LOOP%NUMBER_OF_ITERATIONS) THEN
            TIME_LOOP%TIME_INCREMENT = (TIME_LOOP%STOP_TIME-TIME_LOOP%START_TIME)/TIME_LOOP%NUMBER_OF_ITERATIONS
          ENDIF
          
        ELSE
          CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_OF_ITERATIONS_SET
  
  !
  !================================================================================================================================
  !

  !>Gets the number of iterations for a time type control loop. If the value is not set to something /=0, it will be computed the first time the loop is executed. If it is retrieved earlier and 0 is returned, this means the value was not yet computed. \see OPENCMISS_CMISSControlLoopNumberOfIterationsGet
  SUBROUTINE CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET(CONTROL_LOOP,NUMBER_OF_ITERATIONS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to time control loop to set the number of iterations for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_ITERATIONS !<The number of iterations for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
 
    ENTERS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
        IF(ASSOCIATED(TIME_LOOP)) THEN
          NUMBER_OF_ITERATIONS=TIME_LOOP%NUMBER_OF_ITERATIONS
        ELSE
          CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_OF_ITERATIONS_GET
  
  !
  !================================================================================================================================
  !

  !>Gets the number of sub loops for a control loop. \see OpenCMISS::cmfe_ControlLoopNumberOfSubLoopsGet
  SUBROUTINE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET(CONTROL_LOOP,NUMBER_OF_SUB_LOOPS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to get the number of sub loops for
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_SUB_LOOPS !<On return, the number of sub loops for the specified control loop
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        NUMBER_OF_SUB_LOOPS=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
      ENDIF      
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_GET

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of sub loops in a control loop. \see OpenCMISS::cmfe_ControlLoop_NumberOfSubLoopsSet
  SUBROUTINE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,NUMBER_OF_SUB_LOOPS,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to control loop to set the number of sub loops for
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_SUB_LOOPS !<The number of sub loops to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_PTR_TYPE), ALLOCATABLE :: OLD_SUB_LOOPS(:)
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(NUMBER_OF_SUB_LOOPS>=0) THEN
          IF(NUMBER_OF_SUB_LOOPS/=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
            IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
              ALLOCATE(OLD_SUB_LOOPS(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS),STAT=ERR)
              IF(ERR/=0) CALL FlagError("Could not allocate old sub loops.",ERR,ERROR,*999)
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                OLD_SUB_LOOPS(loop_idx)%PTR=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DEALLOCATE(CONTROL_LOOP%SUB_LOOPS)
            ENDIF
            ALLOCATE(CONTROL_LOOP%SUB_LOOPS(NUMBER_OF_SUB_LOOPS),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate control loop sub loops.",ERR,ERROR,*999)
            IF(NUMBER_OF_SUB_LOOPS>CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
              DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR=>OLD_SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DO loop_idx=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS+1,NUMBER_OF_SUB_LOOPS
                NULLIFY(CONTROL_LOOP%SUB_LOOPS(loop_idx)%ptr)
                CALL ControlLoop_Initialise(CONTROL_LOOP%SUB_LOOPS(loop_idx)%ptr,err,error,*999)
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%PROBLEM=>CONTROL_LOOP%PROBLEM
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%PARENT_LOOP=>CONTROL_LOOP
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%CONTROL_LOOP_LEVEL=CONTROL_LOOP%CONTROL_LOOP_LEVEL+1
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR%SUB_LOOP_INDEX=loop_idx
              ENDDO !loop_idx
            ELSE
              DO loop_idx=1,NUMBER_OF_SUB_LOOPS
                CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR=>OLD_SUB_LOOPS(loop_idx)%PTR
              ENDDO !loop_idx
              DO loop_idx=NUMBER_OF_SUB_LOOPS+1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
                CALL ControlLoop_Finalise(OLD_SUB_LOOPS(loop_idx)%PTR,ERR,ERROR,*999)
              ENDDO !loop_idx
            ENDIF
            IF(ALLOCATED(OLD_SUB_LOOPS)) DEALLOCATE(OLD_SUB_LOOPS)
            CONTROL_LOOP%NUMBER_OF_SUB_LOOPS=NUMBER_OF_SUB_LOOPS
          ENDIF
        ELSE
          LOCAL_ERROR="The given number of sub loops of "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))// &
            & " is invalid. The number of sub loops must be >= 0."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF      
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET")
    RETURN
999 IF(ALLOCATED(OLD_SUB_LOOPS)) DEALLOCATE(OLD_SUB_LOOPS)
    ERRORSEXITS("CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET

  !
  !================================================================================================================================
  !

  !>Gets the output type for a control loop. \see OpenCMISS::cmfe_ControlLoop_OutputTypeGet
  SUBROUTINE CONTROL_LOOP_OUTPUT_TYPE_GET(CONTROL_LOOP,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to the control loop to get the output type for.
    INTEGER(INTG), INTENT(OUT) :: OUTPUT_TYPE !<On exit, the output type of the control loop \see CONTROL_LOOP_ROUTINES_OutputTypes,CONTROL_LOOP_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_OUTPUT_TYPE_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        OUTPUT_TYPE=CONTROL_LOOP%outputType
      ELSE
        CALL FlagError("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_OUTPUT_TYPE_GET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_OUTPUT_TYPE_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_OUTPUT_TYPE_GET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for a control loop. \see OpenCMISS::cmfe_ControlLoop_OutputTypeSet
  SUBROUTINE CONTROL_LOOP_OUTPUT_TYPE_SET(CONTROL_LOOP,OUTPUT_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer the control loop to set the output type for
    INTEGER(INTG), INTENT(IN) :: OUTPUT_TYPE !<The type of control loop output to be set \see CONTROL_LOOP_ROUTINES_OutputTypes,CONTROL_LOOP_ROUTINES
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CONTROL_LOOP_OUTPUT_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE        
        SELECT CASE(OUTPUT_TYPE)
        CASE(CONTROL_LOOP_NO_OUTPUT)
          CONTROL_LOOP%outputType=CONTROL_LOOP_NO_OUTPUT
        CASE(CONTROL_LOOP_PROGRESS_OUTPUT)
          CONTROL_LOOP%outputType=CONTROL_LOOP_PROGRESS_OUTPUT
        CASE(CONTROL_LOOP_TIMING_OUTPUT)
          CONTROL_LOOP%outputType=CONTROL_LOOP_TIMING_OUTPUT
        CASE DEFAULT
          LOCAL_ERROR="The specified control loop output type of "// &
            & TRIM(NUMBER_TO_VSTRING(OUTPUT_TYPE,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("CONTROL_LOOP_OUTPUT_TYPE_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_OUTPUT_TYPE_SET",ERR,ERROR)    
    RETURN 1
   
  END SUBROUTINE CONTROL_LOOP_OUTPUT_TYPE_SET
        
  !
  !================================================================================================================================
  !

  !>Updates the previous values for dependent variables under a time control loop
  SUBROUTINE ControlLoop_PreviousValuesUpdate(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer the time control loop to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearity,timeDependence,variableIdx
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("ControlLoop_PreviousValuesUpdate",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(controlLoop%fieldVariables)) CALL FlagError("Control loop field variables is not associated.",err,error,*999)

    DO variableIdx=1,controlLoop%fieldVariables%numberOfFieldVariables
      fieldVariable=>controlLoop%fieldVariables%fieldVariables(variableIdx)%fieldVariable
      linearity=controlLoop%fieldVariables%fieldVariables(variableIdx)%linearity
      timeDependence=controlLoop%fieldVariables%fieldVariables(variableIdx)%timeDependence
      IF(ASSOCIATED(fieldVariable)) THEN
        SELECT CASE(timeDependence)
        CASE(CONTROL_LOOP_FIELD_VARIABLE_STATIC)
          !Do nothing additional
        CASE(CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC)
          !Copy field variable values to the previous values
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
            & 1.0_DP,err,error,*999)
          SELECT CASE(linearity)
          CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
            !Do nothing additional
          CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
            !Copy residual to the previous residual
            CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
              & 1.0_DP,err,error,*999)
          CASE DEFAULT
            localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC)
          !Copy field variable values to the previous values
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
            & 1.0_DP,err,error,*999)
          !Copy velocity values to the previous velocity
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VELOCITY_VALUES_SET_TYPE, &
            & FIELD_PREVIOUS_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
          SELECT CASE(linearity)
          CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
            !Do nothing additional
          CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
            !Copy residual to the previous residual
            CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
              & 1.0_DP,err,error,*999)
          CASE DEFAULT
            localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC)
          !Copy field variable values to the previous values
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
            & 1.0_DP,err,error,*999)
          !Copy velocity values to the previous velocity
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VELOCITY_VALUES_SET_TYPE, &
            & FIELD_PREVIOUS_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
          !Copy acceleration values to the previous velocity
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
            & FIELD_PREVIOUS_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
          SELECT CASE(linearity)
          CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
            !Do nothing additional
          CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
            !Copy residual to the previous residual
            CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
              & 1.0_DP,err,error,*999)
          CASE DEFAULT
            localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
          & " is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDDO !variableIdx
    
    EXITS("ControlLoop_PreviousValuesUpdate")
    RETURN
999 ERRORSEXITS("ControlLoop_PreviousValuesUpdate",err,error)
    RETURN 1

  END SUBROUTINE ControlLoop_PreviousValuesUpdate

  !
  !================================================================================================================================
  !

  !>Finalises a simple control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_SIMPLE_FINALISE(SIMPLE_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_SIMPLE_TYPE), POINTER, INTENT(INOUT) :: SIMPLE_LOOP !<A pointer to the simple control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_SIMPLE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(SIMPLE_LOOP)) THEN
      DEALLOCATE(SIMPLE_LOOP)
    ENDIF
       
    EXITS("CONTROL_LOOP_SIMPLE_FINALISE")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_SIMPLE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SIMPLE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a simple loop for a control loop.
  SUBROUTINE CONTROL_LOOP_SIMPLE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to initialise the simple loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("CONTROL_LOOP_SIMPLE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%SIMPLE_LOOP)) THEN
        CALL FlagError("The simple loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%SIMPLE_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate simple loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%SIMPLE_LOOP%CONTROL_LOOP=>CONTROL_LOOP
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("CONTROL_LOOP_SIMPLE_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_SIMPLE_FINALISE(CONTROL_LOOP%SIMPLE_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("CONTROL_LOOP_SIMPLE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SIMPLE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solvers for a control loop and all sub control loops.
  RECURSIVE SUBROUTINE CONTROL_LOOP_SOLVERS_DESTROY(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to destroy the solvers for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2

    ENTERS("CONTROL_LOOP_SOLVERS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      !Destroy the solvers in any sub control loops first
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
          CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
          CALL CONTROL_LOOP_SOLVERS_DESTROY(CONTROL_LOOP2,ERR,ERROR,*999)
        ENDDO !loop_idx
      ENDIF
      !Destroy the solvers in this control loop
      IF(ASSOCIATED(CONTROL_LOOP%SOLVERS)) CALL SOLVERS_DESTROY(CONTROL_LOOP%SOLVERS,ERR,ERROR,*999)
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_SOLVERS_DESTROY")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_SOLVERS_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SOLVERS_DESTROY

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solver equations for a control loop and all sub control loops. \todo Create solvers_solver_equations_destory and call?
  RECURSIVE SUBROUTINE CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to destroy the solver for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: loop_idx,solver_idx
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP2
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
 
    ENTERS("CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      !Destroy the solver equations in any sub control loops first
      IF(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS>0) THEN
        DO loop_idx=1,CONTROL_LOOP%NUMBER_OF_SUB_LOOPS
          CONTROL_LOOP2=>CONTROL_LOOP%SUB_LOOPS(loop_idx)%PTR
          CALL CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY(CONTROL_LOOP2,ERR,ERROR,*999)
        ENDDO !loop_idx
      ENDIF
      !Destroy the solver equations in this control loop
      IF(ASSOCIATED(CONTROL_LOOP%SOLVERS)) THEN
        DO solver_idx=1,CONTROL_LOOP%SOLVERS%NUMBER_OF_SOLVERS
          SOLVER=>CONTROL_LOOP%SOLVERS%SOLVERS(solver_idx)%PTR
          IF(ASSOCIATED(SOLVER)) THEN
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) CALL SOLVER_EQUATIONS_DESTROY(SOLVER_EQUATIONS,ERR,ERROR,*999)
          ELSE
            CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
          ENDIF
        ENDDO !solver_idx
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SOLVER_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Gets/returns a pointer to the sub loops as specified by the sub loop index for a control loop.
  SUBROUTINE CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,SUB_LOOP_INDEX,SUB_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to get the sub loop for
    INTEGER(INTG), INTENT(IN) :: SUB_LOOP_INDEX !<The sub loop index to get
!    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(OUT) :: SUB_LOOP !<On return, a pointer to the specified sub loop. Must not be associated on entry.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: SUB_LOOP !<On return, a pointer to the specified sub loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
 
    ENTERS("CONTROL_LOOP_SUB_LOOP_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SUB_LOOP)) THEN
        CALL FlagError("Sub loop is already associated.",ERR,ERROR,*999)
      ELSE
        NULLIFY(SUB_LOOP)
        IF(SUB_LOOP_INDEX>0.AND.SUB_LOOP_INDEX<=CONTROL_LOOP%NUMBER_OF_SUB_LOOPS) THEN
          SUB_LOOP=>CONTROL_LOOP%SUB_LOOPS(SUB_LOOP_INDEX)%PTR
        ELSE
          LOCAL_ERROR="The specified sub loop index of "//TRIM(NUMBER_TO_VSTRING(SUB_LOOP_INDEX,"*",ERR,ERROR))// &
            & " is invalid. The sub loop index must be > 0 and <= "// &
            & TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%NUMBER_OF_SUB_LOOPS,"*",ERR,ERROR))//"."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        ENDIF
      ENDIF      
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_SUB_LOOP_GET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_SUB_LOOP_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_SUB_LOOP_GET

  !
  !================================================================================================================================
  !

  !>Finalises a time control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_TIME_FINALISE(TIME_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER, INTENT(INOUT) :: TIME_LOOP !<A pointer to the time control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_TIME_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(TIME_LOOP)) THEN
      DEALLOCATE(TIME_LOOP)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIME_FINALISE")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TIME_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIME_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a time loop for a control loop.
  SUBROUTINE CONTROL_LOOP_TIME_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to initialise the time loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("CONTROL_LOOP_TIME_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%TIME_LOOP)) THEN
        CALL FlagError("The time loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%TIME_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate time loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%TIME_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%TIME_LOOP%NUMBER_OF_ITERATIONS=0
        CONTROL_LOOP%TIME_LOOP%GLOBAL_ITERATION_NUMBER=0
        CONTROL_LOOP%TIME_LOOP%CURRENT_TIME=0.0_DP
        CONTROL_LOOP%TIME_LOOP%START_TIME=0.0_DP
        CONTROL_LOOP%TIME_LOOP%STOP_TIME=1.0_DP
        CONTROL_LOOP%TIME_LOOP%TIME_INCREMENT=0.01_DP
        CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER=0
        CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER=0
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIME_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_TIME_FINALISE(CONTROL_LOOP%TIME_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("CONTROL_LOOP_TIME_INITIALISE",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE CONTROL_LOOP_TIME_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Gets the current time parameters for a time control loop. \see OpenCMISS::cmfe_ControlLoop_CurrentTimesGet
  SUBROUTINE CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
    & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP
    REAL(DP), INTENT(OUT) :: START_TIME
    REAL(DP), INTENT(OUT) :: STOP_TIME
    REAL(DP), INTENT(OUT) :: CURRENT_TIME
    REAL(DP), INTENT(OUT) :: TIME_INCREMENT
    INTEGER(INTG), INTENT(OUT) :: CURRENT_LOOP_ITERATION
    INTEGER(INTG), INTENT(OUT) :: OUTPUT_ITERATION_NUMBER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables    
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(CONTROL_LOOP_TYPE), POINTER :: PARENT_LOOP
    INTEGER(INTG), POINTER :: CONTROL_LOOP_LEVEL
    INTEGER(INTG) :: I

    ENTERS("CONTROL_LOOP_TIMES_GET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CONTROL_LOOP_LEVEL=>CONTROL_LOOP%CONTROL_LOOP_LEVEL
        PARENT_LOOP=>CONTROL_LOOP
        DO I=CONTROL_LOOP_LEVEL,1,-1
          IF(CONTROL_LOOP_LEVEL==0) THEN
            CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
          ELSE
            IF(PARENT_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
              TIME_LOOP=>PARENT_LOOP%TIME_LOOP
              IF(ASSOCIATED(TIME_LOOP)) THEN
                START_TIME=TIME_LOOP%START_TIME
                STOP_TIME=TIME_LOOP%STOP_TIME
                CURRENT_TIME=TIME_LOOP%CURRENT_TIME
                TIME_INCREMENT=TIME_LOOP%TIME_INCREMENT
                CURRENT_LOOP_ITERATION=TIME_LOOP%ITERATION_NUMBER
                OUTPUT_ITERATION_NUMBER=TIME_LOOP%OUTPUT_NUMBER
              ELSE
                CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
              ENDIF
              EXIT
            ELSE
              PARENT_LOOP=>PARENT_LOOP%PARENT_LOOP
            ENDIF
          ENDIF
        ENDDO
      ELSE
        CALL FlagError("Control loop has not been finished.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIMES_GET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TIMES_GET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIMES_GET

  !
  !================================================================================================================================
  !

  !>Sets the time parameters for a time control loop. \see OpenCMISS::cmfe_ControlLoop_TimesSet
  SUBROUTINE CONTROL_LOOP_TIMES_SET(CONTROL_LOOP,START_TIME,STOP_TIME,TIME_INCREMENT,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to set the times for
    REAL(DP), INTENT(IN) :: START_TIME !<The start time for the time control loop.
    REAL(DP), INTENT(IN) :: STOP_TIME !<The stop time for the time control loop.
    REAL(DP), INTENT(IN) :: TIME_INCREMENT !<The time increment for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CONTROL_LOOP_TIMES_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
        IF(ASSOCIATED(TIME_LOOP)) THEN
          IF(ABS(TIME_INCREMENT)<=ZERO_TOLERANCE) THEN
            LOCAL_ERROR="The specified time increment of "//TRIM(NUMBER_TO_VSTRING(TIME_INCREMENT,"*",ERR,ERROR))// &
              & " is invalid. The time increment must not be zero."          
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          ELSE
            IF(TIME_INCREMENT>0.0_DP) THEN
              IF(STOP_TIME<=START_TIME) THEN
                LOCAL_ERROR="The specified stop time of "//TRIM(NUMBER_TO_VSTRING(STOP_TIME,"*",ERR,ERROR))// &
                  & " is incompatiable with a specified start time of "//TRIM(NUMBER_TO_VSTRING(START_TIME,"*",ERR,ERROR))// &
                  & ". For a positive time increment the stop time must be > than the start time."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ELSE
              IF(START_TIME<=STOP_TIME) THEN
                LOCAL_ERROR="The specified start time of "//TRIM(NUMBER_TO_VSTRING(START_TIME,"*",ERR,ERROR))// &
                  & " is incompatiable with a specified stop time of "//TRIM(NUMBER_TO_VSTRING(STOP_TIME,"*",ERR,ERROR))// &
                  & ". For a negative time increment the stop time must be < than the start time."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              ENDIF
            ENDIF
          ENDIF
          TIME_LOOP%START_TIME=START_TIME
          TIME_LOOP%STOP_TIME=STOP_TIME
          TIME_LOOP%TIME_INCREMENT=TIME_INCREMENT
          TIME_LOOP%NUMBER_OF_ITERATIONS=0    ! reset number of iterations
        ELSE
          CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIMES_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TIMES_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIMES_SET

  !
  !================================================================================================================================
  !

  !>Sets the output parameters for a time control loop. \see OpenCMISS::cmfe_ControlLoop_TimeOutputSet
  SUBROUTINE CONTROL_LOOP_TIME_OUTPUT_SET(CONTROL_LOOP,OUTPUT_FREQUENCY,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to set the times for
    INTEGER(INTG) :: OUTPUT_FREQUENCY !<The output frequency modulo to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    
    ENTERS("CONTROL_LOOP_TIME_OUTPUT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            IF(OUTPUT_FREQUENCY>=0) THEN
              TIME_LOOP%OUTPUT_NUMBER=OUTPUT_FREQUENCY
            ELSE
              CALL FlagError("Invalid output frequency. The frequency should be greater than or equal to zero, but is "// &
                & TRIM(NUMBER_TO_VSTRING(OUTPUT_FREQUENCY,"*",ERR,ERROR))//".",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIME_OUTPUT_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TIME_OUTPUT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIME_OUTPUT_SET

  !
  !================================================================================================================================
  !

  !>Sets the input parameters for a time control loop. \see OpenCMISS::cmfe_ControlLoop_TimeInputSet
  SUBROUTINE CONTROL_LOOP_TIME_INPUT_SET(CONTROL_LOOP,INPUT_OPTION,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: CONTROL_LOOP !<A pointer to control loop to set the times for
    INTEGER(INTG) :: INPUT_OPTION !<The input option modulo to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables    
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP
    
    ENTERS("CONTROL_LOOP_TIME_INPUT_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has been finished.",ERR,ERROR,*999)
      ELSE
        IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            TIME_LOOP%INPUT_NUMBER=INPUT_OPTION
          ELSE
            CALL FlagError("Control loop time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        ELSE
          CALL FlagError("The specified control loop is not a time control loop.",ERR,ERROR,*999)
        ENDIF
      ENDIF          
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_TIME_INPUT_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TIME_INPUT_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TIME_INPUT_SET
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the control loop type.
  SUBROUTINE CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,LOOP_TYPE,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to control loop to set the type of
    INTEGER(INTG), INTENT(IN) :: LOOP_TYPE !<The type of loop type to set \see PROBLEM_CONSTANTS_ControlLoopTypes,PROBLEM_CONSTANTS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("CONTROL_LOOP_TYPE_SET",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%CONTROL_LOOP_FINISHED) THEN
        CALL FlagError("Control loop has already been finished.",ERR,ERROR,*999)
      ELSE
        IF(LOOP_TYPE/=CONTROL_LOOP%LOOP_TYPE) THEN
          !Initialise the new loop type
          SELECT CASE(LOOP_TYPE)
          CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
            CALL CONTROL_LOOP_SIMPLE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
            CALL CONTROL_LOOP_FIXED_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_TIME_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_WHILE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL CONTROL_LOOP_LOAD_INCREMENT_INITIALISE(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The loop type of "//TRIM(NUMBER_TO_VSTRING(LOOP_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          !Finialise the old loop type
          SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
          CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
            CALL CONTROL_LOOP_SIMPLE_FINALISE(CONTROL_LOOP%SIMPLE_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
            CALL CONTROL_LOOP_FIXED_FINALISE(CONTROL_LOOP%FIXED_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
            CALL CONTROL_LOOP_TIME_FINALISE(CONTROL_LOOP%TIME_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
            CALL CONTROL_LOOP_WHILE_FINALISE(CONTROL_LOOP%WHILE_LOOP,ERR,ERROR,*999)
          CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
            CALL CONTROL_LOOP_LOAD_INCREMENT_FINALISE(CONTROL_LOOP%LOAD_INCREMENT_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The control loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%LOOP_TYPE,"*",ERR,ERROR))//" is invalid."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CONTROL_LOOP%LOOP_TYPE=LOOP_TYPE
        ENDIF
      ENDIF      
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("CONTROL_LOOP_TYPE_SET")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_TYPE_SET",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_TYPE_SET

  !
  !================================================================================================================================
  !

  !>Finalises a while control loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_WHILE_FINALISE(WHILE_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER, INTENT(INOUT) :: WHILE_LOOP !<A pointer to the while control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_WHILE_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(WHILE_LOOP)) THEN
      DEALLOCATE(WHILE_LOOP)
    ENDIF
       
    EXITS("CONTROL_LOOP_WHILE_FINALISE")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_WHILE_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_WHILE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a while loop for a control loop.
  SUBROUTINE CONTROL_LOOP_WHILE_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to initialise the while loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("CONTROL_LOOP_WHILE_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%WHILE_LOOP)) THEN
        CALL FlagError("The while loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%WHILE_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate while loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%WHILE_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%WHILE_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=100
        CONTROL_LOOP%WHILE_LOOP%ABSOLUTE_TOLERANCE=1.0E-5_DP
        CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("CONTROL_LOOP_WHILE_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_WHILE_FINALISE(CONTROL_LOOP%WHILE_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("CONTROL_LOOP_WHILE_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_WHILE_INITIALISE

  !
  !================================================================================================================================
  !
  
  !>Finalises a load increment loop and deallocates all memory.
  SUBROUTINE CONTROL_LOOP_LOAD_INCREMENT_FINALISE(LOAD_INCREMENT_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_LOAD_INCREMENT_TYPE), POINTER, INTENT(INOUT) :: LOAD_INCREMENT_LOOP !<A pointer to the load increment control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
 
    ENTERS("CONTROL_LOOP_LOAD_INCREMENT_FINALISE",ERR,ERROR,*999)

    IF(ASSOCIATED(LOAD_INCREMENT_LOOP)) THEN
      DEALLOCATE(LOAD_INCREMENT_LOOP)
    ENDIF
       
    EXITS("CONTROL_LOOP_LOAD_INCREMENT_FINALISE")
    RETURN
999 ERRORSEXITS("CONTROL_LOOP_LOAD_INCREMENT_FINALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_LOAD_INCREMENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises a load increment loop for a control loop.
  SUBROUTINE CONTROL_LOOP_LOAD_INCREMENT_INITIALISE(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(INOUT) :: CONTROL_LOOP !<A pointer to the control loop to initialise the load increment loop for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: DUMMY_ERR
    TYPE(VARYING_STRING) :: DUMMY_ERROR
 
    ENTERS("CONTROL_LOOP_LOAD_INCREMENT_INITIALISE",ERR,ERROR,*998)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%LOAD_INCREMENT_LOOP)) THEN
        CALL FlagError("The load increment loop is already associated for this control loop.",ERR,ERROR,*998)
      ELSE
        ALLOCATE(CONTROL_LOOP%LOAD_INCREMENT_LOOP,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate load increment loop for the control loop.",ERR,ERROR,*999)
        CONTROL_LOOP%LOAD_INCREMENT_LOOP%CONTROL_LOOP=>CONTROL_LOOP
        CONTROL_LOOP%LOAD_INCREMENT_LOOP%ITERATION_NUMBER=0
        CONTROL_LOOP%LOAD_INCREMENT_LOOP%MAXIMUM_NUMBER_OF_ITERATIONS=1 ! default is full load in one step
        CONTROL_LOOP%LOAD_INCREMENT_LOOP%OUTPUT_NUMBER=0
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*998)
    ENDIF
       
    EXITS("CONTROL_LOOP_LOAD_INCREMENT_INITIALISE")
    RETURN
999 CALL CONTROL_LOOP_LOAD_INCREMENT_FINALISE(CONTROL_LOOP%LOAD_INCREMENT_LOOP,DUMMY_ERR,DUMMY_ERROR,*998)
998 ERRORSEXITS("CONTROL_LOOP_LOAD_INCREMENT_INITIALISE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE CONTROL_LOOP_LOAD_INCREMENT_INITIALISE

  !
  !================================================================================================================================
  !

END MODULE CONTROL_LOOP_ROUTINES

