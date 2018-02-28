!> \file
!> \author Chris Bradley
!> \brief This module contains all the low-level base routines e.g., all debug, control, and low-level communication routines.
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

!> This module contains all the low-level base routines e.g., all debug, control, and low-level communication routines.
MODULE BaseRoutines

  USE Constants
  USE Kinds
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING
  USE MACHINE_CONSTANTS 

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  INTEGER(INTG), PARAMETER :: MAX_OUTPUT_LINES=500 !<Maximum number of lines that can be output \see BaseRoutines::WriteStr
  INTEGER(INTG), PARAMETER :: MAX_OUTPUT_WIDTH=132 !<Maximum width of output line \see BaseRoutines::WriteStr

  !> \addtogroup BaseRoutines_OutputType BaseRoutines::OutputType
  !> \brief Output type parameter
  !> \see BaseRoutines
  !>@{  
  INTEGER(INTG), PARAMETER :: GENERAL_OUTPUT_TYPE=1 !<General output type \see BaseRoutines_OutputType,BaseRoutines
  INTEGER(INTG), PARAMETER :: DIAGNOSTIC_OUTPUT_TYPE=2 !<Diagnostic output type \see BaseRoutines_OutputType,BaseRoutines
  INTEGER(INTG), PARAMETER :: TIMING_OUTPUT_TYPE=3 !<Timing output type \see BaseRoutines_OutputType,BaseRoutines
  INTEGER(INTG), PARAMETER :: ERROR_OUTPUT_TYPE=4 !<Error output type \see BaseRoutines_OutputType,BaseRoutines
  INTEGER(INTG), PARAMETER :: WARNING_OUTPUT_TYPE=5 !<Warning output type \see BaseRoutines_OutputType,BaseRoutines
  INTEGER(INTG), PARAMETER :: HELP_OUTPUT_TYPE=6 !<Help output type \see BaseRoutines_OutputType,BaseRoutines
  !>@}

  !> \addtogroup BaseRoutines_FileUnits BaseRoutines::FileUnits
  !> \brief File unit parameters
  !> \see BaseRoutines
  !>@{  
  INTEGER(INTG), PARAMETER :: ECHO_FILE_UNIT=10 !<File unit for echo files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: DIAGNOSTICS_FILE_UNIT=11 !<File unit for diagnostic files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: TIMING_FILE_UNIT=12 !<File unit for timing files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: LEARN_FILE_UNIT=13 !<File unit for learn files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: IO1_FILE_UNIT=21 !<File unit for general IO 1 files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: IO2_FILE_UNIT=22 !<File unit for general IO 2 files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: IO3_FILE_UNIT=23 !<File unit for general IO 3 files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: IO4_FILE_UNIT=24 !<File unit for general IO 4 files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: IO5_FILE_UNIT=25 !<File unit for general IO 5 files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: TEMPORARY_FILE_UNIT=80 !<File unit for temporary files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: OPEN_COMFILE_UNIT=90 !<File unit for open command files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: START_READ_COMFILE_UNIT=90 !<First file unit for read command files \see BaseRoutines_FileUnits,BaseRoutines
  INTEGER(INTG), PARAMETER :: STOP_READ_COMFILE_UNIT=99 !<Last file unit for read command files \see BaseRoutines_FileUnits,BaseRoutines
  !>@}

  !> \addtogroup BaseRoutines_DiagnosticTypes BaseRoutines::DiagnosticTypes
  !> \brief Diganostic type parameters
  !> \see BaseRoutines,OpenCMISS_DiagnosticTypes
  !>@{  
  INTEGER(INTG), PARAMETER :: ALL_DIAG_TYPE=1 !<Type for setting diagnostic output in all routines \see BaseRoutines_DiagnosticTypes,BaseRoutines
  INTEGER(INTG), PARAMETER :: IN_DIAG_TYPE=2 !<Type for setting diagnostic output in one routine \see BaseRoutines_DiagnosticTypes,BaseRoutines
  INTEGER(INTG), PARAMETER :: FROM_DIAG_TYPE=3 !<Type for setting diagnostic output from one routine downwards \see BaseRoutines_DiagnosticTypes,BaseRoutines
  !>@}

  !> \addtogroup BaseRoutines_TimingTypes BaseRoutines::TimingTypes
  !> \brief Timing type parameters
  !> \see BaseRoutines,OpenCMISS_TimingTypes
  !>@{  
  INTEGER(INTG), PARAMETER :: ALL_TIMING_TYPE=1 !<Type for setting timing output in all routines \see BaseRoutines_TimingTypes,BaseRoutines
  INTEGER(INTG), PARAMETER :: IN_TIMING_TYPE=2 !<Type for setting timing output in one routine \see BaseRoutines_TimingTypes,BaseRoutines
  INTEGER(INTG), PARAMETER :: FROM_TIMING_TYPE=3 !<Type for setting timing output from one routine downwards \see BaseRoutines_TimingTypes,BaseRoutines
  !>@}

  !Module types

  !>Contains information for an item in the routine list for diagnostics or timing
  TYPE RoutineListItemType
    CHARACTER(LEN=63) :: name !<Name of the routine
    INTEGER(INTG) :: numberOfInvocations !<Number of times the routine has been invocted
    REAL(SP) :: totalInclusiveCPUTime !<Total User CPU time spent in the routine inclusive of calls
    REAL(SP) :: totalInclusiveSystemTime !<Total System CPU time spent in the routine inclusive of calls
    REAL(SP) :: totalExclusiveCPUTime !<Total User CPU time spent in the routine exclusive of calls
    REAL(SP) :: totalExclusiveSystemTime !<Total System CPU time spent in the routine exclusive of calls
    TYPE(RoutineListItemType), POINTER :: nextRoutine !<Pointer to the next routine item in the routine list
  END TYPE RoutineListItemType

  !>Contains information for the routine list for diagnostics or timing
  TYPE RoutineListType
    TYPE(RoutineListItemType), POINTER :: head !<A pointer to the head of the routine list.
  END TYPE RoutineListType

  !>Contains information for an item in the routine invocation stack
  TYPE RoutineStackItemType
    CHARACTER(LEN=63) :: name !<Name of the routine
    REAL(SP) :: inclusiveCPUTime !<User CPU time spent in the routine inclusive of calls 
    REAL(SP) :: inclusiveSystemTime !<System CPU time spent in the routine inclusive of calls 
    REAL(SP) :: exclusiveCPUTime !<User CPU time spent in the routine exclusive of calls 
    REAL(SP) :: exclusiveSystemTime !<System CPU time spent in the routine exclusive of calls 
    LOGICAL :: diagnostics !<.TRUE. if diagnostics are active in the routine
    LOGICAL :: timing !<.TRUE. if timing is active in the routine
    TYPE(RoutineListItemType), POINTER :: routineListItem !<Pointer to the routine list item for diagnostics or timing
    TYPE(RoutineStackItemType), POINTER :: previousRoutine !<Pointer to the previous routine in the routine stack
  END TYPE RoutineStackItemType

  !>Contains information for the routine invocation stack
  TYPE RoutineStackType
    TYPE(RoutineStackItemType), POINTER :: stackPointer !<Pointer to the top of the stack
  END TYPE RoutineStackType

  !Module variables

  INTEGER(INTG), SAVE :: myComputationalNodeNumber !<The computational rank for this node
  INTEGER(INTG), SAVE :: numberOfComputationalNodes !<The number of computational nodes
  INTEGER(INTG), ALLOCATABLE :: cmissRandomSeeds(:) !<The current error handling seeds for OpenCMISS
  LOGICAL, SAVE :: diagnostics !<.TRUE. if diagnostic output is required in any routines.
  LOGICAL, SAVE :: diagnostics1 !<.TRUE. if level 1 diagnostic output is active in the current routine
  LOGICAL, SAVE :: diagnostics2 !<.TRUE. if level 2 diagnostic output is active in the current routine
  LOGICAL, SAVE :: diagnostics3 !<.TRUE. if level 3 diagnostic output is active in the current routine
  LOGICAL, SAVE :: diagnostics4 !<.TRUE. if level 4 diagnostic output is active in the current routine
  LOGICAL, SAVE :: diagnostics5 !<.TRUE. if level 5 diagnostic output is active in the current routine
  LOGICAL, SAVE :: diagnosticsLevel1 !<.TRUE. if the user has requested level 1 diagnostic output to be active
  LOGICAL, SAVE :: diagnosticsLevel2 !<.TRUE. if the user has requested level 2 diagnostic output to be active
  LOGICAL, SAVE :: diagnosticsLevel3 !<.TRUE. if the user has requested level 3 diagnostic output to be active
  LOGICAL, SAVE :: diagnosticsLevel4 !<.TRUE. if the user has requested level 4 diagnostic output to be active
  LOGICAL, SAVE :: diagnosticsLevel5 !<.TRUE. if the user has requested level 5 diagnostic output to be active
  LOGICAL, SAVE :: diagAllSubroutines !<.TRUE. if diagnostic output is required in all routines
  LOGICAL, SAVE :: diagFromSubroutine !<.TRUE. if diagnostic output is required from a particular routine
  LOGICAL, SAVE :: diagFileOpen !<.TRUE. if the diagnostic output file is open
  LOGICAL, SAVE :: diagOrTiming !<.TRUE. if diagnostics or time is .TRUE.
  LOGICAL, SAVE :: echoOutput !<.TRUE. if all output is to be echoed to the echo file
  LOGICAL, SAVE :: timing !<.TRUE. if timing output is required in any routines.
  LOGICAL, SAVE :: timingSummary !<.TRUE. if timing output will be summary form via a TimingSummaryOutput call otherwise timing will be output for routines when the routine exits \see BaseRoutines::TimingSummaryOutput
  LOGICAL, SAVE :: timingAllSubroutines !<.TRUE. if timing output is required in all routines
  LOGICAL, SAVE :: timingFromSubroutine !<.TRUE. if timing output is required from a particular routine
  LOGICAL, SAVE :: timingFileOpen !<.TRUE. if the timing output file is open
  CHARACTER(LEN=MAXSTRLEN), SAVE :: outputString(MAX_OUTPUT_LINES) !<The array of lines to output
  TYPE(RoutineListType), SAVE :: diagRoutineList !<The list of routines for which diagnostic output is required
  TYPE(RoutineListType), SAVE :: timingRoutineList !<The list of routines for which timing output is required
  TYPE(RoutineStackType), SAVE :: routineStack !<The routime invocation stack

  !Interfaces

  INTERFACE

!!!!NOTE: This module needs to call the c cputime function directly in order to avoid a circular module loop when timer uses
!!!!      base_routines.

    SUBROUTINE CPUTimer(returnTime, timeType, err, cError) BIND(C,NAME="CPUTimer")
      USE ISO_C_BINDING
      REAL(C_DOUBLE), INTENT(OUT) :: returnTime
      INTEGER(C_INT), INTENT(IN) :: timeType
      INTEGER(C_INT), INTENT(OUT) :: err
      CHARACTER(C_CHAR), INTENT(OUT) :: cError(*)
    END SUBROUTINE CPUTimer

  END INTERFACE

  !>Extracts the error message part of an error string
  INTERFACE ExtractErrorMessage
    MODULE PROCEDURE ExtractErrorMessageC
    MODULE PROCEDURE ExtractErrorMessageVS
  END INTERFACE ExtractErrorMessage

  !>Extracts the error stack part of the error string
  INTERFACE ExtractErrorStack
    MODULE PROCEDURE ExtractErrorStackC
    MODULE PROCEDURE ExtractErrorStackVS
  END INTERFACE ExtractErrorStack

  !>Flags an error condition 
  INTERFACE FLAG_ERROR
    MODULE PROCEDURE FlagErrorC
    MODULE PROCEDURE FlagErrorVS
  END INTERFACE FLAG_ERROR
  
  !>Flags an error condition 
  INTERFACE FlagError
    MODULE PROCEDURE FlagErrorC
    MODULE PROCEDURE FlagErrorVS
  END INTERFACE FlagError
  
  !>Flags a warning to the user
  INTERFACE FLAG_WARNING
    MODULE PROCEDURE FlagWarningC
    MODULE PROCEDURE FlagWarningVS
  END INTERFACE FLAG_WARNING

  !>Flags a warning to the user
  INTERFACE FlagWarning
    MODULE PROCEDURE FlagWarningC
    MODULE PROCEDURE FlagWarningVS
  END INTERFACE FlagWarning

  PUBLIC GENERAL_OUTPUT_TYPE,DIAGNOSTIC_OUTPUT_TYPE,TIMING_OUTPUT_TYPE,ERROR_OUTPUT_TYPE,HELP_OUTPUT_TYPE

  PUBLIC ALL_DIAG_TYPE,IN_DIAG_TYPE,FROM_DIAG_TYPE

  PUBLIC OPEN_COMFILE_UNIT,START_READ_COMFILE_UNIT,STOP_READ_COMFILE_UNIT,TEMPORARY_FILE_UNIT

  PUBLIC ALL_TIMING_TYPE,IN_TIMING_TYPE,FROM_TIMING_TYPE

  PUBLIC LEARN_FILE_UNIT,IO1_FILE_UNIT,IO2_FILE_UNIT,IO3_FILE_UNIT,IO4_FILE_UNIT,IO5_FILE_UNIT

  PUBLIC diagnostics1,diagnostics2,diagnostics3,diagnostics4,diagnostics5

  PUBLIC cmissRandomSeeds
  
  PUBLIC outputString

  PUBLIC BaseRoutines_Finalise,BaseRoutines_Initialise

  PUBLIC ComputationalNodeNumbersSet

  PUBLIC DiagnosticsSetOn,DiagnosticsSetOff

  PUBLIC Enters,Errors,Exits

  PUBLIC ExtractErrorMessage
  
  PUBLIC ExtractErrorStack
  
  PUBLIC FLAG_ERROR,FLAG_WARNING

  PUBLIC FlagError,FlagWarning
  
  PUBLIC OutputSetOff,OutputSetOn

  PUBLIC RandomSeedsGet,RandomSeedsSizeGet,RandomSeedsSet
  
  PUBLIC TimingSetOn,TimingSetOff
  
  PUBLIC TimingSummaryOutput
   
  PUBLIC WriteError
  
  PUBLIC WriteStr

CONTAINS

  !
  !================================================================================================================================
  !

  !>Records the entry into the named procedure and initialises the error code \see BaseRoutines::EXITS
  SUBROUTINE Enters(name,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: name !<The name of the routine being entered
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    CHARACTER(C_CHAR) :: cError(MAXSTRLEN)
    REAL(DP) :: entersCPUTime,entersSystemTime
    LOGICAL :: finished
    TYPE(RoutineListItemType), POINTER :: listRoutinePtr
    TYPE(RoutineStackItemType), POINTER :: newRoutinePtr,routinePtr

    IF(diagOrTiming) THEN
      !$OMP CRITICAL(ENTERS_1)
      ALLOCATE(newRoutinePtr,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new routine stack item.",err,error,*999)
      newRoutinePtr%diagnostics=.FALSE.
      newRoutinePtr%timing=.FALSE.
      newRoutinePtr%name=name(1:LEN_TRIM(name))
      IF(ASSOCIATED(routineStack%stackPointer)) THEN
        newRoutinePtr%previousRoutine=>routineStack%stackPointer
        routineStack%stackPointer=>newRoutinePtr
      ELSE
        NULLIFY(newRoutinePtr%previousRoutine)
        routineStack%stackPointer=>newRoutinePtr
      ENDIF
      routinePtr=>routineStack%stackPointer
      NULLIFY(routinePtr%routineListItem)
      IF(diagnostics) THEN
        IF(diagAllSubroutines) THEN !turn diagnostics on in all subroutines
          routinePtr%diagnostics=.TRUE.
        ELSE !diagnostics on in selected subroutines
          finished=.FALSE.
          listRoutinePtr=>diagRoutineList%head
          DO WHILE(ASSOCIATED(listRoutinePtr).AND..NOT.finished)
            IF(listRoutinePtr%name(1:LEN_TRIM(listRoutinePtr%name))==routinePtr%name(1:LEN_TRIM(routinePtr%name))) THEN
              routinePtr%diagnostics=.TRUE.
              routinePtr%routineListItem=>listRoutinePtr
              finished=.TRUE.
            ELSE
              listRoutinePtr=>listRoutinePtr%nextRoutine
            ENDIF
          ENDDO
          IF(diagFromSubroutine) THEN
            IF(ASSOCIATED(routinePtr%previousRoutine)) THEN
              IF(routinePtr%previousRoutine%diagnostics) routinePtr%diagnostics=.TRUE.
            ENDIF
          ENDIF
        ENDIF
        IF(routinePtr%diagnostics) THEN
          diagnostics1=diagnosticsLevel1
          diagnostics2=diagnosticsLevel2
          diagnostics3=diagnosticsLevel3
          diagnostics4=diagnosticsLevel4
          diagnostics5=diagnosticsLevel5
        ELSE
          diagnostics1=.FALSE.
          diagnostics2=.FALSE.
          diagnostics3=.FALSE.
          diagnostics4=.FALSE.
          diagnostics5=.FALSE.
        ENDIF
        IF(routinePtr%diagnostics) THEN
          WRITE(outputString,'("*** Enters: ",A)') name(1:LEN_TRIM(name))
          CALL WriteStr(DIAGNOSTIC_OUTPUT_TYPE,err,error,*999)
        ELSE IF(ASSOCIATED(routinePtr%previousRoutine)) THEN
          !CPB 16/05/2007 Only show the calls if we have level 3 diagnostics or higher
          IF(diagnostics3) THEN
            IF(routinePtr%previousRoutine%diagnostics) THEN
              WRITE(outputString,'("*** Calls : ",A)') name(1:LEN_TRIM(name))
              CALL WriteStr(DIAGNOSTIC_OUTPUT_TYPE,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(timing) THEN
        CALL CPUTimer(entersCPUTime,1,ERR,cError)
        CALL CPUTimer(entersSystemTime,2,ERR,cError)
        routinePtr%inclusiveCPUTime=REAL(entersCPUTime,SP)
        routinePtr%inclusiveSystemTime=REAL(entersSystemTime,SP)
        routinePtr%exclusiveCPUTime=0.0_SP
        routinePtr%exclusiveSystemTime=0.0_SP
        IF(timingAllSubroutines) THEN
          routinePtr%timing=.TRUE.
        ELSE
          finished=.FALSE.
          listRoutinePtr=>timingRoutineList%head
          DO WHILE(ASSOCIATED(listRoutinePtr).AND..NOT.finished)
            IF(listRoutinePtr%name(1:LEN_TRIM(listRoutinePtr%name))==routinePtr%name(1:LEN_TRIM(routinePtr%name))) THEN
              routinePtr%timing=.TRUE.
              routinePtr%routineListItem=>listRoutinePtr
              finished=.TRUE.
            ELSE
              listRoutinePtr=>listRoutinePtr%nextRoutine
            ENDIF
          ENDDO
          IF(timingFromSubroutine) THEN
            IF(ASSOCIATED(routinePtr%previousRoutine)) THEN
              IF(routinePtr%previousRoutine%timing) routinePtr%timing=.TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      !$OMP END CRITICAL(ENTERS_1)
    ENDIF

    RETURN
999 RETURN 1
    
  END SUBROUTINE Enters

  !
  !================================================================================================================================
  !

  !>Records the exiting error of the subroutine 
  SUBROUTINE Errors(name,err,error)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: name !<The name of the routine with an error condition
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
 
    IF(err==0) err=1
    !CPB 20/02/07 aix compiler does not like varying strings so split the concatenate statement up into two statements
    localError=error//ERROR_SEPARATOR_CONSTANT
    error=localError//name(1:LEN_TRIM(name))

    RETURN

  END SUBROUTINE Errors

  !
  !================================================================================================================================
  !

  !>Records the exit out of the named procedure \see BaseRoutines::ENTERS
  SUBROUTINE Exits(name)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: name !<The name of the routine exiting
    !Local variables
    CHARACTER(C_CHAR) :: cError(MAXSTRLEN)
    INTEGER(INTG) :: err
    REAL(DP) :: exitsCPUTime,exitsSystemTime
    TYPE(VARYING_STRING) :: error
    TYPE(RoutineStackItemType), POINTER :: previousRoutinePtr,routinePtr

    IF(diagOrTiming) THEN
      !$OMP CRITICAL(EXITS_1)
      routinePtr=>routineStack%stackPointer
      IF(ASSOCIATED(routinePtr)) THEN
        previousRoutinePtr=>routinePtr%previousRoutine
        IF(diagnostics) THEN
          IF(routinePtr%diagnostics) THEN
            WRITE(outputString,'("*** Exits : ",A)') name(1:LEN_TRIM(name))
            CALL WriteStr(DIAGNOSTIC_OUTPUT_TYPE,err,error,*999)
          ENDIF
          IF(ASSOCIATED(previousRoutinePtr)) THEN
            IF(previousRoutinePtr%diagnostics) THEN
              diagnostics1=diagnosticsLevel1
              diagnostics2=diagnosticsLevel2
              diagnostics3=diagnosticsLevel3
              diagnostics4=diagnosticsLevel4
              diagnostics5=diagnosticsLevel5
            ELSE
              diagnostics1=.FALSE.
              diagnostics2=.FALSE.
              diagnostics3=.FALSE.
              diagnostics4=.FALSE.
              diagnostics5=.FALSE.
            ENDIF
          ENDIF
        ENDIF

        IF(timing) THEN
          CALL CPUTimer(exitsCPUTime,1,ERR,cError)
          CALL CPUTimer(exitsSystemTime,2,ERR,cError)
          routinePtr%inclusiveCPUTime=ABS(REAL(exitsCPUTime,SP)-routinePtr%inclusiveCPUTime)
          routinePtr%inclusiveSystemTime=ABS(REAL(exitsSystemTime,SP)-routinePtr%inclusiveSystemTime)
          IF(ASSOCIATED(previousRoutinePtr)) THEN
            previousRoutinePtr%exclusiveCPUTime=previousRoutinePtr%exclusiveCPUTime+routinePtr%inclusiveCPUTime
            previousRoutinePtr%exclusiveSystemTime=previousRoutinePtr%exclusiveSystemTime+routinePtr%inclusiveSystemTime
          ENDIF
          IF(ASSOCIATED(routinePtr%routineListItem)) THEN
            routinePtr%routineListItem%numberOfInvocations=routinePtr%routineListItem%numberOfInvocations+1
            routinePtr%routineListItem%totalInclusiveCPUTime=routinePtr%routineListItem%totalInclusiveCPUTime+ &
              & routinePtr%inclusiveCPUTime
            routinePtr%routineListItem%totalInclusiveSystemTime=routinePtr%routineListItem%totalInclusiveSystemTime+ &
              & routinePtr%inclusiveSystemTime
            IF(ASSOCIATED(previousRoutinePtr)) THEN
              IF(ASSOCIATED(previousRoutinePtr%routineListItem)) THEN
                previousRoutinePtr%routineListItem%totalExclusiveCPUTime=previousRoutinePtr%routineListItem% &
                  & totalExclusiveCPUTime+previousRoutinePtr%exclusiveCPUTime
                previousRoutinePtr%routineListItem%totalExclusiveSystemTime=previousRoutinePtr%routineListItem% &
                  & totalExclusiveSystemTime+previousRoutinePtr%exclusiveSystemTime
              ENDIF
            ENDIF
          ENDIF
          IF(routinePtr%timing) THEN
            IF(.NOT.timingSummary) THEN
              WRITE(outputString,'("*** Timing : ",A)') name(1:LEN_TRIM(name))
              CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
              IF(ASSOCIATED(routinePtr%routineListItem)) THEN
                WRITE(outputString,'("***    Number of invocations: ",I10)') routinePtr%routineListItem%numberOfInvocations
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
                WRITE(outputString, &
                  & '("***    Routine times:  Call Inclusive   Call Exclusive   Total Inclusive   Average Inclusive")')
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
                WRITE(outputString,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"   ",E15.6,"   ",E17.6)')  &
                  & routinePtr%inclusiveCPUTime,routinePtr%inclusiveCPUTime-routinePtr%exclusiveCPUTime, &
                  & routinePtr%routineListItem%totalInclusiveCPUTime,routinePtr%routineListItem% &
                  & totalInclusiveCPUTime/REAL(routinePtr%routineListItem%numberOfInvocations,SP)
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
                WRITE(outputString,'("***    System    (s):  ",E14.6,"   ",E14.6,"   ",E15.6,"   ",E17.6)')  &
                  & routinePtr%inclusiveSystemTime,routinePtr%inclusiveSystemTime-routinePtr%exclusiveSystemTime, &
                  & routinePtr%routineListItem%totalInclusiveSystemTime,routinePtr%routineListItem% &
                  & totalInclusiveSystemTime/REAL(routinePtr%routineListItem%numberOfInvocations,SP)
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
              ELSE
                WRITE(outputString,'("***    Routine times:  Call Inclusive   Call Exclusive")')
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
                WRITE(outputString,'("***    CPU       (s):  ",E14.6,"   ",E14.6)')  &
                  & routinePtr%inclusiveCPUTime,routinePtr%inclusiveCPUTime-routinePtr%exclusiveCPUTime
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
                WRITE(outputString,'("***    System    (s):  ",E14.6,"   ",E14.6)')  &
                  & routinePtr%inclusiveSystemTime,routinePtr%inclusiveSystemTime-routinePtr%exclusiveSystemTime
                CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        IF(ASSOCIATED(previousRoutinePtr)) THEN
          routineStack%stackPointer=>previousRoutinePtr
        ELSE
          NULLIFY(routineStack%stackPointer)
        ENDIF

        !Delete the routine pointer
        DEALLOCATE(routinePtr)

        !ELSE ERROR????
      ENDIF
      !$OMP END CRITICAL(EXITS_1)
    ENDIF

999 RETURN
    
  END SUBROUTINE Exits

  !
  !================================================================================================================================
  !

#include "macros.h"

  !>Set the computational node numbers. Note: this is done as a subroutine as ComputationalEnvironment depends on BaseRoutines.
  SUBROUTINE ComputationalNodeNumbersSet(myNodeNumber,numberOfNodes,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: myNodeNumber !<The node number for this rank.
    INTEGER(INTG), INTENT(IN) :: numberOfNodes !<The number of computational nodes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("ComputationalNodeNumbersSet",err,error,*999)

    IF(numberOfNodes>0) THEN
      IF(myNodeNumber>=0.AND.myNodeNumber<=numberOfNodes-1) THEN        
        myComputationalNodeNumber=myNodeNumber
        numberOfComputationalNodes=numberOfNodes        
      ELSE
        CALL FlagError("Invalid node number.",err,error,*999)
      ENDIF
    ELSE
       CALL FlagError("Invalid number of nodes.",err,error,*999)
    ENDIF
    
    EXITS("ComputationalNodeNumbersSet")
    RETURN 
999 ERRORSEXITS("ComputationalNodeNumbersSet",err,error)
    RETURN 1
    
  END SUBROUTINE ComputationalNodeNumbersSet

  !
  !================================================================================================================================
  !

  !>Extracts the error message from a CMISS error string and returns it as a varying string
  SUBROUTINE ExtractErrorMessageVS(errorMessage,err,error,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(OUT) :: errorMessage !<The extracted error message
    INTEGER(INTG), INTENT(IN) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(IN) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: position

    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    errorMessage=EXTRACT(error,1,position-1)

    RETURN
    
  END SUBROUTINE ExtractErrorMessageVS

  !
  !================================================================================================================================
  !

  !>Extracts the error message from a CMISS error string and returns it as a character array
  SUBROUTINE ExtractErrorMessageC(errorMessage,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(OUT) :: errorMessage !<The extracted error message
    INTEGER(INTG), INTENT(IN) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(IN) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: position

    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    errorMessage=EXTRACT(error,1,position-1)

    RETURN
    
  END SUBROUTINE ExtractErrorMessageC

  !
  !================================================================================================================================
  !

  !>Extracts the error stack from a CMISS error string and returns it as a varying string
  SUBROUTINE ExtractErrorStackVS(errorStack,err,error,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(OUT) :: errorStack !<The extracted error stack
    INTEGER(INTG), INTENT(IN) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(IN) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: position

    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    errorStack=EXTRACT(error,position+1,LEN_TRIM(error))

    RETURN
    
  END SUBROUTINE ExtractErrorStackVS

  !
  !================================================================================================================================
  !

  !>Extracts the error stack from a CMISS error string and returns it as a character array
  SUBROUTINE ExtractErrorStackC(errorStack,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(OUT) :: errorStack !<The extracted error stack
    INTEGER(INTG), INTENT(IN) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(IN) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: position

    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    errorStack=EXTRACT(error,position+1,LEN_TRIM(error))

    RETURN
    
  END SUBROUTINE ExtractErrorStackC

  !
  !================================================================================================================================
  !

  !>Sets the error string specified by a character string and flags an error 
  SUBROUTINE FlagErrorC(string,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The error condition string
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: stringLength

    IF(err==0) err=1
    stringLength=LEN_TRIM(string)
    ERROR=string(1:stringLength)

    RETURN 1
    
  END SUBROUTINE FlagErrorC

  !
  !================================================================================================================================
  !

  !>Sets the error string specified by a varying string and flags an error.
  SUBROUTINE FlagErrorVS(string,err,error,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The error condition string
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    IF(err==0) err=1
    error=string

    RETURN 1
    
  END SUBROUTINE FlagErrorVS

  !
  !================================================================================================================================
  !

  !>Writes a warning message specified by a character string to the user.
  SUBROUTINE FlagWarningC(string,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: string !<The warning string
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    IF(numberOfComputationalNodes>1) THEN
      WRITE(outputString,'(">>WARNING (",I0,"): ",A)') myComputationalNodeNumber,string
    ELSE
      WRITE(outputString,'(">>WARNING: ",A)') string
    ENDIF
    CALL WriteStr(WARNING_OUTPUT_TYPE,err,error,*999)

    RETURN 
999 ERRORS("FlagWarningC",err,error)
    RETURN 1
    
  END SUBROUTINE FlagWarningC

  !
  !================================================================================================================================
  !

  !>Writes a warning message specified by a varying string to the user.
  SUBROUTINE FlagWarningVS(string,err,error,*)

    !Argument variables
    TYPE(VARYING_STRING), INTENT(IN) :: string !<The warning string
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    IF(numberOfComputationalNodes>1) THEN
      WRITE(outputString,'(">>WARNING (",I0,"): ",A)') myComputationalNodeNumber,CHAR(string)
    ELSE
      WRITE(outputString,'(">>WARNING: ",A)') CHAR(string)
    ENDIF
    CALL WriteStr(WARNING_OUTPUT_TYPE,err,error,*999)

    RETURN 
999 ERRORS("FlagWarningVS",err,error)
    RETURN 1
    
  END SUBROUTINE FlagWarningVS

  !
  !================================================================================================================================
  !

  !>Finalises the base_routines module and deallocates all memory. \todo Finish this routine and deallocate memory.
  SUBROUTINE BaseRoutines_Finalise(err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    err=0
    error=""
    !Deallocate the random seeds
    IF(ALLOCATED(cmissRandomSeeds)) DEALLOCATE(cmissRandomSeeds)
    
    RETURN 
999 RETURN 1
  END SUBROUTINE BaseRoutines_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises the variables required for the base_routines module.
  SUBROUTINE BaseRoutines_Initialise(err,error,*)
    
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,randomSeedsSize,time(8)

    err=0
    error=""
    myComputationalNodeNumber=0
    numberOfComputationalNodes=1
    diagnostics=.FALSE.
    diagnostics1=.FALSE.
    diagnostics2=.FALSE.
    diagnostics3=.FALSE.
    diagnostics4=.FALSE.
    diagnostics5=.FALSE.
    diagnosticsLevel1=.FALSE.
    diagnosticsLevel2=.FALSE.
    diagnosticsLevel3=.FALSE.
    diagnosticsLevel4=.FALSE.
    diagnosticsLevel5=.FALSE.
    diagAllSubroutines=.TRUE.
    diagFromSubroutine=.FALSE.
    diagFileOpen=.FALSE.
    diagOrTiming=.FALSE.
    echoOutput=.FALSE.
    timing=.FALSE.
    timingSummary=.FALSE.
    timingAllSubroutines=.TRUE.
    timingFromSubroutine=.FALSE.
    timingFileOpen=.FALSE.
    !Initialise loose tolerance here rather than in constants.f90
    LOOSE_TOLERANCE=SQRT(EPSILON(1.0_DP))
    LOOSE_TOLERANCE_SP=SQRT(EPSILON(1.0_SP))
    !Setup the random seeds based on the time
    CALL RANDOM_SEED(SIZE=randomSeedsSize)
    ALLOCATE(cmissRandomSeeds(randomSeedsSize),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate random seeds.",err,error,*999)
    cmissRandomSeeds(1:randomSeedsSize)=[(i,i=1,randomSeedsSize)]
    CALL DATE_AND_TIME(VALUES=time)
    cmissRandomSeeds(1)=3600000*time(5)+60000*time(6)+1000*time(7)+time(8)
    CALL RANDOM_SEED(PUT=cmissRandomSeeds)

    !Initialise outputString
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      DO i=1,MAX_OUTPUT_LINES
        outputString(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      DO i=1,MAX_OUTPUT_LINES
        DO j=1,MAXSTRLEN
          outputString(i)(j:j)=' '
        ENDDO !j
      ENDDO !i
    CASE(WINDOWS_OS)
      DO i=1,MAX_OUTPUT_LINES
        outputString(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Operating system not implemented.",err,error,*999)
    END SELECT

    !Initialise diagnostics and tracing
    NULLIFY(routineStack%stackPointer)
    NULLIFY(diagRoutineList%head)
    NULLIFY(timingRoutineList%head)

    RETURN 
999 RETURN 1
    
  END SUBROUTINE BaseRoutines_Initialise

  !
  !================================================================================================================================
  !

  !>Sets diagnositics off. \see BaseRoutines::DiagnosticsSetOn,OpenCMISS::Iron::DiagnosticsSetOn
  SUBROUTINE DiagnosticsSetOff(err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(RoutineListItemType), POINTER :: nextRoutine,routine

    ENTERS("DiagnosticsSetOff",err,error,*999)

    IF(diagnostics) THEN
      IF(diagFileOpen) THEN
        diagFileOpen=.FALSE.
        CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      ENDIF
      IF(diagAllSubroutines) THEN
        diagAllSubroutines=.FALSE.
      ELSE
        routine=>diagRoutineList%head
        DO WHILE(ASSOCIATED(routine))
          nextRoutine=>routine%nextRoutine
          DEALLOCATE(routine)
          routine=>nextRoutine
        ENDDO
        NULLIFY(diagRoutineList%head)
        diagFromSubroutine=.FALSE.
      ENDIF
      diagnosticsLevel1=.FALSE.
      diagnosticsLevel2=.FALSE.
      diagnosticsLevel3=.FALSE.
      diagnosticsLevel4=.FALSE.
      diagnosticsLevel5=.FALSE.
      diagnostics1=.FALSE.
      diagnostics2=.FALSE.
      diagnostics3=.FALSE.
      diagnostics4=.FALSE.
      diagnostics5=.FALSE.
      diagnostics=.FALSE.
      diagOrTiming=timing
    ELSE
      CALL FlagError("Diagnositics is not on.",err,error,*999)
    ENDIF

    EXITS("DiagnosticsSetOff")
    RETURN
999 ERRORSEXITS("DiagnosticsSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE DiagnosticsSetOff

  !
  !================================================================================================================================
  !

  !>Sets diagnositics on. \see BaseRoutines::DiagnosticsSetOff,OpenCMISS::Iron::DiagnosticsSetOff
  SUBROUTINE DiagnosticsSetOn(diagType,levelList,diagFilename,routineList,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: diagType !<The type of diagnostics to set on \see BaseRoutines_DiagnosticTypes
    INTEGER(INTG), INTENT(IN) :: levelList(:) !<The list of diagnostic levels to set on
    CHARACTER(LEN=*), INTENT(IN) :: diagFilename !<If present the name of the file to output diagnostic information to. If omitted the diagnostic output is sent to the screen
    CHARACTER(LEN=*), INTENT(IN) :: routineList(:) !<The list of routines to set diagnostics on in.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,level
    CHARACTER(LEN=MAXSTRLEN) :: filename
    TYPE(RoutineListItemType), POINTER :: nextRoutine,previousRoutine,routine
 
    NULLIFY(routine)
    
    ENTERS("DiagnosticsSetOn",err,error,*999)

    IF(LEN_TRIM(diagFilename)>=1) THEN
      IF(diagFileOpen) CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      IF(numberOfComputationalNodes>1) THEN
        WRITE(filename,'(A,".diag.",I0)') diagFilename(1:LEN_TRIM(diagFilename)),myComputationalNodeNumber
      ELSE
        filename=diagFilename(1:LEN_TRIM(diagFilename))//".diag"
      ENDIF
      OPEN(UNIT=DIAGNOSTICS_FILE_UNIT,FILE=filename(1:LEN_TRIM(filename)),STATUS="UNKNOWN",IOSTAT=err)
      IF(err/=0) CALL FlagError("Could not open diagnostics file.",err,error,*999)
      diagFileOpen=.TRUE.
    ENDIF
    SELECT CASE(diagType)
    CASE(ALL_DIAG_TYPE)
      diagAllSubroutines=.TRUE.
    CASE(IN_DIAG_TYPE,FROM_DIAG_TYPE)
      diagAllSubroutines=.FALSE.
      diagFromSubroutine=diagType==FROM_DIAG_TYPE
      IF(ASSOCIATED(diagRoutineList%head)) THEN
        routine=>diagRoutineList%head
        DO WHILE(ASSOCIATED(routine))
          nextRoutine=>routine%nextRoutine
          DEALLOCATE(routine)
          routine=>nextRoutine
        ENDDO
        NULLIFY(diagRoutineList%head)
      ENDIF
      ALLOCATE(routine,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate routine list item.",err,error,*999)
      routine%name=routineList(1)
      previousRoutine=>routine
      NULLIFY(routine%nextRoutine)
      diagRoutineList%head=>routine
      DO i=2,SIZE(routineList,1)
        ALLOCATE(routine,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate routine list item.",err,error,*999)
        routine%name=routineList(i)
        NULLIFY(routine%nextRoutine)
        previousRoutine%nextRoutine=>routine
        previousRoutine=>routine
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Invalid diagnostic type.",err,error,*999)
    END SELECT
    DO i=1,SIZE(levelList,1)
      level=levelList(i)
      SELECT CASE(level)
      CASE(1)
        diagnosticsLevel1=.TRUE.
      CASE(2)
        diagnosticsLevel2=.TRUE.
      CASE(3)
        diagnosticsLevel3=.TRUE.
      CASE(4)
        diagnosticsLevel4=.TRUE.
      CASE(5)
        diagnosticsLevel5=.TRUE.
      CASE DEFAULT
        CALL FlagError("Invalid diagnostic level.",err,error,*999)
      END SELECT
    ENDDO !i
    diagnostics=.TRUE.
    diagOrTiming=.TRUE.

    EXITS("DiagnosticsSetOn")
    RETURN
999 IF(diagFileOpen) THEN
      CLOSE(UNIT=DIAGNOSTICS_FILE_UNIT)
      diagFileOpen=.FALSE.
    ENDIF
    routine=>diagRoutineList%head
    DO WHILE(ASSOCIATED(routine))
      nextRoutine=>routine%nextRoutine
      DEALLOCATE(routine)
      routine=>nextRoutine
    ENDDO
    NULLIFY(diagRoutineList%head)
    diagAllSubroutines=.FALSE.
    diagFromSubroutine=.FALSE.
    diagnosticsLevel1=.FALSE.
    diagnosticsLevel2=.FALSE.
    diagnosticsLevel3=.FALSE.
    diagnosticsLevel4=.FALSE.
    diagnosticsLevel5=.FALSE.
    diagnostics=.FALSE.
    diagOrTiming=timing
    ERRORSEXITS("DiagnosticsSetOn",err,error)
    RETURN 1
    
  END SUBROUTINE DiagnosticsSetOn

  !
  !================================================================================================================================
  !

  !>Sets writes file echo output off. \see BaseRoutines::OutputSetOn,OpenCMISS::Iron::OutputSetOff
  SUBROUTINE OutputSetOff(err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables

    ENTERS("OutputSetOff",err,error,*999)

    IF(echoOutput) THEN
      echoOutput=.FALSE.
      CLOSE(UNIT=ECHO_FILE_UNIT)
    ELSE
      CALL FlagError("Write output is not on.",err,error,*999)
    ENDIF

    EXITS("OutputSetOff")
    RETURN
999 ERRORSEXITS("OutputSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE OutputSetOff

  !
  !================================================================================================================================
  !

  !>Sets writes file echo output on. \see BaseRoutines::OutputSetOff,OpenCMISS::Iron::OutputSetOn
  SUBROUTINE OutputSetOn(echoFilename,err,error,*)

    !Argument variables
    CHARACTER(LEN=*), INTENT(IN) :: echoFilename !<The filename of the file to echo output to
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    CHARACTER(LEN=MAXSTRLEN) :: filename

    ENTERS("OutputSetOn",err,error,*999)

    IF(echoOutput) THEN
      CALL FlagError("Write output is already on.",err,error,*999)
    ELSE
      IF(numberOfComputationalNodes>1) THEN
        WRITE(filename,'(A,".out.",I0)') echoFilename(1:LEN_TRIM(echoFilename)),myComputationalNodeNumber        
      ELSE
        filename=echoFilename(1:LEN_TRIM(echoFilename))//".out"
      ENDIF
      OPEN(UNIT=ECHO_FILE_UNIT,FILE=filename(1:LEN_TRIM(filename)),STATUS="UNKNOWN",IOSTAT=err)
      IF(err/=0) CALL FlagError("Could not open write output file.",err,error,*999)
      echoOutput=.TRUE.
    ENDIF

    EXITS("OutputSetOn")
    RETURN
999 ERRORSEXITS("OutputSetOn",err,error)
    RETURN 1
    
  END SUBROUTINE OutputSetOn

  !
  !================================================================================================================================
  !

  !>Returns the random seeds for CMISS \see OpenCMISS::Iron::RandomSeedsGet
  SUBROUTINE RandomSeedsGet(randomSeeds,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: randomSeeds(:) !<On return, the random seeds.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    CHARACTER(LEN=MAXSTRLEN) :: localError
    
    ENTERS("RandomSeedsGet",err,error,*999)

    IF(SIZE(randomSeeds,1)>=SIZE(cmissRandomSeeds,1)) THEN
      randomSeeds(1:SIZE(cmissRandomSeeds,1))=cmissRandomSeeds(1:SIZE(cmissRandomSeeds,1))
    ELSE
      WRITE(localError,'("The size of the supplied random seeds array of ",I2," is too small. The size must be >= ",I2,".")') &
        & SIZE(randomSeeds,1),SIZE(cmissRandomSeeds,1)
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("RandomSeedsGet")
    RETURN
999 ERRORSEXITS("RandomSeedsGet",err,error)
    RETURN 1
    
  END SUBROUTINE RandomSeedsGet

  !
  !================================================================================================================================
  !

  !>Returns the size of the random seeds array for CMISS \see OpenCMISS::Iron::RandomSeedsSizeGet
  SUBROUTINE RandomSeedsSizeGet(randomSeedsSize,err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: randomSeedsSize !<On return, the size of the random seeds array.
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables

    ENTERS("RandomSeedsSizeGet",err,error,*999)

    randomSeedsSize=SIZE(cmissRandomSeeds,1)
    
    EXITS("RandomSeedsSizeGet")
    RETURN
999 ERRORSEXITS("RandomSeedsSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE RandomSeedsSizeGet

  !
  !================================================================================================================================
  !

  !>Sets the random seeds for cmiss \see OpenCMISS::Iron::RandomSeedsSet
  SUBROUTINE RandomSeedsSet(randomSeeds,err,error,*)
  
    !Argument variables 
    INTEGER(INTG), INTENT(IN) :: randomSeeds(:) !<The random seeds to set. 
    INTEGER(INTG), INTENT(INOUT) :: err !<The error string
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error code
    !Local Variables
    
    ENTERS("RandomSeedsSet",err,error,*999)

    IF(SIZE(randomSeeds,1)>SIZE(cmissRandomSeeds,1)) THEN
      cmissRandomSeeds(1:SIZE(cmissRandomSeeds,1))=randomSeeds(1:SIZE(cmissRandomSeeds,1))
    ELSE
      cmissRandomSeeds(1:SIZE(randomSeeds,1))=randomSeeds(1:SIZE(randomSeeds,1))
    ENDIF

    EXITS("RandomSeedsSet")
    RETURN
999 ERRORSEXITS("RandomSeedsSet",err,error)
    RETURN 1
    
  END SUBROUTINE RandomSeedsSet

  !
  !================================================================================================================================
  !

  !>Sets timing off. \see BaseRoutines:TimingSetOn,OpenCMISS::Iron::TimingSetOff
  SUBROUTINE TimingSetOff(err,error,*)

   !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(RoutineListItemType), POINTER :: nextRoutine,routine

    ENTERS("TimingSetOff",err,error,*999)

    IF(timing) THEN
      IF(timingFileOpen) THEN
        timingFileOpen=.FALSE.
        CLOSE(UNIT=TIMING_FILE_UNIT)
      ENDIF
      IF(timingAllSubroutines) THEN
        timingAllSubroutines=.FALSE.
      ELSE
        routine=>timingRoutineList%head
        DO WHILE(ASSOCIATED(routine))
          nextRoutine=>routine%nextRoutine
          DEALLOCATE(routine)
          routine=>nextRoutine
        ENDDO
        NULLIFY(timingRoutineList%head)
        timingFromSubroutine=.FALSE.
      ENDIF
      timingSummary=.FALSE.
      timing=.FALSE.
      diagOrTiming=diagnostics
    ELSE
      CALL FlagError("Timing is not on.",err,error,*999)
    ENDIF

    EXITS("TimingSetOff")
    RETURN
999 ERRORSEXITS("TimingSetOff",err,error)
    RETURN 1
    
  END SUBROUTINE TimingSetOff

  !
  !================================================================================================================================
  !

  !>Sets timing on. \see BaseRoutines:TimingSetOff,OpenCMISS::Iron::TimingSetOn
  SUBROUTINE TimingSetOn(timingType,timingSummaryFlag,timingFilename,routineList,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: timingType !<The type of timing to set on \see BaseRoutines_TimingTypes
    LOGICAL, INTENT(IN) :: timingSummaryFlag !<.TRUE. if the timing information will be output with subsequent TimingSummaryOutput calls, .FALSE. if the timing information will be output every time the routine exits 
    CHARACTER(LEN=*), INTENT(IN) :: timingFilename !<If present the name of the file to output timing information to. If omitted the timing output is sent to the screen
    CHARACTER(LEN=*), INTENT(IN) :: routineList(:) !<The list of routines to set timing on in.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i
    CHARACTER(LEN=MAXSTRLEN) :: filename
    TYPE(RoutineListItemType), POINTER :: nextRoutine,previousRoutine,routine
 
    ENTERS("TimingSetOn",err,error,*999)

    NULLIFY(routine)
    IF(LEN_TRIM(timingFilename)>=1) THEN
      IF(timingFileOpen) CLOSE(UNIT=TIMING_FILE_UNIT)
      IF(numberOfComputationalNodes>1) THEN
        WRITE(filename,'(A,".timing.",I0)') timingFilename(1:LEN_TRIM(timingFilename)),myComputationalNodeNumber
      ELSE
        filename=timingFilename(1:LEN_TRIM(timingFilename))//".timing"
      ENDIF
      OPEN(UNIT=TIMING_FILE_UNIT,FILE=filename(1:LEN_TRIM(filename)),STATUS="UNKNOWN",IOSTAT=err)
      IF(err/=0) CALL FlagError("Could not open timing file.",err,error,*999)
      timingFileOpen=.TRUE.
    ENDIF
    SELECT CASE(timingType)
    CASE(ALL_TIMING_TYPE)
      timingAllSubroutines=.TRUE.
    CASE(IN_TIMING_TYPE,FROM_TIMING_TYPE)
      timingAllSubroutines=.FALSE.
      timingFromSubroutine=timingType==FROM_TIMING_TYPE
      IF(ASSOCIATED(timingRoutineList%head)) THEN
        routine=>timingRoutineList%head
        DO WHILE(ASSOCIATED(routine))
          nextRoutine=>routine%nextRoutine
          DEALLOCATE(routine)
          routine=>nextRoutine
        ENDDO
        NULLIFY(timingRoutineList%head)
      ENDIF
      ALLOCATE(routine,STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate routine list item.",err,error,*999)
      routine%name=routineList(1)
      previousRoutine=>routine
      NULLIFY(routine%nextRoutine)
      timingRoutineList%head=>routine
      routine%numberOfInvocations=0
      routine%totalInclusiveCPUTime=0.0_SP
      routine%totalInclusiveSystemTime=0.0_SP
      routine%totalExclusiveCPUTime=0.0_SP
      routine%totalExclusiveSystemTime=0.0_SP
      DO i=2,SIZE(routineList,1)
        ALLOCATE(routine,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate routine list item.",err,error,*999)
        routine%name=routineList(i)
        NULLIFY(routine%nextRoutine)
        previousRoutine%nextRoutine=>routine
        previousRoutine=>routine
        routine%numberOfInvocations=0
        routine%totalInclusiveCPUTime=0.0_SP
        routine%totalInclusiveCPUTime=0.0_SP
        routine%totalExclusiveCPUTime=0.0_SP
        routine%totalExclusiveCPUTime=0.0_SP
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Invalid timing type.",err,error,*999)
    END SELECT
    timingSummary=timingSummaryFlag
    timing=.TRUE.
    diagOrTiming=.TRUE.

    EXITS("TimingSetOn")
    RETURN
999 IF(timingFileOpen) THEN
      CLOSE(UNIT=TIMING_FILE_UNIT)
      timingFileOpen=.FALSE.
    ENDIF
    routine=>timingRoutineList%head
    DO WHILE(ASSOCIATED(routine))
      nextRoutine=>routine%nextRoutine
      DEALLOCATE(routine)
      routine=>nextRoutine
    ENDDO
    NULLIFY(timingRoutineList%head)
    timingAllSubroutines=.FALSE.
    timingFromSubroutine=.FALSE.
    timing=.FALSE.
    diagOrTiming=diagnostics
    ERRORSEXITS("TimingSetOn",err,error)
    RETURN 1
    
  END SUBROUTINE TimingSetOn

  !
  !================================================================================================================================
  !

  !>Outputs the timing summary. \see OpenCMISS::Iron::TimingSummaryOutput
  SUBROUTINE TimingSummaryOutput(err,error,*)    

    !Argument variables
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(RoutineListItemType), POINTER :: routinePtr

    NULLIFY(routinePtr)
    
    ENTERS("TimingSummaryOutput",err,error,*999)

    IF(timing) THEN
      WRITE(outputString,'("*** Timing Summary: ")') 
      CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
      routinePtr=>timingRoutineList%head
      DO WHILE(ASSOCIATED(routinePtr))
        WRITE(outputString,'("*** Routine : ",A)') TRIM(routinePtr%name)
        CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
        WRITE(outputString,'("***    Number of invocations: ",I10)') routinePtr%numberOfInvocations
        CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
        WRITE(outputString,'("***    Routine times: Total Exclusive  Total Inclusive  Average Exclusive  Average Inclusive")')
        CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
        IF(routinePtr%numberOfInvocations==0) THEN
          WRITE(outputString,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & routinePtr%totalExclusiveCPUTime,routinePtr%totalInclusiveCPUTime, &
            & REAL(routinePtr%numberOfInvocations,SP),REAL(routinePtr%numberOfInvocations,SP)
          CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
          WRITE(outputString,'("***    System    (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & routinePtr%totalExclusiveSystemTime,routinePtr%totalInclusiveSystemTime, &
            & REAL(routinePtr%numberOfInvocations,SP),REAL(routinePtr%numberOfInvocations,SP)
          CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
        ELSE
          WRITE(outputString,'("***    CPU       (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & routinePtr%totalExclusiveCPUTime,routinePtr%totalInclusiveCPUTime, &
            & routinePtr%totalExclusiveCPUTime/REAL(routinePtr%numberOfInvocations,SP), &
            & routinePtr%totalInclusiveCPUTime/REAL(routinePtr%numberOfInvocations,SP)
          CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
          WRITE(outputString,'("***    System    (s):  ",E14.6,"   ",E14.6,"     ",E14.6,"     ",E14.6)')  &
            & routinePtr%totalExclusiveSystemTime,routinePtr%totalInclusiveSystemTime, &
            & routinePtr%totalExclusiveSystemTime/REAL(routinePtr%numberOfInvocations,SP), &
            & routinePtr%totalInclusiveSystemTime/REAL(routinePtr%numberOfInvocations,SP)
          CALL WriteStr(TIMING_OUTPUT_TYPE,err,error,*999)
        ENDIF
        routinePtr=>routinePtr%nextRoutine
      ENDDO
    ELSE
      CALL FlagError("Timing is not on.",err,error,*999)
    ENDIF

    EXITS("TimingSummaryOutput")
    RETURN
999 ERRORSEXITS("TimingSummaryOutput",err,error)
    RETURN 1
    
  END SUBROUTINE TimingSummaryOutput

  !
  !================================================================================================================================
  !

  !>Writes the error string.
  SUBROUTINE WriteError(err,error,*)
  
    !Argument variables
    INTEGER(INTG), INTENT(INOUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(INOUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: endPosition,errorStringLength,indent,lastSpacePosition,localErr,position,startStringLength
    CHARACTER(LEN=MAXSTRLEN) :: indentString=">>"
    CHARACTER(LEN=MAXSTRLEN) :: startString
    TYPE(VARYING_STRING) :: localError,localError2

    indent=2
    IF(numberOfComputationalNodes>1) THEN
      WRITE(startString,'(A,A,I0,A,X,I0,A)') indentString(1:indent),"ERROR (",myComputationalNodeNumber,"):", &
        & ERR,":"
      startStringLength=LEN_TRIM(startString)
    ELSE
      WRITE(startString,'(A,A,X,I0,A)') indentString(1:indent),"ERROR: ",ERR,":"
      startStringLength=LEN_TRIM(startString)
    ENDIF
    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    errorStringLength=position-1
    localError=EXTRACT(error,1,errorStringLength)
    DO WHILE(errorStringLength+startStringLength+1>MAX_OUTPUT_WIDTH)
      endPosition=MAX_OUTPUT_WIDTH-startStringLength-1
      lastSpacePosition=INDEX(EXTRACT(localError,1,endPosition)," ",BACK=.TRUE.)
      IF(lastSpacePosition/=0) endPosition=lastSpacePosition-1
      WRITE(outputString,'(A,X,A)') startString(1:startStringLength),CHAR(EXTRACT(localError,1,endPosition))
      CALL WriteStr(ERROR_OUTPUT_TYPE,localErr,localError2,*999)
      localError=ADJUSTL(EXTRACT(localError,endPosition+1,LEN_TRIM(localError)))
      errorStringLength=LEN_TRIM(localError)
      startString=" "
    ENDDO !not finished
    WRITE(outputString,'(A,X,A)') startString(1:startStringLength),CHAR(localError)
    CALL WriteStr(ERROR_OUTPUT_TYPE,localErr,localError2,*999)
    !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
    localError=REMOVE(error,1,position)
    error=localError
    position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
    indent=indent+2
    DO WHILE(position/=0)
      WRITE(outputString,'(A)') indentString(1:indent)//CHAR(EXTRACT(error,1,position-1))
      CALL WriteStr(ERROR_OUTPUT_TYPE,localErr,localError2,*999)
      !CPB 20/02/07 aix compiler does not like varying strings so split the remove statement up into two statements
      localError=REMOVE(error,1,position)
      error=localError
      position=INDEX(error,ERROR_SEPARATOR_CONSTANT)
      indent=indent+2
    ENDDO
    WRITE(outputString,'(A)') indentString(1:indent)//CHAR(error)
    CALL WriteStr(ERROR_OUTPUT_TYPE,localErr,localError2,*999)

    RETURN
    !Don't return an error code here otherwise we will get into a circular loop
999 RETURN 
  END SUBROUTINE WriteError

  !
  !================================================================================================================================
  !

  !>Writes the output string to a specified output stream.
  SUBROUTINE WriteStr(outputStreamID,err,error,*)


!!!!NOTE: No enters or exits is used here to avoid an infinite loop
!!!!NOTE: This routine is, in general, OS dependent but needs to be defined here so that an module loop is avoided when MACHINE
!!!!      module routines need to use this module.

    !Argument Variables
    INTEGER(INTG), INTENT(IN) :: outputStreamID !<The outputStreamID of the output stream. An outputStreamID of > 9 specifies file output \see BaseRoutines_OutputType,BaseRoutines_FileUnits
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: endLine(MAX_OUTPUT_LINES),i,j,length,numberBlanks,numberRecords

    !Calculate number of records in outputString
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      i=1 
      DO WHILE(outputString(i)(1:1)/=CHAR(0).AND.i<MAX_OUTPUT_LINES)
        i=i+1
      ENDDO
      numberRecords=i-1
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      i=1
      numberBlanks=0
      DO WHILE(i<MAX_OUTPUT_LINES.AND.numberBlanks<2)
        i=i+1
        length=LEN_TRIM(outputString(i))
        IF(length==0) THEN
          numberBlanks=numberBlanks+1
        ELSE
          numberBlanks=0
        ENDIF
      ENDDO
      numberRecords=i-numberBlanks
    CASE(WINDOWS_OS)
      i=1 
      DO WHILE(outputString(i)(1:1)/=CHAR(0).AND.i<MAX_OUTPUT_LINES)
        i=i+1
      ENDDO
      numberRecords=i-1
    CASE DEFAULT
      CALL FlagError("Operating system not implemented.",err,error,*999)
    END SELECT

    DO i=1,numberRecords
      endLine(i)=LEN_TRIM(outputString(i))
    ENDDO !i

    IF(diagFileOpen.AND.outputStreamID==DIAGNOSTIC_OUTPUT_TYPE) THEN
      DO i=1,numberRecords
        IF(endLine(i)<=MAX_OUTPUT_WIDTH) THEN
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') outputString(i)(1:endLine(i))
        ELSE IF(endLine(i)>MAX_OUTPUT_WIDTH.AND.endLine(i)<=MAXSTRLEN) THEN
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:endLine(i))
        ELSE
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
          WRITE(DIAGNOSTICS_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:MAXSTRLEN)
        ENDIF
      ENDDO !i
    ELSE IF(timingFileOpen.AND.outputStreamID==TIMING_OUTPUT_TYPE) THEN
      DO i=1,numberRecords
        IF(endLine(i)<=MAX_OUTPUT_WIDTH) THEN
          WRITE(TIMING_FILE_UNIT,'(A)') outputString(i)(1:endLine(i))
        ELSE IF(endLine(i)>MAX_OUTPUT_WIDTH.AND.endLine(i)<=MAXSTRLEN) THEN
          WRITE(TIMING_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
          WRITE(TIMING_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:endLine(i))
        ELSE
          WRITE(TIMING_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
          WRITE(TIMING_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:MAXSTRLEN)
        ENDIF
      ENDDO !i
    ELSE
      IF(outputStreamID<=9) THEN !not file output
        DO i=1,numberRecords
          IF(endLine(i)<=MAX_OUTPUT_WIDTH) THEN
            WRITE(*,'(A)') outputString(i)(1:endLine(i))
          ELSE IF(endLine(i)>MAX_OUTPUT_WIDTH.AND.endLine(i)<=MAXSTRLEN) THEN
            WRITE(*,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
            WRITE(*,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:endLine(i))
          ELSE
            WRITE(*,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
            WRITE(*,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:MAXSTRLEN)
          ENDIF
        ENDDO !i
      ELSE !file output
        DO i=1,numberRecords
          WRITE(outputStreamID,'(A)') outputString(i)(1:endLine(i))
        ENDDO !i
      ENDIF

      !Echo strings to output file if required

      IF(echoOutput) THEN
        DO i=1,numberRecords
          IF(endLine(i)<=MAX_OUTPUT_WIDTH) THEN
            WRITE(ECHO_FILE_UNIT,'(A)') outputString(i)(1:endLine(i))
          ELSE IF(endLine(i)>MAX_OUTPUT_WIDTH.AND.endLine(i)<=MAXSTRLEN) THEN
            WRITE(ECHO_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
            WRITE(ECHO_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:endLine(i))
          ELSE
            WRITE(ECHO_FILE_UNIT,'(A)') outputString(i)(1:MAX_OUTPUT_WIDTH)
            WRITE(ECHO_FILE_UNIT,'(A)') outputString(i)(MAX_OUTPUT_WIDTH+1:MAXSTRLEN)
          ENDIF
        ENDDO !i
      ENDIF
    ENDIF

    !Reset outputString
    SELECT CASE(MACHINE_OS)
    CASE(VMS_OS)
      DO i=1,numberRecords
        outputString(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE(IRIX_OS,LINUX_OS,AIX_OS)
      DO i=1,numberRecords
        DO j=1,MAXSTRLEN
          outputString(i)(j:j)=' '
        ENDDO !j
      ENDDO !i
    CASE(WINDOWS_OS)
      DO i=1,numberRecords
        outputString(i)(1:1)=CHAR(0)
      ENDDO !i
    CASE DEFAULT
      CALL FlagError("Operating system not implemented.",err,error,*999)
    END SELECT

    RETURN
999 ERRORS("WriteStr",err,error)
    RETURN 1
    
  END SUBROUTINE WriteStr

  !
  !================================================================================================================================
  !

END MODULE BaseRoutines
