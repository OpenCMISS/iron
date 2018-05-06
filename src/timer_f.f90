!> \file
!> \author Chris Bradley
!> \brief This module contains routines for timing the program.
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

!> This module contains routines for timing the program.
MODULE Timer

  USE BaseRoutines
  USE Constants
  USE ISO_C_BINDING
  USE ISO_VARYING_STRING

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module variables
  
  !CPU time parameters
  !> \addtogroup TIMER_TimerType TIMER::TimerType
  !> \brief Timer type parameter
  !> \see TIMER
  !>@{  
  INTEGER(INTG), PARAMETER :: USER_CPU=1 !<User CPU time type \see TIMER_TimerType,TIMER
  INTEGER(INTG), PARAMETER :: SYSTEM_CPU=2 !<System CPU time type \see TIMER_TimerType,TIMER
  INTEGER(INTG), PARAMETER :: TOTAL_CPU=3 !<Total CPU (i.e. User + System) time type \see TIMER_TimerType,TIMER
  !>@}
  
  !Module variables

  !Interfaces

  INTERFACE CPU_TIMER
    MODULE PROCEDURE CPUTimer
  END INTERFACE CPU_TIMER
  
  INTERFACE

    SUBROUTINE CPUTimer_(returnTime, timeType, err, cError) BIND(C,NAME="CPUTimer")
      USE ISO_C_BINDING
      REAL(C_DOUBLE), INTENT(OUT) :: returnTime
      INTEGER(C_INT), INTENT(IN) :: timeType
      INTEGER(C_INT), INTENT(OUT) :: err
      CHARACTER(C_CHAR), INTENT(OUT) :: cerror(*)
    END SUBROUTINE CPUTIMER_

  END INTERFACE

  PUBLIC USER_CPU,SYSTEM_CPU,TOTAL_CPU

  PUBLIC CPU_TIMER

  PUBLIC CPUTimer

CONTAINS

  !
  !============================================================================
  !

  !>CPUTimer returns the CPU time in time(1). timeType indicates the type of time required.
  SUBROUTINE CPUTimer(timeType,time,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: timeType !<The required time type \see Timer_TimerType,Timer
    REAL(SP), INTENT(OUT) :: time(*) !<On return, the requested time.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    REAL(DP) :: returnTime
    CHARACTER(KIND=C_CHAR,LEN=MAXSTRLEN) :: cError

    ENTERS("CPUTimer",err,error,*999)
    
    CALL CPUTimer_(returnTime,timeType,err,cError)
    time(1)=REAL(returnTime,SP)
    IF(err/=0) CALL FlagError(cError,err,error,*999)
    
    EXITS("CPUTimer")
    RETURN
999 ERRORSEXITS("CPUTimer",err,error)
    RETURN 1
    
  END SUBROUTINE CPUTimer

  !
  !============================================================================
  !

END MODULE Timer
