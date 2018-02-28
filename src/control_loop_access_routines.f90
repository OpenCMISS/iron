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
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !> \addtogroup ControlLoopRoutines_ControlLoopIdentifiers ControlLoopRoutines::ControlLoopIdentifiers
  !> \brief The control loop identification parameters
  !>@{
  INTEGER(INTG), PARAMETER :: CONTROL_LOOP_NODE=0 !<The identifier for a each "leaf" node in a control loop. \see ControlLoopRoutines_ControlLoopIdentifiers,ControlLoopRoutines
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
  
  INTERFACE CONTROL_LOOP_SOLVERS_GET
    MODULE PROCEDURE ControlLoop_SolversGet
  END INTERFACE CONTROL_LOOP_SOLVERS_GET

  PUBLIC CONTROL_LOOP_NODE

  PUBLIC ControlLoop_Get

  PUBLIC CONTROL_LOOP_GET

  PUBLIC ControlLoop_ProblemGet

  PUBLIC ControlLoop_SolversGet

  PUBLIC CONTROL_LOOP_SOLVERS_GET

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns the specified control loop as indexed by the control loop identifier from the control loop root.
  SUBROUTINE ControlLoop_Get0(controlLoopRoot,controlLoopIdentifier,controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoopRoot!<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On exit, the specified control loop
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
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoopRoot !<A pointer to the control loop to root
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<controlLoopIdentifiers(identifierIdx). The control loop identifiers
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On exit, the specified control loop. Must not be associated on entry.
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
          & controlLoopIdentifiers(controlLoopIdx)<=controlLoop%NUMBER_OF_SUB_LOOPS) THEN
          controlLoop=>controlLoop%SUB_LOOPS(controlLoopIdentifiers(controlLoopIdx))%ptr
          IF(.NOT.ASSOCIATED(controlLoop)) THEN
            localError="Control sub loop number "//TRIM(NumberToVString(controlLoopIdentifiers(controlLoopIdx),"*",err,error))// &
              & " at identifier index "//TRIM(NumberToVString(controlLoopIdx,"*",err,error))//" is not associated."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ELSE
          localError="Invalid control loop identifier. The identifier at index "// &
            & TRIM(NumberToVString(controlLoopIdx,"*",err,error))//" is "// &
            & TRIM(NumberToVString(controlLoopIdentifiers(controlLoopIdx),"*",err,error))// &
            & ". The identifier must be between 1 and "//TRIM(NumberToVString(controlLoop%NUMBER_OF_SUB_LOOPS,"*",err,error))//"."
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

  !>Returns a pointer to the problem for a control loop.
  SUBROUTINE ControlLoop_ProblemGet(controlLoop,problem,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the problem for.
    TYPE(PROBLEM_TYPE), POINTER :: problem !<On exit, a pointer to the control loop problem. Must not be associated on entry.
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

  !>Returns a pointer to the solvers for a control loop.
  SUBROUTINE ControlLoop_SolversGet(controlLoop,solvers,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to get the solvers for.
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

END MODULE ControlLoopAccessRoutines
