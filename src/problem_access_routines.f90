!> \file
!> \author Chris Bradley
!> \brief This module contains all problem access method routines.
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

!> This module contains all problem access method routines.
MODULE ProblemAccessRoutines
  
  USE BaseRoutines
  USE ControlLoopAccessRoutines
  USE Kinds
  USE SolverAccessRoutines
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  TYPE(PROBLEMS_TYPE), TARGET :: problems
  
  !Interfaces

  INTERFACE Problem_CellMLEquationsGet
    MODULE PROCEDURE Problem_CellMLEquationsGet0
    MODULE PROCEDURE Problem_CellMLEquationsGet1
  END INTERFACE Problem_CellMLEquationsGet

  INTERFACE PROBLEM_CELLML_EQUATIONS_GET   
    MODULE PROCEDURE Problem_CellMLEquationsGet0
    MODULE PROCEDURE Problem_CellMLEquationsGet1
  END INTERFACE PROBLEM_CELLML_EQUATIONS_GET

  INTERFACE Problem_ControlLoopGet
    MODULE PROCEDURE Problem_ControlLoopGet0
    MODULE PROCEDURE Problem_ControlLoopGet1
  END INTERFACE Problem_ControlLoopGet

  INTERFACE PROBLEM_CONTROL_LOOP_GET
    MODULE PROCEDURE Problem_ControlLoopGet0
    MODULE PROCEDURE Problem_ControlLoopGet1
  END INTERFACE PROBLEM_CONTROL_LOOP_GET

  INTERFACE Problem_SolverGet
    MODULE PROCEDURE Problem_SolverGet0
    MODULE PROCEDURE Problem_SolverGet1
  END INTERFACE Problem_SolverGet
  
  INTERFACE PROBLEM_SOLVER_GET
    MODULE PROCEDURE Problem_SolverGet0
    MODULE PROCEDURE Problem_SolverGet1
  END INTERFACE PROBLEM_SOLVER_GET
  
  INTERFACE Problem_SolverEquationsGet
    MODULE PROCEDURE Problem_SolverEquationsGet0
    MODULE PROCEDURE Problem_SolverEquationsGet1
  END INTERFACE Problem_SolverEquationsGet

  INTERFACE PROBLEM_SOLVER_EQUATIONS_GET
    MODULE PROCEDURE Problem_SolverEquationsGet0
    MODULE PROCEDURE Problem_SolverEquationsGet1
  END INTERFACE PROBLEM_SOLVER_EQUATIONS_GET

  INTERFACE PROBLEM_USER_NUMBER_FIND
    MODULE PROCEDURE Problem_UserNumberFind
  END INTERFACE PROBLEM_USER_NUMBER_FIND

  PUBLIC problems

  PUBLIC Problem_CellMLEquationsGet

  PUBLIC PROBLEM_CELLML_EQUATIONS_GET

  PUBLIC Problem_ControlLoopGet

  PUBLIC PROBLEM_CONTROL_LOOP_GET
  
  PUBLIC Problem_ControlLoopRootGet

  PUBLIC Problem_Get

  PUBLIC Problem_SolverGet

  PUBLIC PROBLEM_SOLVER_GET

  PUBLIC Problem_SolverEquationsGet

  PUBLIC PROBLEM_SOLVER_EQUATIONS_GET

  PUBLIC Problem_UserNumberFind

  PUBLIC PROBLEM_USER_NUMBER_FIND

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_CellMLEquationsGet
  SUBROUTINE Problem_CellMLEquationsGet0(problem,controlLoopIdentifier,solverIndex,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier to get the solver CellML equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<On exit, a pointer to the specified solver CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_CellMLEquationsGet0",err,error,*999)

    CALL Problem_CellMLEquationsGet1(Problem,[controlLoopIdentifier],solverIndex,cellMLEquations,err,error,*999)
    
    EXITS("Problem_CellMLEquationsGet0")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsGet0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver CellML equations defined with a solver. \see OPENCMISS::CMISSProblemSolverCellMLEquationsGet
  SUBROUTINE Problem_CellMLEquationsGet1(problem,controlLoopIdentifiers,solverIndex,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier to get the CellML equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,controlLoopRoot
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("Problem_CellMLEquationsGet1",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*999)

    NULLIFY(controlLoopRoot)
    NULLIFY(controlLoop)
    NULLIFY(solvers)
    NULLIFY(solver)

    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    
    EXITS("Problem_CellMLEquationsGet1")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsGet1",err,error)
    RETURN 1
  END SUBROUTINE Problem_CellMLEquationsGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopGet
  SUBROUTINE Problem_ControlLoopGet0(problem,controlLoopIdentifier,controlLoop,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_ControlLoopGet0",err,error,*999)

    CALL Problem_ControlLoopGet1(problem,[controlLoopIdentifier],controlLoop,err,error,*999) 
       
    EXITS("Problem_ControlLoopGet0")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopGet0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control_loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopGet
  SUBROUTINE Problem_ControlLoopGet1(problem,controlLoopIdentifiers,controlLoop,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get the control loop for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On return, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_ControlLoopGet1",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(problem%CONTROL_LOOP)) CALL FlagError("Problem control loop is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL ControlLoop_Get(problem%CONTROL_LOOP,controlLoopIdentifiers,controlLoop,err,error,*999)
    
    EXITS("Problem_ControlLoopGet1")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the root control loop for a problem. 
  SUBROUTINE Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<a pointer to the problem to get the control loop root for.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoopRoot !<On return, a pointer to the control loop root for the problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_ControlLoopRootGet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(controlLoopRoot)) CALL FlagError("Control loop root is already associated.",err,error,*999)

    controlLoopRoot=>problem%CONTROL_LOOP
    IF(.NOT.ASSOCIATED(controlLoopRoot)) THEN
      localError="The problem control loop root is not associated for problem number "// &
        & TRIM(NumberToVString(problem%USER_NUMBER,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Problem_ControlLoopRootGet")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopRootGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopRootGet

  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the problem with the given user number. 
  SUBROUTINE Problem_Get(userNumber,problem,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the problem to find
    TYPE(PROBLEM_TYPE), POINTER :: problem !<On exit, a pointer to the problem with the specified user number if it exists. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_Get",err,error,*999)

    CALL Problem_UserNumberFind(userNumber,problem,err,error,*999)
    IF(.NOT.ASSOCIATED(problem)) THEN
      localError="A problem with an user number of "//TRIM(NumberToVString(userNumber,"*",err,error))//" does not exist."
      CALL FlagError(localError,err,error,*999)
    ENDIF
  
    EXITS("Problem_Get")
    RETURN
999 ERRORSEXITS("Problem_Get",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Get

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OpenCMISS::Iron::cmfe_Problem_SolverGet
  SUBROUTINE Problem_SolverGet0(problem,controlLoopIdentifier,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: solver !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_SolverGet0",err,error,*999)

    CALL Problem_SolverGet1(problem,[controlLoopIdentifier],solverIndex,solver,err,error,*999) 
       
    EXITS("Problem_SolverGet0")
    RETURN
999 ERRORSEXITS("Problem_SolverGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverGet0
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a problem control loop. \see OpenCMISS::Iron::cmfe_Problem_SolverGet
  SUBROUTINE Problem_SolverGet1(problem,controlLoopIdentifiers,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get the solver for.
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<The control loop identifier to get the solver for.
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index to get the solver for.
    TYPE(SOLVER_TYPE), POINTER :: solver !<On return, a pointer to the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,controlLoopRoot
    TYPE(SOLVERS_TYPE), POINTER :: solvers
 
    ENTERS("Problem_SolverGet1",err,error,*998)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*999)

    NULLIFY(controlLoopRoot)
    NULLIFY(controlLoop)
    NULLIFY(solvers)
    NULLIFY(solver)

    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
     
    EXITS("Problem_SolverGet1")
    RETURN
999 NULLIFY(SOLVER)
998 ERRORSEXITS("Problem_SolverGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverGet1
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsGet
  SUBROUTINE Problem_SolverEquationsGet0(problem,controlLoopIdentifier,solverIndex,solverEquations,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifier !<The control loop identifier to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Problem_SolverEquationsGet0",err,error,*999)

    CALL Problem_SolverEquationsGet1(problem,[controlLoopIdentifier],solverIndex,solverEquations,err,error,*999)
    
    EXITS("Problem_SolverEquationsGet0")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsGet0",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsGet0

  !
  !================================================================================================================================
  !

  !>Returns a pointer to a solver equations defined with a solver. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsGet
  SUBROUTINE Problem_SolverEquationsGet1(problem,controlLoopIdentifiers,solverIndex,solverEquations,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to get solver equations for
    INTEGER(INTG), INTENT(IN) :: controlLoopIdentifiers(:) !<controlLoopIdentifiers(identifierIdx). The control loop identifiers to get the solver equations for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in the solvers to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,controlLoopRoot
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_SolverEquationsGet1",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*999)

    NULLIFY(controlLoopRoot)
    NULLIFY(controlLoop)
    NULLIFY(solvers)
    NULLIFY(solver)

    CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
    CALL ControlLoop_Get(controlLoopRoot,controlLoopIdentifiers,controlLoop,err,error,*999)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    CALL Solvers_SolverGet(solvers,solverIndex,solver,err,error,*999)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        
    EXITS("Problem_SolverEquationsGet1")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsGet1",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsGet1
  
  !
  !================================================================================================================================
  !

  !>Finds and returns a pointer to the problem identified by a user number. If no problem with that user number exists problem is left nullified.
  SUBROUTINE Problem_UserNumberFind(userNumber,problem,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number to find.
    TYPE(PROBLEM_TYPE), POINTER :: problem !<On return, a pointer to the problem with the given user number. If no problem with that user number exists then the pointer is returned as NULL. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_UserNumberFind",err,error,*999)

    IF(ASSOCIATED(problem)) CALL FlagError("Problem is already associated.",err,error,*999)
   
    NULLIFY(problem)
    IF(ASSOCIATED(problems%problems)) THEN
      DO problemIdx=1,problems%NUMBER_OF_PROBLEMS
        IF(ASSOCIATED(problems%problems(problemIdx)%ptr)) THEN
          IF(problems%problems(problemIdx)%ptr%USER_NUMBER==userNumber) THEN
            problem=>problems%problems(problemIdx)%ptr
            EXIT
          ENDIF
        ELSE
          localError="The problem pointer in problems is not associated for problem index "// &
            & TRIM(NumberToVString(problemIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !problemIdx
    ENDIF
    
    EXITS("Problem_UserNumberFind")
    RETURN
999 ERRORSEXITS("Problem_UserNumberFind",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_UserNumberFind

  !
  !================================================================================================================================
  !
 
END MODULE ProblemAccessRoutines
