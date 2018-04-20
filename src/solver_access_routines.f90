!> \file
!> \author Chris Bradley
!> \brief This module contains all solver access method routines.
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

!> This module contains all solver access method routines.
MODULE SolverAccessRoutines
  
  USE BaseRoutines
  USE Kinds
  USE Strings
  USE Types
  
#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  INTERFACE SOLVER_CELLML_EQUATIONS_GET
    MODULE PROCEDURE Solver_CellMLEquationsGet
  END INTERFACE SOLVER_CELLML_EQUATIONS_GET

  INTERFACE SOLVER_SOLVER_EQUATIONS_GET
    MODULE PROCEDURE Solver_SolverEquationsGet
  END INTERFACE SOLVER_SOLVER_EQUATIONS_GET

  INTERFACE SOLVERS_SOLVER_GET
    MODULE PROCEDURE Solvers_SolverGet
  END INTERFACE SOLVERS_SOLVER_GET
  
  INTERFACE SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET
    MODULE PROCEDURE SolverEquations_BoundaryConditionsGet
  END INTERFACE SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET

  PUBLIC CellMLEquations_SolverGet
  
  PUBLIC Solver_CellMLEquationsGet

  PUBLIC SOLVER_CELLML_EQUATIONS_GET

  PUBLIC Solver_ControlLoopGet
  
  PUBLIC Solver_SolverEquationsGet

  PUBLIC SOLVER_SOLVER_EQUATIONS_GET

  PUBLIC Solver_SolversGet

  PUBLIC Solvers_ControlLoopGet

  PUBLIC Solvers_SolverGet

  PUBLIC SOLVERS_SOLVER_GET

  PUBLIC SolverEquations_BoundaryConditionsGet

  PUBLIC SOLVER_EQUATIONS_BOUNDARY_CONDITIONS_GET
  
  PUBLIC SolverEquations_SolverGet

  PUBLIC SolverEquations_SolverMappingGet

  PUBLIC SolverEquations_SolverMatricesGet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver for a CellML equations.
  SUBROUTINE CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*)

    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<A pointer to the CellML equations to get the solver for
    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("CellMLEquations_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)

    solver=>cellMLEquations%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("CellML equations solver is not associated.",err,error,*999)
      
    EXITS("CellMLEquations_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("CellMLEquations_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_SolverGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the CellML equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_CellMLEquationsGet
  SUBROUTINE Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the CellML equations for
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<On exit, a pointer to the specified CellML equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_CellMLEquationsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    !IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)

    cellMLEquations=>solver%CELLML_EQUATIONS
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("Solver CellML equations is not associated.",err,error,*999)
      
    EXITS("Solver_CellMLEquationsGet")
    RETURN
999 NULLIFY(cellMLEquations)
998 ERRORSEXITS("Solver_CellMLEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_CellMLEquationsGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the control loop for a solver.
  SUBROUTINE Solver_ControlLoopGet(solver,controlLoop,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the control loop for
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On exit, a pointer to the control loop for the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers
 
    ENTERS("Solver_ControlLoopGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)

    solvers=>solver%solvers
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solver solvers is not associated.",err,error,*999)
    controlLoop=>solvers%CONTROL_LOOP
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Solvers control loop is not associated.",err,error,*999)
    
    EXITS("Solver_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solver_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_ControlLoopGet
  
  !
  !================================================================================================================================
  !

  !>Returns a pointer to the solver equations for a solver. \see OpenCMISS::Iron::cmfe_Solver_SolverEquationsGet
  SUBROUTINE Solver_SolverEquationsGet(solver,solverEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<On exit, a pointer to the specified solver equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolverEquationsGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(.NOT.solver%SOLVER_FINISHED) CALL FlagError("Solver has not been finished.",err,error,*998)
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)

    solverEquations=>solver%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver solver equations is not associated.",err,error,*999)
      
    EXITS("Solver_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("Solver_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolverEquationsGet
  
  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the solvers for a solver.
  SUBROUTINE Solver_SolversGet(solver,solvers,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to get the solvers for.
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<On exit, A pointer to the solvers for the solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solver_SolversGet",err,error,*998)

    IF(ASSOCIATED(solvers)) CALL FlagError("Solvers is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
      
    solvers=>solver%solvers
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("The solver solvers is not associated.",err,error,*999)
       
    EXITS("Solver_SolversGet")
    RETURN
999 NULLIFY(solvers)
998 ERRORSEXITS("Solver_SolversGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolversGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the control loop for a solvers.
  SUBROUTINE Solvers_ControlLoopGet(solvers,controlLoop,err,error,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<A pointer to the solvers to get the control loop for.
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<On exit, A pointer to the control loop for the solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("Solvers_ControlLoopGet",err,error,*998)

    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
      
    controlLoop=>solvers%CONTROL_LOOP
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("The solvers control loop is not associated.",err,error,*999)
       
    EXITS("Solvers_ControlLoopGet")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("Solvers_ControlLoopGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_ControlLoopGet

  !
  !================================================================================================================================
  !
  
  !>Returns a pointer to the specified solver in the list of solvers.
  SUBROUTINE Solvers_SolverGet(solvers,solverIndex,solver,err,error,*)

    !Argument variables
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<A pointer to the solvers to get the solver for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The specified solver to get
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the specified solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solvers_SolverGet",err,error,*998)

    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*998)
    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(solverIndex<=0.OR.solverIndex>solvers%NUMBER_OF_SOLVERS) THEN
      localError="The specified solver index of "//TRIM(NumberToVString(solverIndex,"*",err,error))// &
        & " is invalid. The solver index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solvers%NUMBER_OF_SOLVERS,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(solvers%solvers)) CALL FlagError("Solvers solvers is not associated.",err,error,*998)
      
    solver=>solvers%solvers(solverIndex)%ptr
    IF(.NOT.ASSOCIATED(solver)) THEN
      localError="The solvers solver is not associated for solver index "// &
        & TRIM(NumberToVString(solverIndex,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
       
    EXITS("Solvers_SolverGet")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("Solvers_SolverGet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_SolverGet

  !
  !================================================================================================================================
  !

  !>Gets the boundary conditions for solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_BoundaryConditionsGet
  SUBROUTINE SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions !<On exit, a pointer to the boundary conditions for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_BoundaryConditionsGet",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)
    IF(ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is already associated.",err,error,*999)

    boundaryConditions=>solverEquations%BOUNDARY_CONDITIONS
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Solver equations boundary conditions is not associated.", &
      & err,error,*999)
 
    EXITS("SolverEquations_BoundaryConditionsGet")
    RETURN
999 ERRORSEXITS("SolverEquations_BoundaryConditionsGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_BoundaryConditionsGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver for solver equations. 
  SUBROUTINE SolverEquations_SolverGet(solverEquations,solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver for
    TYPE(SOLVER_TYPE), POINTER :: solver !<On exit, a pointer to the solver for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverGet",err,error,*998)

    IF(ASSOCIATED(solver)) CALL FlagError("Solver is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    solver=>solverEquations%solver
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver equations solver is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverGet")
    RETURN
998 NULLIFY(solver)
999 ERRORSEXITS("SolverEquations_SolverGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver mapping for solver equations. 
  SUBROUTINE SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver mapping for
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<On exit, a pointer to the solver mapping for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMappingGet",err,error,*998)

    IF(ASSOCIATED(solverMapping)) CALL FlagError("Solver mapping is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)

    solverMapping=>solverEquations%SOLVER_MAPPING
    IF(.NOT.ASSOCIATED(solverMapping)) CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverMappingGet")
    RETURN
998 NULLIFY(solverMapping)
999 ERRORSEXITS("SolverEquations_SolverMappingGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMappingGet
     
  !
  !================================================================================================================================
  !

  !>Gets the solver matrices for solver equations. 
  SUBROUTINE SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the solver matrices for
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices !<On exit, a pointer to the solver matrices for the specified solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_SolverMatricesGet",err,error,*998)

    IF(ASSOCIATED(solverMatrices)) CALL FlagError("Solver matrices is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    IF(.NOT.solverEquations%SOLVER_EQUATIONS_FINISHED) CALL FlagError("Solver equations has not been finished.",err,error,*999)

    solverMatrices=>solverEquations%SOLVER_MATRICES
    IF(.NOT.ASSOCIATED(solverMatrices)) CALL FlagError("Solver equations solver matrices is not associated.",err,error,*999)
 
    EXITS("SolverEquations_SolverMatricesGet")
    RETURN
998 NULLIFY(solverMatrices)
999 ERRORSEXITS("SolverEquations_SolverMatricesGet",err,error)
    RETURN 1

  END SUBROUTINE SolverEquations_SolverMatricesGet
     
  !
  !================================================================================================================================
  !

END MODULE SolverAccessRoutines
