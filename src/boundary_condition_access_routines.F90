!> \file
!> \author Chris Bradley
!> \brief This module contains all boundary condition access method routines.
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

!> This module contains all boundary condition access method routines.
MODULE BoundaryConditionAccessRoutines
  
  USE BaseRoutines
  USE ISO_VARYING_STRING
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

  PUBLIC BoundaryConditions_SolverEquationsGet
  
  PUBLIC BoundaryConditionsVariable_DirichletConditionsGet

  PUBLIC BoundaryConditionsVariable_NeumannConditionsGet

  PUBLIC BoundaryConditionsVariable_PressureIncConditionsGet  
 
CONTAINS

  !
  !================================================================================================================================
  !

  !>Gets the solver equations for boundary conditions.
  SUBROUTINE BoundaryConditions_SolverEquationsGet(boundaryConditions,solverEquations,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions !<A pointer to the boundary conditions to get the solver equations for
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<On exit, a pointer to the solver equations in the specified boundary conditions. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditions_SolverEquationsGet",err,error,*998)

    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

    solverEquations=>boundaryConditions%SOLVER_EQUATIONS
    IF(.NOT.ASSOCIATED(solverEquations)) &
      & CALL FlagError("Solver equations is not associated for the boundary conditions.",err,error,*999)
       
    EXITS("BoundaryConditions_SolverEquationsGet")
    RETURN
999 NULLIFY(solverEquations)
998 ERRORSEXITS("BoundaryConditions_SolverEquationsGet",err,error)
    RETURN 1
    
  END SUBROUTINE BoundaryConditions_SolverEquationsGet

  !
  !================================================================================================================================
  !

  !>Gets the dirichlet boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,dirichletBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the dirichlet boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: dirichletBoundaryCOnditions !<On exit, a pointer to the dirichlet boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_DirichletConditionsGet",err,error,*998)

    IF(ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    dirichletBoundaryConditions=>boundaryConditionsVariable%DIRICHLET_BOUNDARY_CONDITIONS
    IF(.NOT.ASSOCIATED(dirichletBoundaryConditions)) &
      & CALL FlagError("Dirichlet boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
       
    EXITS("BoundaryConditionsVariable_DirichletConditionsGet")
    RETURN
999 NULLIFY(dirichletBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_DirichletConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_DirichletConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_DirichletConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the neumann boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,neumannBoundaryConditions, &
    & err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the neumann boundary conditions for
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryCOnditions !<On exit, a pointer to the neumann boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_NeumannConditionsGet",err,error,*998)

    IF(ASSOCIATED(neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    neumannBoundaryConditions=>boundaryConditionsVariable%neumannBoundaryConditions
    IF(.NOT.ASSOCIATED(neumannBoundaryConditions)) &
      & CALL FlagError("Neumann boundary conditions is not associated for the boundary conditions variable.",err,error,*999)
       
    EXITS("BoundaryConditionsVariable_NeumannConditionsGet")
    RETURN
999 NULLIFY(neumannBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_NeumannConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_NeumannConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_NeumannConditionsGet

  !
  !================================================================================================================================
  !

  !>Gets the pressure incremented boundary conditions for a boundary condition variable.
  SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsGet(boundaryConditionsVariable, &
    & pressureIncrementedBoundaryConditions,err,error,*)

    !Argument variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable !<A pointer to the boundary conditions variable to get the pressure incremented boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: pressureIncrementedBoundaryConditions !<On exit, a pointer to the pressure incremented boundary conditions in the specified boundary conditions variable. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("BoundaryConditionsVariable_PressureIncConditionsGet",err,error,*998)

    IF(ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(boundaryConditionsVariable)) &
      & CALL FlagError("Boundary conditions variable is not associated.",err,error,*999)

    pressureIncrementedBoundaryConditions=>boundaryConditionsVariable%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
    IF(.NOT.ASSOCIATED(pressureIncrementedBoundaryConditions)) &
      & CALL FlagError("Pressure incremented boundary conditions is not associated for the boundary conditions variable.", &
      & err,error,*999)
       
    EXITS("BoundaryConditionsVariable_PressureIncConditionsGet")
    RETURN
999 NULLIFY(pressureIncrementedBoundaryConditions)
998 ERRORS("BoundaryConditionsVariable_PressureIncConditionsGet",err,error)
    EXITS("BoundaryConditionsVariable_PressureIncConditionsGet")
    RETURN 1
    
  END SUBROUTINE BoundaryConditionsVariable_PressureIncConditionsGet

  !
  !================================================================================================================================
  !

END MODULE BoundaryConditionAccessRoutines
