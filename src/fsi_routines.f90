!> \file
!> \authors Andreas Hessenthaler
!> \brief This module handles all routines pertaining to finite elasticity coupled with navier stokes for fsi problems
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

!>This module handles all routines pertaining to finite elasticity coupled with navier stokes for fsi problems
MODULE FSIRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE CMISS_CELLML
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE NAVIER_STOKES_EQUATIONS_ROUTINES
  USE ProblemAccessRoutines
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PUBLIC FSI_ProblemSetup
  
  PUBLIC FSI_ProblemSpecificationSet

  PUBLIC FSI_PreSolve
  PUBLIC FSI_PostSolve

  PUBLIC FSI_PreLoop
  PUBLIC FSI_PostLoop

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a finite elasticity Navier-Stokes equation type.
  SUBROUTINE FSI_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("FSI_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) &
      & CALL FlagError("Finite elasticity Navier-Stokes problem specificaion must have a least three entries.",err,error,*999)
    
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
      ALLOCATE(problem%specification(3),stat=err)
      IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
      problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS, &
        & PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE, &
        & problemSubtype]
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a finite elasticity Navier-Stokes type of a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FSI_ProblemSpecificationSet")
    RETURN
999 ERRORS("FSI_ProblemSpecificationSet",err,error)
    EXITS("FSI_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FSI_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity navier stokes equations problem.
  SUBROUTINE FSI_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("FSI_ProblemSetup",err,error,*999)

    IF(.NOT.ASSOCIATED(PROBLEM)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have >= 3 entries for a finite elasticity-ALE NS problem.",err,error,*999)
        
    SELECT CASE(problem%specification(3))       
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
      !Finite Elasticity with Navier Stokes ALE
      SELECT CASE(problemSetup%SETUP_TYPE)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(problemSetup%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an finite elasticity ALE navier stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop as parent loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_OutputTypeSet(controlLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)       
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          NULLIFY(controlLoop)
          CALL Problem_ControlLoopGet(problem,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity navier stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        !Get the control loop
        NULLIFY(controlLoop)
        CALL Problem_ControlLoopGet(problem,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        SELECT CASE(problemSetup%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Start the solvers creation
          NULLIFY(solvers)
          CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
          SELECT CASE(problem%specification(3))       
          CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            CALL Solvers_NumberSet(solvers,3,err,error,*999)
            !Set the first solver to be an CellML Evaluator for time varying boundary conditions
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI boundary condition CellML evaluation solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a first order dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI dynamic nonlinear solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be a linear solver for the Laplace mesh movement problem
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI mesh movement linear solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            CALL Solvers_NumberSet(solvers,5,err,error,*999)
            !Set the first solver to be an CellML Evaluator for time varying boundary conditions
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI fluid boundary condition CellML evaluation solver",err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI solid boundary condition CellML evaluation solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the third solver to be an CellML Integrator for growth rate integration
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI growth CellML integration solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the fourth solver to be a first order dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI dynamic nonlinear solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the fifth solver to be a linear solver for the Laplace mesh movement problem
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,5,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"FSI mesh movement linear solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
              & " does not equal a standard finite elasticity navier stokes equation subtype."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity navier stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoop)
          CALL Problem_ControlLoopGet(problem,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          SELECT CASE(problem%specification(3))       
          CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            !Get the dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            !Create the solver equations
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !Get the linear moving mesh solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            !Create the solver equations
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            !Get the dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            !Create the solver equations
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !Get the linear moving mesh solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,5,solver,err,error,*999)
            !Create the solver equations
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
              & " does not equal a standard finite elasticity navier stokes equation subtype."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoop)
          CALL Problem_ControlLoopGet(problem,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          SELECT CASE(problem%specification(3))       
          CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            !Get the dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the dynamic solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)    
            !Get the moving mesh solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the moving mesh solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
          CASE(PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
            & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
            !Get the dynamic solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,4,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the dynamic solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)    
            !Get the moving mesh solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,5,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the moving mesh solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
          CASE DEFAULT
            localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
              & " does not equal a standard finite elasticity navier stokes equation subtype."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity navier stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        !Get the control loop
        NULLIFY(controlLoop)
        CALL Problem_ControlLoopGet(problem,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(problemSetup%ACTION_TYPE)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the fluid CellML BC evaluation solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the CellML equations
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          !Get the solid CellML BC evaluation solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Create the CellML equations
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          IF(problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            !Get the CellML growth integration solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            !Create the CellML equations
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the fluid CellML BC evaluator solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Get the CellML equations for the CellML evaluator solver
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          !Get the solid CellML BC evaluator solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Get the CellML equations for the CellML evaluator solver
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          IF(problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            !Get the CellML growth integration solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,3,solver,err,error,*999)
            !Get the CellML equations for the CellML integration solver
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ENDIF
        CASE DEFAULT            
          localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a finite elasticity equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
          & " is invalid for a finite elasticity ALE navier stokes equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " does not equal a standard finite elasticity navier stokes equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("FSI_ProblemSetup")
    RETURN
999 ERRORSEXITS("FSI_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FSI_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity navier stokes problem pre-solve.
  SUBROUTINE FSI_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FSI_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      CALL FlagError("Problem specification must have three entries for an elasticity Navier-Stokes problem.",err,error,*999)
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
      IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE.OR.solver%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
        IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          CALL NavierStokes_PreSolveALEUpdateMesh(solver,err,error,*999)
        ENDIF
        !Pre solve for ALE NavierStokes equations set
        CALL NAVIER_STOKES_PRE_SOLVE(solver,err,error,*999)
        !Pre solve for FiniteElasticity equations set
        !Nothing to be done???
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("FSI_PreSolve")
    RETURN
999 ERRORSEXITS("FSI_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FSI_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Performs finite elasticity navier stokes problem post solve.
  SUBROUTINE FSI_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: Err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: Error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    LOGICAL :: fluidEquationsSetFound
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_TYPE), POINTER :: dynamicSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("FSI_PostSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      CALL FlagError("Problem specification must have three entries for an elasticity Navier-Stokes problem.",err,error,*999)
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
      & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, & 
      & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      IF(solver%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
        !Post solve for the linear solver
        NULLIFY(dynamicSolver)
        IF(problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
          & problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
          CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
        ELSE
          CALL Solvers_SolverGet(solvers,4,dynamicSolver,err,error,*999)
        ENDIF
        IF(ASSOCIATED(dynamicSolver%DYNAMIC_SOLVER)) THEN
          dynamicSolver%DYNAMIC_SOLVER%ALE=.TRUE.
        ELSE  
          CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
        END IF
      ELSE IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
        !Post solve for the dynamic solver
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        equationsSetIdx=1
        fluidEquationsSetFound=.FALSE.
        DO WHILE(.NOT.fluidEquationsSetFound.AND.equationsSetIdx<=solverMapping%NUMBER_OF_EQUATIONS_SETS)
          equations=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equations
          NULLIFY(equationsSet)
          CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
          IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
            & .AND.equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
            & .AND.(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
            & .OR.equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
            fluidEquationsSetFound=.TRUE.
          ELSE
            equationsSetIdx=equationsSetIdx+1
          END IF
        END DO
        IF(.NOT.fluidEquationsSetFound) THEN
          localError="Fluid equations set not found."
          CALL FlagError(localError,Err,Error,*999)
        END IF
        !CALL NavierStokes_WallShearStressCalculate(equationsSet,err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FSI_PostSolve")
    RETURN
999 ERRORSEXITS("FSI_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FSI_PostSolve

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE FSI_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables

    ENTERS("FSI_PreLoop",err,error,*999)
    
    CALL FlagError("FSI_PreLoop not implemented.",err,error,*999)

    EXITS("FSI_PreLoop")
    RETURN
999 ERRORSEXITS("FSI_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FSI_PreLoop

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration. Updates interface and fluid geometric fields and exports fields.
  SUBROUTINE FSI_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,currentIteration,derivativeIdx,equationsSetIndex,nodeIdx,numberOfComponents,outputIteration, &
      & solidNode,versionIdx
    LOGICAL :: fluidEquationsSetFound=.FALSE.,solidEquationsSetFound=.FALSE.
    REAL(DP) :: startTime,currentTime,stopTime,timeIncrement,VALUE
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(EQUATIONS_SET_TYPE), POINTER :: solidEquationsSet,fluidEquationsSet,equationsSet
    TYPE(FIELD_TYPE), POINTER :: solidGeometricField,interfaceGeometricField,solidDependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: geometricVariable
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: interfaceCondition
    TYPE(INTERFACE_TYPE), POINTER :: fsInterface
    TYPE(NODES_TYPE), POINTER :: interfaceNodes
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_TYPE), POINTER :: dynamicSolver,linearSolver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: dynamicSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: dynamicSolverMapping
    TYPE(VARYING_STRING) :: method,solidFileName,fluidFileName,interfaceFileName
 
    ENTERS("FSI_PostLoop",err,error,*999)
    
    !Check pointers
    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop not associated.",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    NULLIFY(solvers)
    CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
    !Get times
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,err,error,*999)    
    !Get solvers for FSI
    IF(problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
      & problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
      NULLIFY(dynamicSolver)
      CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
      NULLIFY(linearSolver)
      CALL Solvers_SolverGet(solvers,3,linearSolver,err,error,*999)
    ELSE IF(problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
      & problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
      NULLIFY(dynamicSolver)
      CALL Solvers_SolverGet(solvers,4,dynamicSolver,err,error,*999)
      NULLIFY(linearSolver)
      CALL Solvers_SolverGet(solvers,5,linearSolver,err,error,*999)
    ELSE
      NULLIFY(dynamicSolver)
      CALL Solvers_SolverGet(solvers,3,dynamicSolver,err,error,*999)
      NULLIFY(linearSolver)
      CALL Solvers_SolverGet(solvers,4,linearSolver,err,error,*999)
    ENDIF
    !==============================================================================================================================
    !First update mesh and calculate boundary velocity values
    !CALL NavierStokes_PreSolveALEUpdateMesh(dynamicSolver,err,error,*999)
    !==============================================================================================================================
    
    !Update interface geometric field and export results
    NULLIFY(dynamicSolverEquations)
    CALL Solver_SolverEquationsGet(dynamicSolver,dynamicSolverEquations,err,error,*999)
    NULLIFY(dynamicSolverMapping)
    CALL SolverEquations_SolverMappingGet(dynamicSolverEquations,dynamicSolverMapping,err,error,*999)
    
    equationsSetIndex=1
    fluidEquationsSetFound=.FALSE.
    solidEquationsSetFound=.FALSE.
    DO WHILE((equationsSetIndex<=dynamicSolverMapping%NUMBER_OF_EQUATIONS_SETS.AND..NOT.solidEquationsSetFound) &
      & .OR.(equationsSetIndex<=dynamicSolverMapping%NUMBER_OF_EQUATIONS_SETS.AND..NOT.fluidEquationsSetFound))
      equationsSet=>dynamicSolverMapping%EQUATIONS_SETS(equationsSetIndex)%ptr
      IF(equationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS &
        & .AND.equationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE &
        & .AND.((equationsSet%specification(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE).OR. &
        & (equationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE).OR. &
        & (equationsSet%specification(3)==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE))) THEN
        solidEquationsSet=>equationsSet
        solidEquationsSetFound=.TRUE.
      ELSE IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
        & .AND.equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
        & .AND.(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
        & .OR.equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
        fluidEquationsSet=>equationsSet
        fluidEquationsSetFound=.TRUE.
      ELSE
        CALL FlagError("Invalid equations sets associated with dynamic solver for FSI.",err,error,*999)
      ENDIF
      equationsSetIndex=equationsSetIndex+1
    ENDDO
    IF(.NOT.solidEquationsSetFound) CALL FlagError("Could not find solid equations set for FSI.",err,error,*999)
    IF(.NOT.fluidEquationsSetFound) CALL FlagError("Could not find fluid equations set for FSI.",err,error,*999)
    NULLIFY(solidGeometricField)
    CALL EquationsSet_GeometricFieldGet(solidEquationsSet,solidGeometricField,err,error,*999)
    CALL Field_NumberOfComponentsGet(solidGeometricField,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
    IF(DynamicSolverMapping%NUMBER_OF_INTERFACE_CONDITIONS>1) &
      & CALL FlagError("Invalid number of interface conditions. Must be 1 for FSI.",err,error,*999)
    NULLIFY(solidDependentField)
    CALL EquationsSet_DependentFieldGet(solidEquationsSet,solidDependentField,err,error,*999)
    NULLIFY(interfaceCondition)
    CALL SolverMapping_InterfaceConditionGet(dynamicSolverMapping,1,interfaceCondition,err,error,*999)
    NULLIFY(fsInterface)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,fsInterface,err,error,*999)
    NULLIFY(interfaceNodes)
    CALL Interface_NodesGet(fsInterface,interfaceNodes,err,error,*999)
    NULLIFY(interfaceGeometricField)
    CALL InterfaceCondition_GeometricFieldGet(interfaceCondition,interfaceGeometricField,err,error,*999)
    !===============================================================================================================
    !Update interface geometric field

!!TODO: Make this a interface subroutine that will update the interface and do it all properly
    
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    DO componentIdx=1,numberOfComponents
      SELECT CASE(geometricVariable%components(componentIdx)%INTERPOLATION_TYPE)
      CASE(FIELD_NODE_BASED_INTERPOLATION)
        NULLIFY(domain)
        CALL FieldVariable_DomainGet(geometricVariable,componentIdx,domain,err,error,*999)
        DO nodeIdx=1,domain%topology%nodes%TOTAL_NUMBER_OF_NODES
          solidNode=interfaceNodes%COUPLED_NODES(1,nodeIdx)
          DO derivativeIdx=1,domain%topology%nodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
            DO versionIdx=1,domain%topology%nodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
              CALL Field_ParameterSetGetNode(solidDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,solidNode,componentIdx,VALUE,err,error,*999)
              CALL Field_ParameterSetUpdateNode(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,&
                & versionIdx,derivativeIdx,nodeIdx,componentIdx,VALUE,err,error,*999)
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
        
      CASE DEFAULT
        CALL FlagError("Interface geometric component does not have node based interpolation.",err,error,*999)
      END SELECT
    ENDDO !componentIdx
    CALL Field_ParameterSetUpdateStart(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateFinish(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    
    !===============================================================================================================
    
    IF(outputIteration/=0) THEN
      IF(MOD(currentIteration,outputIteration)==0) THEN
        !Export fields
        solidFileName="./output/Solid/Solid"//TRIM(NumberToVString(INT(currentIteration),"*",err,error))
        fluidFileName="./output/Fluid/Fluid"//TRIM(NumberToVString(INT(currentIteration),"*",err,error))
        interfaceFileName="./output/Interface/Interface"//TRIM(NumberToVString(INT(currentIteration),"*",err,error))
        method="FORTRAN"
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Export fields... ",err,error,*999)
        ENDIF
        !Export solid fields
        IF(.NOT.ASSOCIATED(solidEquationsSet%region)) CALL FlagError("Solid region not associated.",err,error,*999)
        IF(.NOT.ASSOCIATED(solidEquationsSet%region%fields)) CALL FlagError("Solid fields not associated.",err,error,*999)
        CALL FIELD_IO_NODES_EXPORT(solidEquationsSet%region%fields,solidFileName,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(solidEquationsSet%region%fields,solidFileName,method,err,error,*999)
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,solidFileName,err,error,*999)
        IF(.NOT.ASSOCIATED(fluidEquationsSet%REGION)) CALL FlagError("Fluid region not associated.", err,error,*999)
        IF(.NOT.ASSOCIATED(fluidEquationsSet%REGION%FIELDS)) CALL FlagError("Fluid fields not associated.",err,error,*999)
        !Export fluid fields
        IF(.NOT.ASSOCIATED(fluidEquationsSet%region)) CALL FlagError("Fluid region not associated.",err,error,*999)
        IF(.NOT.ASSOCIATED(fluidEquationsSet%region%fields)) CALL FlagError("Fluid fields not associated.",err,error,*999)
        CALL FIELD_IO_NODES_EXPORT(fluidEquationsSet%region%fields,fluidFileName,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(fluidEquationsSet%region%fields,fluidFileName,method,err,error,*999)
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,fluidFileName,err,error,*999)
        !Export interface fields
        IF(.NOT.ASSOCIATED(fsInterface%fields)) CALL FlagError("Interface fields not associated.",err,error,*999)
        CALL FIELD_IO_NODES_EXPORT(fsInterface%fields,interfaceFileName,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(fsInterface%fields,interfaceFileName,method,err,error,*999)
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,interfaceFileName,err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
        ENDIF
      ENDIF
    ENDIF
        
    EXITS("FSI_PostLoop")
    RETURN
999 ERRORSEXITS("FSI_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FSI_PostLoop

  !
  !================================================================================================================================
  !

END MODULE FSIRoutines
