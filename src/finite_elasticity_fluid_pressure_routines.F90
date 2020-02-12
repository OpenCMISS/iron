!> \file
!> \authors Adam Reeve
!> \brief This module handles all routines pertaining to finite elasticity coupled with fluid pressure for poroelasticity.
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

!>This module handles all routines pertaining to finite elasticity coupled with fluid pressure for poroelasticity.


MODULE FINITE_ELASTICITY_FLUID_PRESSURE_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE EquationsRoutines
  USE EquationsSetAccessRoutines
  USE FINITE_ELASTICITY_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PUBLIC FinElasticityFluidPressure_EquationsSetSetup
  PUBLIC FinElasticityFluidPressure_EquationsSetSolnMethodSet
  PUBLIC FinElasticityFluidPressure_EquationsSetSpecificationSet

  PUBLIC ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP
  PUBLIC FinElasticityFluidPressure_ProblemSpecificationSet
  
  PUBLIC FinElasticityFluidPressure_FiniteElementCalculate

  PUBLIC ELASTICITY_FLUID_PRESSURE_PRE_SOLVE
  PUBLIC ELASTICITY_FLUID_PRESSURE_POST_SOLVE

  PUBLIC FinElasticityFluidPressure_ControlLoopPreLoop

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity fluid pressure equation type of a multi physics equations set class.
  SUBROUTINE FinElasticityFluidPressure_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("FinElasticityFluidPressure_EquationsSetSolnMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a "// &
          & "finite elasticity-fluid pressure class equations set.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a finite elasticity fluid pressure equation type of a multi physics equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("FinElasticityFluidPressure_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSolnMethodSet",ERR,ERROR)
    EXITS("FinElasticityFluidPressure_EquationsSetSolnMethodSet")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a finite elasticity fluid pressure equation type of a fluid mechanics equations set class.
  SUBROUTINE FinElasticityFluidPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FinElasticityFluidPressure_EquationsSetSpecificationSet",ERR,ERROR,*999)

    CALL FlagError("FinElasticityFluidPressure_EquationsSetSpecificationSet is not implemented.",ERR,ERROR,*999)

    EXITS("FinElasticityFluidPressure_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSpecificationSet",ERR,ERROR)
    EXITS("FinElasticityFluidPressure_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equation.
  SUBROUTINE FinElasticityFluidPressure_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string


    ENTERS("FinElasticityFluidPressure_EquationsSetSetup",ERR,ERROR,*999)

    CALL FlagError("FinElasticityFluidPressure_EquationsSetSetup is not implemented.",ERR,ERROR,*999)

    EXITS("FinElasticityFluidPressure_EquationsSetSetup")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSetup",ERR,ERROR)
    EXITS("FinElasticityFluidPressure_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity fluid pressure equation finite element equations set.
  SUBROUTINE FinElasticityFluidPressure_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    ENTERS("FinElasticityFluidPressure_FiniteElementCalculate",ERR,ERROR,*999)

    CALL FlagError("FinElasticityFluidPressure_FiniteElementCalculate is not implemented.",ERR,ERROR,*999)

    EXITS("FinElasticityFluidPressure_FiniteElementCalculate")
    RETURN
999 ERRORS("FinElasticityFluidPressure_FiniteElementCalculate",ERR,ERROR)
    EXITS("FinElasticityFluidPressure_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a finite elasticity fluid pressure equation type.
  SUBROUTINE FinElasticityFluidPressure_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("FinElasticityFluidPressure_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a finite elasticity fluid pressure type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE, &
          & problemSubtype]
      ELSE
        CALL FlagError("Finite elasticity fluid pressure problem specificaion must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("FinElasticityFluidPressure_ProblemSpecificationSet")
    RETURN
999 ERRORS("FinElasticityFluidPressure_ProblemSpecificationSet",err,error)
    EXITS("FinElasticityFluidPressure_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equations problem.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS
    TYPE(SolversType), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(ALLOCATED(problem%specification)) THEN
        IF(.NOT.ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(problem%specification,1)<3) THEN
          CALL FlagError("Problem specification must have three entries for a finite elasticity-Darcy problem.", &
            & err,error,*999)
        END IF
      ELSE
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))

      !--------------------------------------------------------------------
      !   Standard finite elasticity fluid pressure
      !--------------------------------------------------------------------
      CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%setupType)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for an finite elasticity ALE fluid pressure  equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a load increment loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,CONTROL_LOAD_INCREMENT_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity fluid pressure equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation for the solver
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL Solvers_NumberOfSolversSet(SOLVERS,1,ERR,ERROR,*999)
            !
            !Set the first solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
                & " is invalid for a finite elasticity fluid pressure equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop and solvers
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            ! 
            !Get the nonlinear solver and create the solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !
            !Finish the creation of the solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a finite elasticity fluid pressure equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
            & " is invalid for a finite elasticity ALE fluid pressure equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT

      CASE DEFAULT
        LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a standard finite elasticity fluid pressure equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)

      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_PROBLEM_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity fluid pressure problem pre-solve.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for an elasticity fluid pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            IF(CONTROL_LOOP%loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
              CALL FiniteElasticity_PreSolve(solver,err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a fluid pressure fluid type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("ELASTICITY_FLUID_PRESSURE_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_PRE_SOLVE
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure  problem post solve.
  SUBROUTINE ELASTICITY_FLUID_PRESSURE_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for an elasticity fluid pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            CALL FiniteElasticity_PostSolve(solver,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
              & " is not valid for a finite elasticity fluid pressure type of a multi physics problem class."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE")
    RETURN
999 ERRORSEXITS("ELASTICITY_FLUID_PRESSURE_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE ELASTICITY_FLUID_PRESSURE_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE FinElasticityFluidPressure_ControlLoopPreLoop(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string

    !Local Variables
    TYPE(SolverType), POINTER :: SOLVER_FLUID_PRESSURE
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP_FLUID_PRESSURE

    ENTERS("FinElasticityFluidPressure_ControlLoopPreLoop",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP_FLUID_PRESSURE)
    NULLIFY(SOLVER_FLUID_PRESSURE)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
        ! Eventually we may want to do different things depending on problem type/subtype
        ! too but for now we can just check the loop type.
        SELECT CASE(CONTROL_LOOP%loopType)
        CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          IF(CONTROL_LOOP%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"-- Starting load increment --",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(GENERAL_OUTPUT_TYPE,"LOAD INCREMENT NUMBER =           ", &
              & CONTROL_LOOP%loadIncrementLoop%iterationNumber,ERR,ERROR,*999)
            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
          ENDIF
          IF(DIAGNOSTICS1) THEN
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting load increment --",ERR,ERROR,*999)
            CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"LOAD INCREMENT NUMBER =           ", &
              & CONTROL_LOOP%loadIncrementLoop%iterationNumber,ERR,ERROR,*999)
            CALL WRITE_STRING(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",ERR,ERROR,*999)
          ENDIF

        CASE DEFAULT
          !do nothing
        END SELECT
      ELSE
        CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("FinElasticityFluidPressure_ControlLoopPreLoop")
    RETURN
999 ERRORS("FinElasticityFluidPressure_ControlLoopPreLoop",ERR,ERROR)
    EXITS("FinElasticityFluidPressure_ControlLoopPreLoop")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_ControlLoopPreLoop

  !
  !================================================================================================================================
  !

END MODULE FINITE_ELASTICITY_FLUID_PRESSURE_ROUTINES
