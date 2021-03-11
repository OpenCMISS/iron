!> \file
!> \authors 
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
!> Contributor(s): Adam Reeve, Chris Bradley
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
MODULE FiniteElasticityFluidPressureRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE EquationsRoutines
  USE EquationsSetAccessRoutines
  USE FiniteElasticityRoutines
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

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FinElasticityFluidPressure_EquationsSetSetup
  
  PUBLIC FinElasticityFluidPressure_EquationsSetSolnMethodSet
  
  PUBLIC FinElasticityFluidPressure_EquationsSetSpecificationSet

  PUBLIC FinElasticityFluidPressure_ProblemSetup
  
  PUBLIC FinElasticityFluidPressure_ProblemSpecificationSet
  
  PUBLIC FinElasticityFluidPressure_FiniteElementCalculate

  PUBLIC FinElasticityFluidPressure_PreSolve
  
  PUBLIC FinElasticityFluidPressure_PostSolve

  PUBLIC FinElasticityFluidPressure_PreLoop

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity fluid pressure equation type of a multi physics equations set class.
  SUBROUTINE FinElasticityFluidPressure_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FinElasticityFluidPressure_EquationsSetSolnMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a finite elasticity fluid pressure equation type of a multi physics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FinElasticityFluidPressure_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSolnMethodSet",err,error)
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

    ENTERS("FinElasticityFluidPressure_EquationsSetSpecificationSet",err,error,*999)

    CALL FlagError("FinElasticityFluidPressure_EquationsSetSpecificationSet is not implemented.",err,error,*999)

    EXITS("FinElasticityFluidPressure_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSpecificationSet",err,error)
    EXITS("FinElasticityFluidPressure_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure equation.
  SUBROUTINE FinElasticityFluidPressure_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string


    ENTERS("FinElasticityFluidPressure_EquationsSetSetup",err,error,*999)

    CALL FlagError("FinElasticityFluidPressure_EquationsSetSetup is not implemented.",err,error,*999)

    EXITS("FinElasticityFluidPressure_EquationsSetSetup")
    RETURN
999 ERRORS("FinElasticityFluidPressure_EquationsSetSetup",err,error)
    EXITS("FinElasticityFluidPressure_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity fluid pressure equation finite element equations set.
  SUBROUTINE FinElasticityFluidPressure_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("FinElasticityFluidPressure_FiniteElementCalculate",err,error,*999)

    CALL FlagError("FinElasticityFluidPressure_FiniteElementCalculate is not implemented.",err,error,*999)

    EXITS("FinElasticityFluidPressure_FiniteElementCalculate")
    RETURN
999 ERRORS("FinElasticityFluidPressure_FiniteElementCalculate",err,error)
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
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("FinElasticityFluidPressure_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error))// &
        & " is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a finite elasticity fluid pressure type of a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_FINITE_ELASTICITY_FLUID_PRESSURE_TYPE,problemSubtype]

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
  SUBROUTINE FinElasticityFluidPressure_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("FinElasticityFluidPressure_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " does not equal a standard finite elasticity fluid pressure equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for an finite elasticity ALE fluid pressure  equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a load increment loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRoot(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity fluid pressure equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRoot(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation for the solver
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !
        !Set the first solver to be a nonlinear solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity fluid pressure equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !Get the control loop and solvers
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRoot(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      !Get the nonlinear solver 
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)       
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the creation of the solver equations
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity fluid pressure equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a finite elasticity ALE fluid pressure equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("FinElasticityFluidPressure_ProblemSetup")
    RETURN
999 ERRORSEXITS("FinElasticityFluidPressure_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity fluid pressure problem pre-solve.
  SUBROUTINE FinElasticityFluidPressure_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopType,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FinElasticityFluidPressure_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      IF(loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE) CALL FiniteElasticity_PreSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a fluid pressure fluid type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FinElasticityFluidPressure_PreSolve")
    RETURN
999 ERRORSEXITS("FinElasticityFluidPressure_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity fluid pressure problem post solve.
  SUBROUTINE FinElasticityFluidPressure_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FinElasticityFluidPressure_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      CALL FiniteElasticity_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a finite elasticity fluid pressure type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FinElasticityFluidPressure_PostSolve")
    RETURN
999 ERRORSEXITS("FinElasticityFluidPressure_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_PostSolve

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE FinElasticityFluidPressure_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iterationNumber,loopType,outputType
    TYPE(ProblemType), POINTER :: problem

    ENTERS("FinElasticityFluidPressure_PreLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    
    ! Eventually we may want to do different things depending on problem type/subtype
    ! too but for now we can just check the loop type.
    SELECT CASE(loopType)
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
      CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"-- Starting load increment --",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Load Increment Number =           ",iterationNumber,err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting load increment --",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Load Increment Number =           ",iterationNumber,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF      
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("FinElasticityFluidPressure_PreLoop")
    RETURN
999 ERRORS("FinElasticityFluidPressure_PreLoop",err,error)
    EXITS("FinElasticityFluidPressure_PreLoop")
    RETURN 1
    
  END SUBROUTINE FinElasticityFluidPressure_PreLoop

  !
  !================================================================================================================================
  !

END MODULE FiniteElasticityFluidPressureRoutines
