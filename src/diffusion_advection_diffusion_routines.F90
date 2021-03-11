!> \file
!> \authors Andrew Cookson
!> \brief This module handles all routines pertaining to diffusion coupled to diffusion.
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
!> Contributor(s):  Andrew Cookson, Chris Bradley
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

!>This module handles all routines pertaining to diffusion coupled to diffusion.
MODULE DiffusionAdvectionDiffusionRoutines

  USE AdvectionDiffusionEquationsRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines  
  USE DiffusionEquationsRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths  
  USE MatrixVector
  USE MeshRoutines
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"


  IMPLICIT NONE

  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSetup
  
  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSpecSet
  
  PUBLIC DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet

  PUBLIC DiffusionAdvectionDiffusion_ProblemSetup
  
  PUBLIC DiffusionAdvectionDiffusion_ProblemSpecificationSet
  
  PUBLIC DiffusionAdvectionDiffusion_FiniteElementCalculate

  PUBLIC DiffusionAdvectionDiffusion_PreSolve
  
  PUBLIC DiffusionAdvectionDiffusion_PostSolve

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
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
        & " is not valid for a diffusion & advection-diffusion equation type of a multi physics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet",err,error)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion & advection-diffusion coupled equation.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSetup",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)
             
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSetup",err,error)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSetup")
    RETURN 1

  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a coupled diffusion & advection-diffusion equation finite element equations set.
  SUBROUTINE DiffusionAdvectionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("DiffusionAdvectionDiffusion_FiniteElementCalculate",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)
      
    EXITS("DiffusionAdvectionDiffusion_FiniteElementCalculate")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_FiniteElementCalculate",err,error)
    EXITS("DiffusionAdvectionDiffusion_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSpecSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("DiffusionAdvectionDiffusion_EquationsSetSpecSet",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)

    EXITS("DiffusionAdvectionDiffusion_EquationsSetSpecSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_EquationsSetSpecSet",err,error)
    EXITS("DiffusionAdvectionDiffusion_EquationsSetSpecSet")
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_EquationsSetSpecSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a coupled diffusion & advection-diffusion problem.
  SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a coupled diffusion & advection-diffusion type of a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_DIFFUSION_ADVECTION_DIFFUSION_TYPE,problemSubtype]

    EXITS("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_ProblemSpecificationSet",err,error)
    EXITS("DiffusionAdvectionDiffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the coupled diffusion-diffusion equations problem.
  SUBROUTINE DiffusionAdvectionDiffusion_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solverDiffusion,solverAdvectionDiffusion
    TYPE(SolverEquationsType), POINTER :: solverEquationsDiffusion,solverEquationsAdvectionDiffusion
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DiffusionAdvectionDiffusion_ProblemSetup",err,error,*999)

    CALL Problem_Specification(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " does not equal a coupled source diffusion & advection-diffusion equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a diffusion & advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a time control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ContolLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a coupled diffusion & advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
        !Set the first solver to be a linear solver for the advection-diffusion problem
        NULLIFY(solverAdvectionDiffusion)
        CALL Solvers_SolverGet(solvers,1,solverAdvectionDiffusion,err,error,*999)
        CALL Solver_TypeSet(solverAdvectionDiffusion,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solverAdvectionDiffusion,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solverAdvectionDiffusion,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solverAdvectionDiffusion,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_LibraryTypeSet(solverAdvectionDiffusion,SOLVER_CMISS_LIBRARY,err,error,*999)
        !Set the second solver to be a linear solver for the diffusion problem
        NULLIFY(solverDiffusion)
        CALL Solvers_SolverGet(solvers,2,solverDiffusion,err,error,*999)
        CALL Solver_TypeSet(solverDiffusion,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solverDiffusion,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solverDiffusion,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solverDiffusion,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_LibraryTypeSet(solverDiffusion,SOLVER_CMISS_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a coupled diffusion & advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !Get the control loop and solvers
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoop)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the advection-diffusion solver and create the advection-diffusion solver equations
        CALL Solvers_SolverGet(solvers,1,solverAdvectionDiffusion,err,error,*999)
        CALL SolverEquations_CreateStart(solverAdvectionDiffusion,solverEquationsAdvectionDiffusion,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquationsAdvectionDiffusion,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquationsAdvectionDiffusion, & 
          & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquationsAdvectionDiffusion,SOLVER_SPARSE_MATRICES,err,error,*999)
        !Get the diffusion solver and create the diffusion solver equations
        CALL Solvers_SolverGet(solvers,2,solverDiffusion,err,error,*999)
        CALL SolverEquations_CreateStart(solverDiffusion,solverEquationsDiffusion,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquationsDiffusion,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquationsDiffusion, & 
          & SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquationsDiffusion,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the creation of the advection-diffusion solver equations
        NULLIFY(solverAdvectionDiffusion)
        CALL Solvers_SolverGet(solvers,1,solverAdvectionDiffusion,err,error,*999)
        NULLIFY(solverEquationsAdvectionDiffusion)
        CALL Solver_SolverEquationsGet(solverAdvectionDiffusion,solverEquationsAdvectionDiffusion,err,error,*999)
        CALL SolverEquations_CreateFinish(solverEquationsAdvectionDiffusion,err,error,*999)             
        !
        !Finish the creation of the diffusion solver equations
        NULLIFY(solverDiffusion)
        CALL Solvers_SolverGet(solvers,2,solverDiffusion,err,error,*999)
        NULLIFY(solverEquationsDiffusion)
        CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
        CALL SolverEquations_CreateFinish(solverEquationsDiffusion,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a coupled diffusion & advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a coupled diffusion & advection-diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("DiffusionAdvectionDiffusion_ProblemSetup")
    RETURN
999 ERRORSEXITS("DiffusionAdvectionDiffusion_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion-diffusion problem pre-solve.
  SUBROUTINE DiffusionAdvectionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("DiffusionAdvectionDiffusion_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      IF(solverNumber==1) THEN
        !copy current value of concentration_one to another variable
        CALL AdvectionDiffusion_PreSolveStoreCurrentSoln(solver,err,error,*999)
        !Set source term to be updated value of concentration_two
        !CALL AdvectionDiffusion_PreSolveGetSourceValue(solver,err,error,*999)        
        !Update indpendent data fields
        !CALL AdvectionDiffusion_PreSolveUpdateInputData(solver,err,error,*999) 
      ELSE IF(solverNumber==2) THEN
        !copy current value of concentration_one to another variable
        CALL Diffusion_PreSolveStoreCurrentSolution(solver,err,error,*999)
        !compute value of constant source term - evaluated from lamdba*(0.5*(c_1^{t+1}+c_1^{t}) - c_2^{t})
        !CALL Diffusion_PreSolveGetSourceValue(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DiffusionAdvectionDiffusion_PreSolve")
    RETURN
999 ERRORSEXITS("DiffusionAdvectionDiffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Sets up the diffusion-diffusion problem post solve.
  SUBROUTINE DiffusionAdvectionDiffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("DiffusionAdvectionDiffusion_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_Specification(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      IF(solverNumber==1) THEN
        !CALL AdvectionDiffusion_PostSolve(solver,err,error,*999)
      ELSE IF(solverNumber==2) THEN
        !CALL Diffusion_PostSolve(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DiffusionAdvectionDiffusion_PostSolve")
    RETURN
999 ERRORSEXITS("DiffusionAdvectionDiffusion_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE DiffusionAdvectionDiffusion_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("DiffusionAdvectionDiffusion_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      !CALL AdvectionDiffusion_PostSolveOutputData(controlLoop,solver,err,error,*999)
      !CALL Diffusion_PostSolveOutputData(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(controlLoop%problem%SPECIFICATION(3),"*",err,error))// &
        & " is not valid for a diffusion & advection-diffusion type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("DiffusionAdvectionDiffusion_PostSolveOutputData")
    RETURN
999 ERRORS("DiffusionAdvectionDiffusion_PostSolveOutputData",err,error)
    EXITS("DiffusionAdvectionDiffusion_PostSolveOutputData")
    RETURN 1
    
  END SUBROUTINE DiffusionAdvectionDiffusion_PostSolveOutputData
      
  !   
  !================================================================================================================================
  !

END MODULE DiffusionAdvectionDiffusionRoutines
