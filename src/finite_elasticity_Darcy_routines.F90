!> \file
!> \authors Christian Michler, Jack Lee
!> \brief This module handles all routines pertaining to finite elasticity coupled with Darcy.
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
!> Contributor(s): Christian Michler, Jack Lee, Chris Bradley
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

!>This module handles all routines pertaining to finite elasticity coupled with Darcy.
MODULE FiniteElasticityDarcyRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines  
  USE DarcyEquationsRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FINITE_ELASTICITY_ROUTINES
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

  PUBLIC FiniteElasticityDarcy_EquationsSetSetup
  
  PUBLIC FiniteElasticityDarcy_EquationsSetSpecificationSet
  
  PUBLIC FiniteElasticityDarcy_EquationsSetSolutionMethodSet

  PUBLIC FiniteElasticityDarcy_ProblemSetup
  
  PUBLIC FiniteElasticityDarcy_ProblemSpecificationSet
  
  PUBLIC FiniteElasticityDarcy_FiniteElementCalculate

  PUBLIC FiniteElasticityDarcy_PreSolve
  
  PUBLIC FiniteElasticityDarcy_PostSolve

  PUBLIC FiniteElasticityDarcy_PreLoop
  
  PUBLIC FiniteElasticityDarcy_PostLoop

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a finite elasticity Darcy equation type of a multi physics equations set class.
  SUBROUTINE FiniteElasticityDarcy_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FiniteElasticityDarcy_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_ELASTICITY_DARCY_SUBTYPE)
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
        & " is not valid for a finite elasticity Darcy equation type of a multi physics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("FiniteElasticityDarcy_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("FiniteElasticityDarcy_EquationsSetSolutionMethodSet",err,error)
    EXITS("FiniteElasticityDarcy_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy equation.
  SUBROUTINE FiniteElasticityDarcy_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FiniteElasticityDarcy_EquationsSetSetup",err,error,*999)

    CALL FlagError("FiniteElasticityDarcy_EquationsSetSetup still needs to be implemented.",err,error,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to setup the equations set of a monolithic
    ! finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since equationsSetSetup of respective equations_set is called.
    !=================================================================
             
    EXITS("FiniteElasticityDarcy_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a finite elasticity Darcy equation finite element equations set.
  SUBROUTINE FiniteElasticityDarcy_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FiniteElasticityDarcy_FiniteElementCalculate",err,error,*999)

    CALL FlagError("FiniteElasticityDarcy_FiniteElementCalculate still needs to be implemented.",err,error,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to calculate the finite-element matrices and vectors
    ! of a monolithic finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since FiniteElementCalculate of respective equations set is called.
    !=================================================================

    EXITS("FiniteElasticityDarcy_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a finite elasticity Darcy  equation type of a multi physics equations set class.
  SUBROUTINE FiniteElasticityDarcy_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("FiniteElasticityDarcy_EquationsSetSpecificationSet",err,error,*999)

    CALL FlagError("FiniteElasticityDarcy_EquationsSetSpecificationSet still needs to be implemented.",err,error,*999)

    !=================================================================
    ! This routine still needs to be implemented.
    ! It will be used to set the equations_set_subtype
    ! of a monolithic finite-elasticity Darcy system.
    ! For the partitioned solution this routine is not called,
    ! since EQUATIONS_SET_SUBTYPE_SET of respective equations_set is called.
    !=================================================================

    EXITS("FiniteElasticityDarcy_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("FiniteElasticityDarcy_EquationsSetSpecificationSet",err,error)
    EXITS("FiniteElasticityDarcy_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a finite elasticity Darcy equation type.
  SUBROUTINE FiniteElasticityDarcy_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a finite elasticity Darcy type of a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_FINITE_ELASTICITY_DARCY_TYPE,problemSubtype]

    EXITS("FiniteElasticityDarcy_ProblemSpecificationSet")
    RETURN
999 ERRORS("FiniteElasticityDarcy_ProblemSpecificationSet",err,error)
    EXITS("FiniteElasticityDarcy_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy equations problem.
  SUBROUTINE FiniteElasticityDarcy_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot,solidSubLoop,fluidSubLoop,subiterationLoop
    TYPE(SolverType), POINTER :: solver, solverMatProperties, solverSolid
    TYPE(SolverEquationsType), POINTER :: solverEquations, solverEquationsMatProperties, solverEquationsSolid
    TYPE(SolversType), POINTER :: solidSolvers,fluidSolvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
      !OK
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      !OK
    CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecifiction(3),"*",err,error))// &
        & " does not equal a standard finite elasticity Darcy equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE)
      !--------------------------------------------------------------------
      !   s t a n d a r d   f i n i t e   e l a s t i c i t y   D a r c y
      !--------------------------------------------------------------------
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
            & " is invalid for an finite elasticity ALE Darcy  equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,2,err,error,*999)
          !Solid, load incremented control loop
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,solidSubLoop,err,error,*999)
          CALL ControlLoop_TypeSet(solidSubLoop,CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
          !Fluid control loop
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,fluidSubLoop,err,error,*999)
          CALL ControlLoop_TypeSet(fluidSubLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
          !Sub-loops are finished when parent is finished
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
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
          !Start the solvers creation for the solid solver
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL Solvers_CreateStart(solidSubLoop,solidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(solidSolvers,1,err,error,*999)
          !
          !Set the first solver to be a nonlinear solver for the finite elasticity
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          CALL Solver_TypeSet(solverSolid,SOLVER_NONLINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solverSolid,SOLVER_PETSC_LIBRARY,err,error,*999)
          !
          !Start the solvers creation for the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL Solvers_CreateStart(fluidSubLoop,fluidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(fluidSolvers,2,err,error,*999)
          !
          !Set the second solver to be a linear solver for the material update
          NULLIFY(solverMatProperties)
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          CALL Solver_TypeSet(solverMatProperties,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solverMatProperties,SOLVER_PETSC_LIBRARY,err,error,*999)
          !
          !Set the third solver to be a linear solver for ALE Darcy
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solid solvers
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solidSolvers,err,error,*999)
          !Get the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(fluidSolvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop and solvers
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Get the finite elasticity solver and create the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL SolverEquations_CreateStart(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquationsSolid,SOLVER_SPARSE_MATRICES,err,error,*999)
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Get the material-properties solver and create the material-properties solver equations
          NULLIFY(solverMatProperties)
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          NULLIFY(solverEquationsMatProperties)
          CALL SolverEquations_CreateStart(solverMatProperties,solverEquationsMatProperties,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquationsMatProperties,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquationsMatProperties,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquationsMatProperties,SOLVER_SPARSE_MATRICES,err,error,*999)
          !
          !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Finish the creation of the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL Solver_SolverEquationsGet(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquationsSolid,err,error,*999)             
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Finish the creation of the material-properties solver equations
          NULLIFY(solverMatProperties)
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          NULLIFY(solverEquationsMatProperties)
          CALL Solver_SolverEquationsGet(solverMatProperties,solverEquationsMatProperties,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquationsMatProperties,err,error,*999)             
          !
          !Finish the creation of the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity ALE Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT      
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      !--------------------------------------------------------------------
      !   q u a s i s t a t i c   f i n i t e   e l a s t i c i t y   t r a n s i e n t   D a r c y
      !--------------------------------------------------------------------
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
            & " is invalid for an finite elasticity ALE Darcy  equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(controlLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up a subiteration loop
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          CALL ControlLoop_LabelSet(subiterationLoop,"Subiteration Loop",err,error,*999)
          CALL ControlLoop_TypeSet(subiterationLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_MaximumIterationsSet(subiterationLoop,9,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(subiterationLoop,2,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(subiterationLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up load incremented control loop for Solid
          nullify(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          CALL ControlLoop_LabelSet(solidSubLoop,"Finite Elasticity Load Increment Loop",err,error,*999)
          CALL ControlLoop_TypeSet(solidSubLoop,CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
          !For problems that require it, the user can get the solid subloop using:
          !CALL CMISSProblemControlLoopGet(Problem,[1,1,CMISSControlLoopNode],ControlLoopSolid,Err)
          !And then set the number of load increments to 3 for example with:
          !CALL CMISSControlLoopMaximumIterationsSet(ControlLoopSolid,3,Err)
          CALL ControlLoop_MaximumIterationsSet(solidSubLoop,1,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(solidSubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up control loop for Fluid
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          CALL ControlLoop_LabelSet(fluidSubLoop,"Darcy Simple Loop",err,error,*999)
          CALL ControlLoop_TypeSet(fluidSubLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(fluidSubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          controlLoopRoot=>PROBLEM%controlLoop
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
          !Sub-loops are finished when parent is finished
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
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
          !Start the solvers creation for the solid solver
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL Solvers_CreateStart(solidSubLoop,solidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(solidSolvers,1,err,error,*999)
          !
          !Set the solid solver to be a nonlinear solver for the finite elasticity
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          CALL Solver_TypeSet(solverSolid,SOLVER_NONLINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solverSolid,SOLVER_PETSC_LIBRARY,err,error,*999)
          !
          !Start the solvers creation for the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL Solvers_CreateStart(fluidSubLoop,fluidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(fluidSolvers,1,err,error,*999)
          !
          !Set the solver to be a first-order dynamic solver for Darcy
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solid solvers
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solidSolvers,err,error,*999)
          
          !Get the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(fluidSolvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop and solvers
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Get the finite elasticity solver and create the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL SolverEquations_CreateStart(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquationsSolid,SOLVER_SPARSE_MATRICES,err,error,*999)
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Finish the creation of the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL Solver_SolverEquationsGet(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquationsSolid,err,error,*999)             
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Finish the creation of the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity ALE Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !--------------------------------------------------------------------
      !   q u a s i s t a t i c   e l a s t i c i t y   t r a n s i e n t   D a r c y   M A T E R I A L   S O L V E
      !--------------------------------------------------------------------
      SELECT CASE(problemSetup%setupType)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing???
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for an finite elasticity ALE Darcy  equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(controlLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up a subiteration loop
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          CALL ControlLoop_TypeSet(subiterationLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_MaximumIterationsSet(subiterationLoop,9,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(subiterationLoop,2,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(subiterationLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up load incremented control loop for Solid
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          CALL ControlLoop_TypeSet(solidSubLoop,CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
          !For problems that require it, the user can get the solid subloop using:
          !CALL cmfe_Problem_ControlLoopGet(problem,[1,1,CMFE_CONTROL_LOOP_NODE],controlLoopSolid,err)
          !And then set the number of load increments to 3 for example with:
          !CALL cmfe_ControlLoop_MaximumIterationsSet(controlLoopSolid,3,err)
          CALL ControlLoop_MaximumIterationsSet(solidSubLoop,1,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(solidSubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
          
          !Set up control loop for Fluid
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          CALL ControlLoop_TypeSet(fluidSubLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
!!TODO: DEFAULT SHOULD BE NO OUPUT???
          CALL ControlLoop_OutputTypeSet(fluidSubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          controlLoopRoot=>problem%controlLoop
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
          !Sub-loops are finished when parent is finished
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
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
          !Start the solvers creation for the solid solver
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL Solvers_CreateStart(solidSubLoop,solidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(solidSolvers,1,err,error,*999)
          !
          !Set the solid solver to be a nonlinear solver for the finite elasticity
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          CALL Solver_TypeSet(solverSolid,SOLVER_NONLINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solverSolid,SOLVER_PETSC_LIBRARY,err,error,*999)
          !
          !Start the solvers creation for the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL Solvers_CreateStart(fluidSubLoop,fluidSolvers,err,error,*999)
          CALL Solvers_NumberOfSolversSet(fluidSolvers,2,err,error,*999)
          !
          !Set the solver to be a linear solver for the material update
          NULLIFY(solverMatProperties)
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          CALL Solver_TypeSet(solverMatProperties,SOLVER_LINEAR_TYPE,err,error,*999)
          CALL Solver_LibraryTypeSet(solverMatProperties,SOLVER_PETSC_LIBRARY,err,error,*999)
          !
          !Set the other solver to be a first-order dynamic solver for Darcy
          nullify(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solid solvers
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solidSolvers,err,error,*999)
          
          !Get the fluid solvers
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(fluidSolvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop and solvers
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Get the finite elasticity solver and create the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL SolverEquations_CreateStart(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquationsSolid,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquationsSolid,SOLVER_SPARSE_MATRICES,err,error,*999)
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Get the material-properties solver and create the material-properties solver equations
          NULLIFY(solverMatProperties)
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          NULLIFY(solverEquationsMatProperties)
          CALL SolverEquations_CreateStart(solverMatProperties,solverEquationsMatProperties,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquationsMatProperties,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquationsMatProperties,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquationsMatProperties,SOLVER_SPARSE_MATRICES,err,error,*999)
          !
          !Get the Darcy-ALE solver and create the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          NULLIFY(subiterationLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subiterationLoop,err,error,*999)
          NULLIFY(solidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,1,solidSubLoop,err,error,*999)
          NULLIFY(solidSolvers)
          CALL ControlLoop_SolversGet(solidSubLoop,solidSolvers,err,error,*999)
          !
          !Finish the creation of the finite elasticity solver equations
          NULLIFY(solverSolid)
          CALL Solvers_SolverGet(solidSolvers,1,solverSolid,err,error,*999)
          NULLIFY(solverEquationsSolid)
          CALL Solver_SolverEquationsGet(solverSolid,solverEquationsSolid,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquationsSolid,err,error,*999)             
          !
          NULLIFY(fluidSubLoop)
          CALL ControlLoop_SubLoopGet(subiterationLoop,2,fluidSubLoop,err,error,*999)
          NULLIFY(fluidSolvers)
          CALL ControlLoop_SolversGet(fluidSubLoop,fluidSolvers,err,error,*999)
          !
          !Finish the creation of the material-properties solver equations
          NULLIFY(solverMatProperties0
          CALL Solvers_SolverGet(fluidSolvers,1,solverMatProperties,err,error,*999)
          NULLIFY(solverEquationsMatProperties)
          CALL Solver_SolverEquationsGet(solverMatProperties,solverEquationsMatProperties,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquationsMatProperties,err,error,*999)             
          !
          !Finish the creation of the Darcy-ALE solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(fluidSolvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a finite elasticity ALE Darcy equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a finite elasticity ALE Darcy equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(problem%SPECIFICATION(3),"*",err,error))// &
        & " does not equal a standard finite elasticity Darcy equation subtype."
      CALL FlagError(localError,err,error,*999)      
    END SELECT
       
    EXITS("FiniteElasticityDarcy_ProblemSetup")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the finite elasticity Darcy problem pre-solve.
  SUBROUTINE FiniteElasticityDarcy_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopType,pSpecification(3),solverNumber,solverOutputType
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
      IF(loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.solverNumber==1) THEN
        CALL FiniteElasticity_PreSolve(solver,err,error,*999)
      ELSE IF(loopType==CONTROL_SIMPLE_TYPE) THEN
        IF(solverNumber==1) THEN
          !IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
          !  CALL WriteString(GENERAL_OUTPUT_TYPE,"Now working on material parameters.",err,error,*999)
          !ENDIF
        ELSE IF(solverNumber==2) THEN
          IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now working on Darcy",err,error,*999)
          ENDIF
        ENDIF
        CALL Darcy_PreSolve(solver,err,error,*999)
      ENDIF
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      IF(loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE.AND.solverNumber==1) THEN
        CALL FiniteElasticity_PreSolve(solver,err,error,*999)
      ELSE IF(loopType==CONTROL_SIMPLE_TYPE) THEN
        IF(solverNumber==1) THEN
          IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now working on material parameters",err,error,*999)
          ENDIF
        ELSE IF(solverNumber==2) THEN
          IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now working on Darcy",err,error,*999)
          ENDIF
        ENDIF
        CALL Darcy_PreSolve(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy fluid type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FiniteElasticityDarcy_PreSolve")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy  problem post solve.
  SUBROUTINE FiniteElasticityDarcy_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      !CALL FiniteElasticityDarcy_PostSolveOutputData(solver,err,error,*999)
      CALL FiniteElasticity_PostSolve(solver,err,error,*999)
      CALL Darcy_PostSolve(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a finite elasticity Darcy type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("FiniteElasticityDarcy_PostSolve")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_PostSolve

  !
  !================================================================================================================================
  !

  !>Runs before each control loop iteration
  SUBROUTINE FiniteElasticityDarcy_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iteratinNumber,loopType,outputType,pSpecification(3)
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoopDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverDarcy
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_PreLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    
    ! Eventually we may want to do different things depending on problem type/subtype
    ! too but for now we can just check the loop type.
    SELECT CASE(loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"==================================================",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"=============== Starting time step ===============",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Current Time          = ",currentTime,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Time Increment        = ",timeIncrement,err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"==================================================",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"==================================================",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"=============== Starting time step ===============",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Current Time          = ",currentTime,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Time Increment        = ",timeIncrement,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"==================================================",err,error,*999)
      ENDIF
      CALL Darcy_PreLoop(CONTROL_LOOP,err,error,*999)
      CALL FiniteElasticity_PreLoop(CONTROL_LOOP,err,error,*999)
      
    CASE(CONTROL_WHILE_LOOP_TYPE)
      !Subiteration loop
      CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"++++ Starting subiteration ++++",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Subiteration Number       =   ",iterationNumber,err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"++++ Starting subiteration ++++",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Subiteration Number       =   ",iterationNumber,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"+++++++++++++++++++++++++++++++",err,error,*999)
      ENDIF
      CALL ControlLoop_Get(controlLoop,[2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solverDarcy)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
        CALL Solvers_SolverGet(solvers,1,solverDarcy,err,error,*999)
      CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
        CALL Solvers_SolverGet(solvers,2,solverDarcy,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for FiniteElasticityDarcy_PreLoop."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL Darcy_PreSolveStorePreviousIterate(solverDarcy,err,error,*999)
    CASE(CONTROL_SIMPLE_TYPE)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"-- Starting fluid solve iteration --",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting fluid solve iteration --",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF      
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
      CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"-- Starting solid solve iteration --",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Load Increment number =           ",iterationNumber,err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"-- Starting solid solve iteration --",err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Load increment number =           ",iterationNumber,err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"------------------------------------",err,error,*999)
      ENDIF      
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("FiniteElasticityDarcy_PreLoop")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_PreLoop

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE FiniteElasticityDarcy_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopType,outputType,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoopDarcy
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverDarcy
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_PostLoop",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    
    SELECT CASE(loopType)
    CASE(CONTROL_TIME_LOOP_TYPE)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"End of time step",err,error,*999)
      ENDIF
    CASE(CONTROL_WHILE_LOOP_TYPE)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE)
        !subiteration
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"End of subiteration",err,error,*999)
        ENDIF
        NULLIFY(controlLoopDarcy)
        CALL ControlLoop_Get(controlLoop,[2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solverDarcy)
        CALL Solvers_SolverGet(solvers,1,solverDarcy,err,error,*999)
        !CALL Darcy_AccelerateConvergence(solverDarcy,err,error,*999)
        CALL Darcy_MonitorConvergence(solverDarcy,err,error,*999)
        !CALL Darcy_PostSolveOutputData(solverDarcy,err,error,*999)
      CASE(PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
        !subiteration
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"End of subiteration",err,error,*999)
        ENDIF
        NULLIFY(controlLoopDarcy)
        CALL ControlLoop_Get(controlLoop,[2,CONTROL_LOOP_NODE],controlLoopDarcy,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solverDarcy)
        CALL Solvers_SolverGet(solvers,2,solverDarcy,err,error,*999)
        !CALL Darcy_AccelerateConvergence(solverDarcy,err,error,*999)
        CALL Darcy_MonitorConvergence(solverDarcy,err,error,*999)
        !CALL Darcy_PostSolveOutputData(solverDarcy,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Darcy fluid type of a multi physics problem class with a while control loop."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(CONTROL_SIMPLE_TYPE)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"End of fluid solve iteration",err,error,*999)
      ENDIF
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"End of solid solve iteration",err,error,*999)
      ENDIF
    CASE DEFAULT
      !do nothing
    END SELECT

    EXITS("FiniteElasticityDarcy_PostLoop")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_PostLoop

  !
  !================================================================================================================================
  !

  !>Sets up the finite elasticity Darcy problem post solve output data.
  SUBROUTINE FiniteElasticityDarcy_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("FiniteElasticityDarcy_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_DARCY_SUBTYPE,PROBLEM_PGM_ELASTICITY_DARCY_SUBTYPE, &
      & PROBLEM_QUASISTATIC_ELASTICITY_TRANSIENT_DARCY_SUBTYPE,PROBLEM_QUASISTATIC_ELAST_TRANS_DARCY_MAT_SOLVE_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      IF(solverNumber==1) THEN
        CALL FiniteElasticity_PostSolveOutputData(solver,err,error,*999)
      ELSE IF(solverNumber==2.OR.solverNumber==3) THEN
        CALL Darcy_PostSolveOutputData(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy fluid type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("FiniteElasticityDarcy_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("FiniteElasticityDarcy_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE FiniteElasticityDarcy_PostSolveOutputData
      
  !   
  !================================================================================================================================
  !

END MODULE FiniteElasticityDarcyRoutines
