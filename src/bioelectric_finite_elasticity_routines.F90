!> \file
!> \authors Thomas Heidlauf
!> \brief This module handles all routines pertaining to bioelectrics coupled with finite elasticity.
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
!> Contributor(s): Thomas Heidlauf, Chris Bradley, Thomas Klotz
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

!> \defgroup OpenCMISS_BioelectricFiniteElasticity OpenCMISS::Iron::BioelectricFiniteElasticity
!>This module handles all routines pertaining to bioelectrics coupled with finite elasticity.
MODULE BioelectricFiniteElasticityRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BioelectricRoutines
  USE BiodomainEquationsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DecompositionAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FiniteElasticityRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PUBLIC BioelectricFiniteElasticity_EquationsSetSetup
  
  PUBLIC BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  PUBLIC BioelectricFiniteElasticity_ProblemSetup
  
  PUBLIC BioelectricFiniteElasticity_ProblemSpecificationSet
 
  PUBLIC BioelectricFiniteElasticity_FiniteElementCalculate

  PUBLIC BioelectricFiniteElasticity_PreSolve
  
  PUBLIC BioelectricFiniteElasticity_PostSolve

  PUBLIC BioelectricFiniteElasticity_PreLoop
  
  PUBLIC BioelectricFiniteElasticity_PostLoop
  
  PUBLIC BioelectricFiniteElasticity_UpdateGeometricField

CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectrics finite elasticity equation type of a multi physics equations set class.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
      & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
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
        & " is not valid for a bioelectrics finite elasticity equation type of a multi physics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet",err,error)
    EXITS("BioelectricFiniteElasticity_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity equation.
  SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("BioelectricFiniteElasticity_EquationsSetSetup",err,error,*999)

    CALL FlagError("BioelectricFiniteElasticity_EquationsSetSetup is not implemented.",err,error,*999)

    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_EquationsSetSetup",err,error)
    EXITS("BioelectricFiniteElasticity_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectrics finite elasticity equation finite element equations set.
  SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("BioelectricFiniteElasticity_FiniteElementCalculate",err,error,*999)

    CALL FlagError("BioelectricFiniteElasticity_FiniteElementCalculate is not implemented.",err,error,*999)

    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_FiniteElementCalculate",err,error)
    EXITS("BioelectricFiniteElasticity_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric finite elasticity problem type .
  SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE, &
      & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a bioelectric finite elasticity type of a multi physics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,problemSubtype]

    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ProblemSpecificationSet",err,error)
    EXITS("BioelectricFiniteElasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric finite elasticity problem.
  SUBROUTINE BioelectricFiniteElasticity_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(ControlLoopType), POINTER :: monodomainSubLoop,elasticitySubLoop
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers,monodomainSolvers,elasticitySolvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE, &
      & PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " does not equal a transient monodomain quasistatic finite elasticity equation subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    !--------------------------------------------------------------------
    !   Transient Gudunov monodomain, simple finite elasticity  
    !--------------------------------------------------------------------
    
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
          & " is invalid for a bioelectrics finite elasticity equation."
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
!!TODO: Default output type should be none???
        CALL ControlLoop_OutputTypeSet(controlLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        
        !Set up the control sub loop for monodomain
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        CALL ControlLoop_LabelSet(monodomainSubLoop,'Monodomain Time Loop',err,error,*999)
        CALL ControlLoop_TypeSet(monodomainSubLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
!!TODO: Default output type should be none???
        CALL ControlLoop_OutputTypeSet(monodomainSubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        
        !Set up the control sub loop for finite elasicity
        IF(pSpecification(3)==PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE) THEN
          NULLIFY(elasticitySubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
          CALL ControlLoop_LabelSet(elasticitySubLoop,'Elasticity While Loop',err,error,*999)
          CALL ControlLoop_TypeSet(elasticitySubLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
!!TODO: Default output type should be none???
          CALL ControlLoop_OutputTypeSet(elasticitySubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        ELSE
          NULLIFY(elasticitySubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
          CALL ControlLoop_LabelSet(elasticitySubLoop,'Elasticity Load increment loop',err,error,*999)
          CALL ControlLoop_TypeSet(elasticitySubLoop,CONTROL_LOAD_INCREMENT_LOOP_TYPE,err,error,*999)
!!TODO: Default output type should be none???
          CALL ControlLoop_OutputTypeSet(elasticitySubLoop,CONTROL_LOOP_PROGRESS_OUTPUT,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)
        !Sub-loops are finished when parent is finished
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectrics finite elasticity equation."
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
        !Get the monodomain sub loop
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        !Start the solvers creation
        NULLIFY(monodomainSolvers)
        CALL Solvers_CreateStart(monodomainSubLoop,monodomainSolvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(monodomainSolvers,2,err,error,*999)
        !Set the first solver to be a differential-algebraic equations solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
        CALL Solver_LabelSet(solver,"ODE Solver",err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        !Set the second solver to be a dynamic solver 
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,2,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        CALL Solver_LabelSet(solver,"Parabolic solver",err,error,*999)
        CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_DynamicRestartSet(solver,.TRUE.,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        
        !Get the finite elasticity sub loop
        NULLIFY(elasticitySubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
        !Start the solvers creation
        NULLIFY(elasticitySolvers)
        CALL Solvers_CreateStart(elasticitySubLoop,elasticitySolvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(elasticitySolvers,1,err,error,*999)
        !Set the finite elasticity solver to be a nonlinear solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(elasticitySolvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the monodomain solvers
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        NULLIFY(monodomainSolvers)
        CALL ControlLoop_SolversGet(monodomainSubLoop,monodomainSolvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(monodomainSolvers,err,error,*999)
        
        !Get the finite elasticity solvers
        NULLIFY(elasticitySubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
        NULLIFY(elasticitySolvers)
        CALL ControlLoop_SolversGet(elasticitySubLoop,elasticitySolvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(elasticitySolvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectrics finite elasticity equation."
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
        
        !Get the monodomain sub loop and solvers
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        NULLIFY(monodomainSolvers)
        CALL ControlLoop_SolversGet(monodomainSubLoop,monodomainSolvers,err,error,*999)
        !Create the solver equations for the second (parabolic) solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,2,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        
        !Get the finite elasticity sub loop and solvers
        NULLIFY(elasticitySubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
        NULLIFY(elasticitySolvers)
        CALL ControlLoop_SolversGet(elasticitySubLoop,elasticitySolvers,err,error,*999)
        !Get the finite elasticity solver and create the finite elasticity solver equations
        NULLIFY(solver)
        CALL Solvers_SolverGet(elasticitySolvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        
        !Get the monodomain sub loop and solvers
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        NULLIFY(monodomainSolvers)
        CALL ControlLoop_SolversGet(monodomainSubLoop,monodomainSolvers,err,error,*999)
        !Get the solver equations for the second (parabolic) solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,2,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
        
        !Get the finite elasticity sub loop and solvers
        NULLIFY(elasticitySubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
        NULLIFY(elasticitySolvers)
        CALL ControlLoop_SolversGet(elasticitySubLoop,elasticitySolvers,err,error,*999)
        !Finish the creation of the finite elasticity solver equations
        NULLIFY(solver)
        CALL Solvers_SolverGet(elasticitySolvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectrics finite elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        !Get the solvers
        NULLIFY(monodomainSolvers)
        CALL ControlLoop_SolversGet(monodomainSubLoop,monodomainSolvers,err,error,*999)
        !Create the CellML equations for the first DAE solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,1,solver,err,error,*999)
        NULLIFY(cellMLEquations)
        CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(monodomainSubLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,monodomainSubLoop,err,error,*999)
        !Get the solvers
        NULLIFY(monodomainSolvers)
        CALL ControlLoop_SolversGet(monodomainSubLoop,monodomainSolvers,err,error,*999)
        !Get the CellML equations for the first DAE solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(monodomainSolvers,1,solver,err,error,*999)
        NULLIFY(cellMLEquations)
        CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
        !Finish the CellML equations creation
        CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a bioelectrics finite elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a bioelectrics finite elasticity equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("BioelectricFiniteElasticity_ProblemSetup")
    RETURN
999 ERRORSEXITS("BioelectricFiniteElasticity_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectrics finite elasticity problem pre-solve.
  SUBROUTINE BioelectricFiniteElasticity_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopType,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
   
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL Biodomain_PreSolve(solver,err,error,*999)
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        CALL FiniteElasticity_PreSolve(solver,err,error,*999)
      CASE DEFAULT
        localError="Control loop type "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL Biodomain_PreSolve(solver,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL FiniteElasticity_PreSolve(solver,err,error,*999)
      CASE DEFAULT
        localError="Control loop type "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("BioelectricFiniteElasticity_PreSolve")
    RETURN
999 ERRORSEXITS("BioelectricFiniteElasticity_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_PreSolve
      
  !   
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post solve.
  SUBROUTINE BioelectricFiniteElasticity_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solveType
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
      & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      SELECT CASE(solveType)
      CASE(SOLVER_DAE_TYPE)
        CALL Bioelectric_PostSolve(solver,err,error,*999)
      CASE(SOLVER_DYNAMIC_TYPE)
        CALL Bioelectric_PostSolve(solver,err,error,*999)
      CASE(SOLVER_NONLINEAR_TYPE)
        CALL FiniteElasticity_PostSolve(solver,err,error,*999)
      CASE DEFAULT
        localError="Solver solve type "//TRIM(NumberToVString(solveType,"*",err,error))// &
          & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("BioelectricFiniteElasticity_PostSolve")
    RETURN
999 ERRORSEXITS("BioelectricFiniteElasticity_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem pre-control loop.
  SUBROUTINE BioelectricFiniteElasticity_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iterationNumber,loopType,numberOfSubLoops,pSpecification(3)
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_PreLoop",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      NULLIFY(problem)
      CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
      CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
      
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
        SELECT CASE(loopType)
        CASE(CONTROL_TIME_LOOP_TYPE)
          !do nothing ???
        CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
            CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(controlLoop,err,error,*999)
          CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
            CALL BioelectricFiniteElasticity_ComputeTitin(controlLoop,err,error,*999)
            CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(controlLoop,err,error,*999)
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for bioelectric finite elasticity problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL FiniteElasticity_PreLoop(controlLoop,err,error,*999)
        CASE(CONTROL_WHILE_LOOP_TYPE)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
            CALL ControlLoop_IterationNumberGet(controlLoop,iterationNumber,err,error,*999)
            IF(iterationNumber==1) CALL BioelectricFiniteElasticity_IndependentFieldInterpolate(controlLoop,err,error,*999)
            CALL BioelectricFiniteElasticity_ComputeFibreStretch(controlLoop,err,error,*999)
            CALL BioelectricFiniteElasticity_ForceLengthVelocityRelation(controlLoop,err,error,*999)
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for bioelectric finite elasticity problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL FiniteElasticity_PreLoop(controlLoop,err,error,*999)
        CASE DEFAULT
          localError="Control loop loop type "//TRIM(NumberToVString(loopType,"*",err,error))// &
            & " is not valid for bioelectric finite elasticity problem type."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second problem specification of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for a multi physics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      !the main time loop - do nothing!
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_PreLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_PreLoop",err,error)
    EXITS("BioelectricFiniteElasticity_PreLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_PreLoop

  !
  !================================================================================================================================
  !

  !>Computes the fibre stretch for bioelectrics finite elasticity problems.
  SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,elementIdx,equationsSetIdx,gaussPointIdx,iterationNumber,loopType,meshComponentNumber, &
      & numberOfElements,numberOfEquationsSets,numberOfGauss,numberOfIterations,numberOfSubLoops,numberOfXi,variableType
    REAL(DP) :: dZdNu(3,3),dZdNuT(3,3),AZL(3,3),Jznu
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,independentField
    TYPE(FieldInterpolationParametersType), POINTER :: geometricInterpParameters,fibreInterpParameters,dependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,dependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics,dependentInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: fieldVariableU1,residualVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_ComputeFibreStretch",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL ControlLoop_IterationNumberGet(controlLoop,numberOfIterations,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Loop over the equations sets associated with the solver
!!TODO: Why do we loop over the equations sets and get all the points if we do not actually do anything with those objects????
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(fibreField)
          CALL EquationsSet_FibreFieldGet(equationsSet,fibreField,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(independentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(nonlinearMapping)
          CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
        ENDDO !equationsSetIdx
        
        CALL Field_VariableGet(independentField,FIELD_U1_VARIABLE_TYPE,fieldVariableU1,err,error,*999)
        
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        NULLIFY(decompositionTopology)
        CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
        NULLIFY(decompositionElements)
        CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
        CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)

        IF(iterationNumber==1) THEN
          !copy the old fibre stretch to the previous values parameter set
          CALL FieldVariable_ParameterSetsCopy(fieldVariableU1,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP, &
            & err,error,*999)
        ENDIF
        
        !Grab interpolation parameters
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
        NULLIFY(residualVariable)
        CALL EquationsMappingResidual_VariableGet(residualMapping,1,residualVariable,err,error,*999)
        CALL FieldVariable_VariableTypeGet(residualVariable,variableType,err,error,*999)
        NULLIFY(equationsInterpolation)
        CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
        NULLIFY(geometricInterpParameters)
        CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpParameters,err,error,*999)
        NULLIFY(geometricInterpPoint)
        CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
          & err,error,*999)
        NULLIFY(geometricInterpPointMetrics)
        CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpPointMetrics,err,error,*999)
        NULLIFY(fibreInterpParameters)
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpParameters,err,error,*999)
        NULLIFY(fibreInterpPoint)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        NULLIFY(dependentInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,variableType, &
          & dependentInterpParameters,err,error,*999)
        NULLIFY(dependentInterpPoint)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,dependentInterpPoint,err,error,*999)
        NULLIFY(dependentInterpPointMetrics)
        CALL EquationsInterpolation_DependentPointMetricsGet(equationsInterpolation,variableType,dependentInterpPointMetrics, &
          & err,error,*999)
        
        !loop over the elements of the finite elasticity mesh (internal and boundary elements)
        DO elementIdx=1,numberOfElements

          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(dependentBasis)
          CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(dependentBasis,numberOfXi,err,error,*999)
          NULLIFY(dependentQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
          
          !Initialise tensors and matrices
          CALL IdentityMatrix(dZdNu,err,error,*999)
  
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,fibreInterpParameters,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,dependentInterpParameters,err,error,*999)
          
          !Loop over gauss points
          DO gaussPointIdx=1,numberOfGauss
            
            !Interpolate dependent, geometric, fibre and materials fields
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,dependentInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
              & err,error,*999)
            
            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(dependentInterpPointMetrics,geometricInterpPointMetrics, &
              & fibreInterpPoint,dZdNu,Jznu,err,error,*999)
              
            !compute C=F^T F
            CALL MatrixTranspose(dZdNu,dZdNuT,err,error,*999)
            CALL MatrixProduct(dZdNuT,dZdNu,AZL,err,error,*999)
            
            !store the fibre stretch lambda_f, i.e., sqrt(C_11) or AZL(1,1)            
            CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldVariableU1,FIELD_VALUES_SET_TYPE,gaussPointIdx,elementIdx, &
              & 1,SQRT(AZL(1,1)),err,error,*999)
            
          ENDDO !gaussPointIdx
        ENDDO !elementIdx

        !now the ghost elements -- get the relevant info from the other computational nodes
        CALL FieldVariable_ParameterSetUpdateStart(fieldVariableU1,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fieldVariableU1,FIELD_VALUES_SET_TYPE,err,error,*999)
        
      CASE DEFAULT
        localError="Control loop type "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE !the control loop contains subloops
      !do nothing
    ENDIF

    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ComputeFibreStretch",err,error)
    EXITS("BioelectricFiniteElasticity_ComputeFibreStretch")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ComputeFibreStretch

  !
  !================================================================================================================================
  !

  !>Computes the bioelectrics finite elasticity force-length and force_velocity relations.
  SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: counter,currentIteration,dofIdx,elementIdx,equationsSetIdx,gaussPointIdx,iterationNumber,loopType, &
      & maximumIterations,meshComponentNumber,numberOfElements,numberOfEquationsSets,numberOfGauss,numberOfSubLoops, &
      & numberOfXi
    REAL(DP) :: absoluteTolerance,averageVelocity,averageStretch,currentTime,oldAverageStretch,lengthHS,lengthHS0,activeStress, &
      & fibreStretch,oldFibreStretch,lengthFactor,velocityFactor,relativeTolerance,sarcomereLength,timeIncrement,velocity, &
      & maxVelocity,timeStep,kappa,A,S,d,c
    LOGICAL :: continueLoop
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(ControlLoopType), POINTER :: controlLoopParent
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,independentField
    TYPE(FieldInterpolationParametersType), POINTER :: geometricInterpParameters,fibreInterpParameters,dependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,dependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics,dependentInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: uIndependentVariable,u1IndependentVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(controlLoop%numberOfSubLoops==0) THEN
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL ControlLoop_CurrentWhileInformationGet(controlLoop,currentIteration,maximumIterations,absoluteTolerance, &
          & relativeTolerance,continueLoop,err,error,*999)
        CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
        NULLIFY(solvers)        
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Loop over the equations sets associated with the solver
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        !!TODO: WHY ARE WE LOOPING OVER EQUATIONS SETS THUS T0 GET POINTERS TO FIELDS IF THEY ARE NOT USED???
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(independentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        ENDDO !equationsSetIdx
 
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,uIndependentVariable,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_U1_VARIABLE_TYPE,u1IndependentVariable,err,error,*999)

        !get the inital half-sarcomere length
        CALL FieldVariable_ParameterSetGetConstant(u1IndependentVariable,FIELD_VALUES_SET_TYPE,2,lengthHS0,err,error,*999)

        !get the maximum contraction velocity
        CALL FieldVariable_ParameterSetGetConstant(u1IndependentVariable,FIELD_VALUES_SET_TYPE,3,maxVelocity,err,error,*999)
        
        !in the first iteration store the unaltered homogenized active stress field 
        IF(iterationNumber==1) THEN
          CALL FieldVariable_ParameterSetsCopy(uIndependentVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP, &
            & err,error,*999)
        ELSE
          !restore the solution of the previous time step
          CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
            & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        ENDIF

        timeStep=timeIncrement

        NULLIFY(decomposition)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        NULLIFY(decompositionTopology)
        CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
        NULLIFY(decompositionElements)
        CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
        CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
        
        averageVelocity=0.0_DP
        averageStretch=0.0_DP
        oldAverageStretch=0.0_DP
        counter=0

        !loop over the elements of the finite elasticity mesh (internal and boundary elements)
        DO elementIdx=1,numberOfElements

          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          NULLIFY(dependentBasis)
          CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(dependentBasis,numberOfXi,err,error,*999)
          NULLIFY(dependentQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)

          !Loop over gauss points
          DO gaussPointIdx=1,numberOfGauss

            !get the unaltered activeStress at the GP
            CALL FieldVariable_ParameterSetGetGaussPoint(uIndependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,gaussPointIdx, &
              & elementIdx,1,activeStress,err,error,*999)
 
            ! FORCE-LENGTH RELATION -------------------------------------------------------------------------------
              
            !get the current fibre stretch at the GP
            CALL FieldVariable_ParameterSetGetGaussPoint(u1IndependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,gaussPointIdx, &
              & elementIdx,1,fibreStretch,err,error,*999)
            
            !compute the current half-sarcomere length at the GP: l_hs = lambda_f * l_hs_0
            lengthHS=lengthHS0*fibreStretch
              
            !compute the scale factor (0,1) due to sarcomere F-l relation of Gordon, A. M., A.F. Huxley, and F.J. Julian. 
            !The variation in isometric tension with sarcomere length in vertebrate muscle fibres. 
            !The Journal of Physiology 184.1 (1966): 170-192.
            sarcomereLength=2.0_DP*lengthHS
            IF(sarcomereLength<=1.27_DP) THEN
              lengthFactor=0.0_DP
            ELSE IF(sarcomereLength<=1.7_DP) THEN
              lengthFactor=1.6047_DP*sarcomereLength-2.0379_DP
            ELSE IF(sarcomereLength<=2.34_DP) THEN
              lengthFactor=0.4844_DP*sarcomereLength-0.1334_DP
            ELSE IF(sarcomereLength<=2.51_DP) THEN
              lengthFactor=1.0_DP
            ELSEIF(sarcomereLength<=3.94_DP) THEN
              lengthFactor=-0.6993_DP*sarcomereLength+2.7552_DP
            ELSE
              lengthFactor=0.0_DP
            ENDIF
            
            !multiply the activeStress with the scale factor
            activeStress=activeStress*lengthFactor
            
            !update the activeStress at GP
            CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(uIndependentVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
              & elementIdx,1,activeStress,err,error,*999)
            
            ! FORCE-VELOCITY RELATION -------------------------------------------------------------------------------
              
            !get fibre stretch at the GP of the previous time step
            CALL FieldVariable_ParameterSetGetLocalGaussPoint(u1IndependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,gaussPointIdx, &
              & elementIdx,1,oldFibreStretch,err,error,*999)

            !compute the contraction velocity
            velocity=(fibreStretch-oldFibreStretch)/timeStep
            
            !NOTE: maxVelocity is the max shortening velocity, and hence negative
            IF(velocity<maxVelocity) THEN
              CALL FlagWarning('Exceeded maximum contraction velocity (shortening).',err,error,*999)
              IF(iterationNumber<(maximumIterations/2)) THEN
                velocity=velocity*1.0_DP/DBLE((maximumIterations/2)-iterationNumber)
              ENDIF
            ELSE IF(velocity>(ABS(maxVelocity))) THEN
              CALL FlagWarning('Exceeded maximum contraction velocity (lengthening).',err,error,*999)
            ENDIF

            !compute scale factor (velocityFactor)
            kappa=0.24_DP !make mat param TODO
            A=1.6_DP !make mat param TODO
            S=2.5_DP !make mat param TODO
            IF(velocity<0.0_DP) THEN
              !shortening contraction              
              velocityFactor=(1.0_DP-velocity/maxVelocity)/(1.0_DP+velocity/maxVelocity/kappa)
            ELSE
              !lengthening contraction
              d=kappa*(1.0_DP-A)
              c=velocity/maxVelocity*S*(kappa+1)
              velocityFactor=1.0_DP+c*(A-1.0_DP)/(d+c)
            ENDIF
              
            !multiply the activeStress with the scale factor
            activeStress=activeStress*velocityFactor

            !update the activeStress at GP
            CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(uIndependentVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
              & elementIdx,1,activeStress,err,error,*999)

          ENDDO !gaussPointIdx
        ENDDO !elementIdx

!!!!tomo --- try averaging              
!!!              averageVelocity=averageVelocity+velocity
!!!              averageStretch=averageStretch+fibreStretch
!!!              oldAverageStretch=oldAverageStretch+oldFibreStretch
!!!              
!!!              counter=counter+1

!!!            ENDDO !gaussPointIdx
!!!          ENDDO !elementIdx

!!!          velocity=averageVelocity/REAL(counter)
!!!          fibreStretch=averageStretch/REAL(counter)
!!!          oldFibreStretch=oldAverageStretch/REAL(counter)

!!!          averageVelocity=(fibreStretch-oldFibreStretch)/timeStep

!!!!          velocity=0.0_DP

!!!!          !damping
!!!!          IF(iterationNumber<(maximumIterations/2)) THEN
!!!!            velocity=velocity*1.0_DP/DBLE((maximumIterations/2)-iterationNumber)
!!!!          ENDIF

!!!          localError="######### velocity: "//TRIM(NumberToVString(velocity,"*",err,error))//" #########"
!!!          CALL FlagWarning(localError,err,error,*999)
!!!          localError="######### velocity AVERAGE: "//TRIM(NumberToVString(averageVelocity,"*",err,error))//" #########"
!!!          CALL FlagWarning(localError,err,error,*999)
!!!          localError="######### STRETCH: "//TRIM(NumberToVString(fibreStretch,"*",err,error))//" #########"
!!!          CALL FlagWarning(localError,err,error,*999)
!!!          localError="######### OLD STRETCH: "//TRIM(NumberToVString(oldFibreStretch,"*",err,error))//" #########"
!!!          CALL FlagWarning(localError,err,error,*999)


!!!          
!!!          
!!!          velocity=averageVelocity
!!!          
!!!          
!!!          

!!!          !compute scale factor (velocityFactor)
!!!          kappa=0.24_DP !make mat param TODO
!!!          A=1.6_DP !make mat param TODO
!!!          S=2.5_DP !make mat param TODO
!!!          IF(velocity<0.0_DP) THEN
!!!            !shortening contraction              
!!!            velocityFactor=(1-velocity/maxVelocity)/(1+velocity/maxVelocity/kappa)
!!!            IF(velocity<maxVelocity) THEN
!!!              CALL FlagWarning('Exceeded maximum contraction velocity (shortening).',err,error,*999)
!!!            ENDIF
!!!          ELSE
!!!            !lengthening contraction
!!!            d=kappa*(1-A)
!!!            c=velocity/maxVelocity*S*(kappa+1)
!!!            velocityFactor=1+c*(A-1)/(d+c)
!!!            IF(velocity>ABS(maxVelocity)) THEN
!!!              CALL FlagWarning('Exceeded maximum contraction velocity (lengthening).',err,error,*999)
!!!            ENDIF
!!!          ENDIF

!!!          DO elementIdx=1,numberOfElements

!!!            dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)%BASIS       
!!!            dependentQuadratureScheme=>dependentBasis%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
!!!            numberOfGauss=dependentQuadratureScheme%numberOfGauss

!!!            !Loop over gauss points
!!!            DO gaussPointIdx=1,numberOfGauss

!!!              !get the activeStress at the GP
!!!              dofIdx=uIndependentVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gaussPointIdx,elementIdx)
!!!              CALL Field_ParameterSetGetLocalDOF(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & dofIdx,activeStress,err,error,*999)

!!!              !multiply the activeStress with the scale factor
!!!              activeStress=activeStress*velocityFactor

!!!              !update the activeStress at GP
!!!              dofIdx=uIndependentVariable%COMPONENTS(1)%PARAM_TO_DOF_MAP%GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(gaussPointIdx,elementIdx)
!!!              CALL Field_ParameterSetUpdateLocalDOF(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
!!!                & dofIdx,activeStress,err,error,*999)
!!!              
!!!            ENDDO !gaussPointIdx
!!!          ENDDO !elementIdx
!tomo end

        !now the ghost elements -- get the relevant info from the other computational nodes
        CALL FieldVariable_ParameterSetUpdateStart(uIndependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(uIndependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        
      CASE DEFAULT
        localError="Control loop type "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is not valid for a bioelectrics finite elasticity type of a multi physics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE !the control loop contains subloops
      !do nothing
    ENDIF

    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ForceLengthVelocityRelation",err,error)
    EXITS("BioelectricFiniteElasticity_ForceLengthVelocityRelation")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ForceLengthVelocityRelation

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectrics finite elasticity problem post-control loop.
  SUBROUTINE BioelectricFiniteElasticity_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,equationsSetIdx,inputIteration,loopType,numberOfEquationsSets,numberOfSubLoops, &
      & outputIteration,pSpecification(3),regionUserNumber
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: elasticitySubLoop,bioelectricSubLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldsType), POINTER :: fields
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region   
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: filename,localError,method

    ENTERS("BioelectricFiniteElasticity_PostLoop",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      NULLIFY(problem)
      CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
      CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        SELECT CASE(pSpecification(2))
        CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
          !the monodomain time loop - output of the monodomain fields
          CALL Biodomain_PostLoop(controlLoop,err,error,*999)
        CASE DEFAULT
          localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
            & " is not valid for a multi physics problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        CALL BioelectricFiniteElasticity_UpdateGeometricField(controlLoop,.FALSE.,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL BioelectricFiniteElasticity_ConvergenceCheck(controlLoop,err,error,*999)
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE
      !the main time loop - output the finite elasticity fields
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      IF(outputIteration/=0) THEN
        IF(MOD(currentIteration,outputIteration)==0) THEN
          NULLIFY(problem)
          CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
          CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
          !Get the solver. The first solver of the second sub loop will contain the finite elasticity dependent field equation set
          NULLIFY(elasticitySubLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,elasticitySubLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(elasticitySubLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          !Loop over the equations sets associated with the solver
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            NULLIFY(region)
            CALL Field_RegionGet(dependentField,region,err,error,*999)
            CALL Region_UserNumberGet(region,regionUserNumber,err,error,*999)
            NULLIFY(fields)
            CALL Region_FieldsGet(region,fields,err,error,*999)
            filename="MainTime_"//TRIM(NumberToVString(regionUserNumber,"*",err,error))// &
              & "_"//TRIM(NumberToVString(currentIteration,"*",err,error))
            method="FORTRAN"
            CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
          ENDDO !equationsSetIdx
          IF((pSpecification(3)==PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE).OR. &
            & (pSpecification(3)==PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE).OR. &
            & (pSpecification(3)==PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE).OR. &
            & (pSpecification(3)==PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)) THEN
            !Get the solver. The second solver of the first sub loop will contain the bioelectrics equation set
            NULLIFY(bioelectricSubLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,bioelectricSubLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(bioelectricSubLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            NULLIFY(solverMapping)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            !Loop over the equations sets associated with the solver
            CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
            DO equationsSetIdx=1,numberOfEquationsSets
              NULLIFY(equationsSet)
              CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              NULLIFY(region)
              CALL Field_RegionGet(dependentField,region,err,error,*999)
              CALL Region_UserNumberGet(region,regionUserNumber,err,error,*999)
              NULLIFY(fields)
              CALL Region_FieldsGet(region,fields,err,error,*999)
              filename="MainTime_M_"//TRIM(NumberToVString(regionUserNumber,"*",err,error))// &
                & "_"//TRIM(NumberToVString(currentIteration,"*",err,error))
              method="FORTRAN"
              CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
            ENDDO !equationsSetIdx
          ENDIF !PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE
        ENDIF
      ENDIF
    ENDIF

    EXITS("BioelectricFiniteElasticity_PostLoop")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_PostLoop",err,error)
    EXITS("BioelectricFiniteElasticity_PostLoop")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_PostLoop


  !
  !================================================================================================================================
  !

  !>Check for the convergence of the bioelectric finite elasticity while loop, i.e., if the force-length and force-velocity relations yielded a different actual configuration.
  SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,dofIdx,equationsSetIdx,loopType,maximumIterations,nodeIdx,numberOfEquationsSets, &
      & numberOfNodes,numberOfSubLoops
    REAL(DP) :: absoluteTolerance,relativeTolerance,x1,x2,x3,y1,y2,y3,mySum
    LOGICAL :: continueLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("BioelectricFiniteElasticity_ConvergenceCheck",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL ControlLoop_CurrentWhileInformationGet(controlLoop,currentIteration,maximumIterations,absoluteTolerance, &
          & relativeTolerance,continueLoop,err,error,*999)
        CALL ControlLoop_IterationNumberGet(controlLoop,currentIteration,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_SolutionUpdate(solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Loop over the equations sets associated with the solver
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
!!TODO: Why are we looping over equations sets to get pointers if we don't use them???
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        ENDDO !equationsSetIdx
          
        IF(currentIteration>1) THEN
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          NULLIFY(decomposition)
          CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
  
          mySum=0.0_DP
          DO nodeIdx=1,numberOfNodes
            
            !get the current node position (x) and the node position of the last iteration (y)
            CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,1,nodeIdx,1,dofIdx,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,dofIdx,x1,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dofIdx,y1, &
              & err,error,*999)
            CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,1,nodeIdx,2,dofIdx,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,dofIdx,x2,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dofIdx,y2, &
              & err,error,*999)
            CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,1,nodeIdx,3,dofIdx,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,dofIdx,x3,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalDOF(fieldVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,dofIdx,y3, &
              & err,error,*999)            
            !sum up
            mySum=mySum+SQRT((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
          ENDDO !nodeIdx
          
          IF(mySum<1.0E-06_DP) THEN !if converged then:
            CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
            CALL BioelectricFiniteElasticity_UpdateGeometricField(controlLoop,.FALSE.,err,error,*999)
          ELSE IF(currentIteration==maximumIterations) THEN
            CALL BioelectricFiniteElasticity_UpdateGeometricField(controlLoop,.FALSE.,err,error,*999)
            CALL FlagWarning('----------- Maximum number of iterations in while loop reached. -----------',err,error,*999)
          ENDIF          
        ENDIF
        
        !copy the current solution to the previous solution
        CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP, &
          & err,error,*999)
        
      CASE DEFAULT
        !do nothing
      END SELECT
    ELSE !the main time loop
      !do nothing
    ENDIF

    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_ConvergenceCheck",err,error)
    EXITS("BioelectricFiniteElasticity_ConvergenceCheck")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_ConvergenceCheck

  !
  !================================================================================================================================
  !

  !>Update the the bioelectric equation geometric field from the finite elasticity dependent field (deformed geometry)
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField(controlLoop,calcClosestGaussPoint,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    LOGICAL, INTENT(IN) :: calcClosestGaussPoint !<If true then the closest finite elasticity Gauss point for each bioelectrics node is calculated
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,dependentFieldInterpolation,dependentUserNumber,dependentVariableType,dofIdx,dofIdx2, &
      & elementIdx,elementIdx2,elementNumber,elementNumber2,fibreIdx,gaussPoint,gaussPointIdx,geometricFieldInterpolation, &
      & geometricUserNumber,geometricVariableType,myElementIdx,n1,n2,n3,n4,nodeIdx,nodeIdx2,nodesInXi1,nodesInXi2,nodesInXi3, &
      & numberOfAdjacentElements,numberOfComponents,numberOfElements,numberOfGauss,numberOfSubLoops,offset,pSpecification(3), &
      & startElem,startElement,startElementIdx
    REAL(DP) :: currentTime,distance,distanceLeft,distanceRight,elasticityXValue,gaussPosition(3),h,initialDistance, &
      & initialSarcomereLength,maxVelocity,monodomainXValue,oldDistance,oldDistance2,oldDistance3,oldDistance4,previousNode(3), &
      & timeIncrement,timeStep,VALUE,valueLeft,valueRight,velocity,xi(3)
    LOGICAL :: outsideNode
    TYPE(BasisType), POINTER :: basis
    TYPE(ControlLoopType), POINTER :: controlLoopRoot,controlLoopParent,elasticityControlLoop,monodomainControlLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: elasticityDependentField,elasticityGeometricField,elasticityIndependentField, &
      & monodomainDependentField,monodomainGeometricField,monodomainIndependentField
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldVariableType), POINTER :: elasticityDependentVariable,elasticityIndependentVariable,monodomainDependentVariable, &
      & monodomainGeometricVariable,monodomainIndependentVariable,monodomainIndependentVariable2
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_UpdateGeometricField",err,error,*999)

    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      NULLIFY(problem)
      CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
      CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      NULLIFY(controlLoopParent)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoopParent,err,error,*999)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE)

          !get the monodomain sub loop, solvers, solver, and finally geometric and field
          NULLIFY(monodomainControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,1,monodomainControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(monodomainControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(monodomainGeometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,monodomainGeometricField,err,error,*999)
          !get the finite elasticity sub loop, solvers, solver, and finally the dependent field
          NULLIFY(elasticityControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,2,elasticityControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(elasticityControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(elasticityDependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,elasticityDependentField,err,error,*999)
          NULLIFY(monodomainGeometricVariable)
          CALL Field_VariableIndexGet(monodomainGeometricField,1,monodomainGeometricVariable,geometricVariableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(monodomainGeometricVariable,numberOfComponents,err,error,*999)
          NULLIFY(elasticityDependentVariable)
          CALL Field_VariableIndexGet(elasticityDependentField,1,elasticityDependentVariable,dependentVariableType,err,error,*999)
          DO componentIdx=1,numberOfComponents
            !check for identical interpolation of the fields
            CALL FieldVariable_ComponentInterpolationGet(monodomainGeometricVariable,componentIdx,geometricFieldInterpolation, &
              & err,error,*999)
            CALL FieldVariable_ComponentInterpolationGet(elasticityDependentVariable,componentIdx,dependentFieldInterpolation, &
              & err,error,*999)
            IF(geometricFieldInterpolation==dependentFieldInterpolation) THEN
              !copy the dependent field components to the geometric field
              CALL Field_ParametersToFieldParametersCopy(elasticityDependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,monodomainGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,err,error,*999)
            ELSE
              CALL Field_UserNumberGet(monodomainGeometricField,geometricUserNumber,err,error,*999)
              CALL Field_UserNumberGet(elasticityDependentField,dependentUserNumber,err,error,*999)
              localError="The interpolation type of component number "// &
                & TRIM(NumberToVString(componentIdx,"*",err,error))//" of field number "// &
                & TRIM(NumberToVString(geometricUserNumber,"*",err,error))// &
                & " does not coincide with the interpolation type of field number "// &
                & TRIM(NumberToVString(dependentUserNumber,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDDO !componentIdx

        CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)

          CALL ControlLoop_CurrentTimesGet(controlLoopParent,currentTime,timeIncrement,err,error,*999)
          !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
          NULLIFY(monodomainControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,1,monodomainControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(monodomainControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(monodomainGeometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,monodomainGeometricField,err,error,*999)
          ! the Field_V_Variable_Type contains the 3D nodal positions
          NULLIFY(monodomainDependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,monodomainDependentField,err,error,*999)
          NULLIFY(monodomainIndependentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,monodomainIndependentField,err,error,*999)
          !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
          NULLIFY(elasticityControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,2,elasticityControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(elasticityControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(elasticityGeometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,elasticityGeometricField,err,error,*999)
          NULLIFY(elasticityDependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,elasticityDependentField,err,error,*999)
          NULLIFY(elasticityIndependentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,elasticityIndependentField,err,error,*999)
          !get the finite elasticity dependent field interpolation parameters of this element
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(equationsInterpolation)
          CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
          NULLIFY(interpolationParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & interpolationParameters,err,error,*999)
          NULLIFY(interpolatedPoint)
          CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,interpolatedPoint, &
            & err,error,*999)

          nodeIdx=0
          nodeIdx2=0
          fibreIdx=0
          NULLIFY(monodomainDependentVariable)
          CALL Field_VariableGet(monodomainDependentField,FIELD_V_VARIABLE_TYPE,monodomainDependentVariable,err,error,*999)
          NULLIFY(monodomainGeometricVariable)
          CALL Field_VariableGet(monodomainGeometricField,FIELD_U_VARIABLE_TYPE,monodomainGeometricVariable,err,error,*999)
          NULLIFY(elasticityIndependentVariable)
          CALL Field_VariableGet(elasticityIndependentField,FIELD_V_VARIABLE_TYPE,elasticityIndependentVariable,err,error,*999)
          NULLIFY(monodomainIndependentVariable)
          CALL Field_VariableGet(monodomainIndependentField,FIELD_U1_VARIABLE_TYPE,monodomainIndependentVariable,err,error,*999)
          NULLIFY(monodomainIndependentVariable2)
          CALL Field_VariableGet(monodomainIndependentField,FIELD_U2_VARIABLE_TYPE,monodomainIndependentVariable2,err,error,*999)

          NULLIFY(decomposition)
          CALL Field_DecompositionGet(elasticityGeometricField,decomposition,err,error,*999)
          NULLIFY(decompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
          NULLIFY(decompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)

          !get the maximum contraction velocity 
          CALL FieldVariable_ParameterSetGetConstant(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE,2,maxVelocity, &
            & err,error,*999)
          !NOTE: maxVelocity is the max shortening velocity, and hence negative!!!
          !The max lengthening velocity is assumed to be abs(maxVelocity)/2.0

          !get the time step of the elasticity problem
          timeStep=timeIncrement

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          !no need to consider ghost elements here since only bioelectrical fields are changed
          CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
          DO elementIdx=1,numberOfElements

!!TODO is elementNumber not just the same as elementIdx????
            elementNumber=decompositionElements%elements(elementIdx)%localNumber
            myElementIdx=elementIdx

            !the Field_V_Variable_Type of the FE independent field contains the number of nodes in each Xi-direction of
            !the bioelectrics grid
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
              & 1,nodesInXi1,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
              & 2,nodesInXi2,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
              & 3,nodesInXi3,err,error,*999)
            !beginning of a fibre in this element: 1=yes, 0=no
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
              & 4,startElem,err,error,*999)            

            !if there is no bioelectrics grid in this finite elasticity element, or the fibres don't begin in this element,
            IF((nodesInXi1==0).OR.(nodesInXi2==0).OR.(nodesInXi3==0).OR.(startElem==0)) CYCLE

            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,interpolationParameters,err,error,*999)
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
            NULLIFY(quadratureScheme)
            CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)

            startElement=elementNumber
            startElementIdx=myElementIdx

            !assume Xi(1) to be normal to the seed surface, i.e. the seed points have xi(1)=0
            xi=[0.0_DP,1.0_DP/(REAL(2*nodesInXi2)),1.0_DP/(REAL(2*nodesInXi3))]

            !assume that the bioelectrics node numbers are increased in order Xi(1), Xi(2), Xi(3) 
            DO n3=1,nodesInXi3
              DO n2=1,nodesInXi2
                fibreIdx=fibreIdx+1

                !loop over the FE elements that contain nodes of the very same fibres
                DO                      
                  DO n1=1,nodesInXi1
                    nodeIdx=nodeIdx+1

                    !store the fibre number this bioelectrics node belongs to.
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE,1,1, &
                      & nodeIdx,3,fibreIdx,err,error,*999) 

                    !find the interpolated position of the bioelectric grid node from the FE dependent field
                    CALL Field_InterpolateXi(NO_PART_DERIV,XI,interpolatedPoint,err,error,*999)
                    !update the bioelectrics dependent field Field_V_Variable_Type
                    !the Field_V_Variable_Type of the monodomain dependent field contains the nodal positions in 3D
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,1,interpolatedPoint%values(1,NO_PART_DERIV),err,error,*999)
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,2,interpolatedPoint%values(2,NO_PART_DERIV),err,error,*999)
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,3,interpolatedPoint%values(3,NO_PART_DERIV),err,error,*999)

                    IF((n1==1).AND.(elementNumber==startElement)) THEN
                      !a new line of bioelectrics grid nodes begins
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,0.0_DP,err,error,*999)
                    ELSE
                      !get the position in 3D of the previous node
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,1,previousNode(1),err,error,*999)
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,2,previousNode(2),err,error,*999)
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,3,previousNode(3),err,error,*999)

                      !compute the distance between the previous node and the actual node
                      distance=SQRT((interpolatedPoint%values(1,NO_PART_DERIV)- &
                        & previousNode(1))*(interpolatedPoint%values(1,NO_PART_DERIV)-previousNode(1))+ &
                        & (interpolatedPoint%values(2,NO_PART_DERIV)- &
                        & previousNode(2))*(interpolatedPoint%values(2,NO_PART_DERIV)-previousNode(2))+ &
                        & (interpolatedPoint%values(3,NO_PART_DERIV)- &
                        & previousNode(3))*(interpolatedPoint%values(3,NO_PART_DERIV)-previousNode(3)))

                      !CONTRACTION VELOCITY CALCULATION                          
                      !get the distance between the 2 nodes in the previous time step
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,oldDistance,err,error,*999)                      
                      !get the distance between the 2 nodes before two time step
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,4,oldDistance2,err,error,*999)                      
                      !get the distance between the 2 nodes before three time step
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,5,oldDistance3,err,error,*999)                      
                      !get the distance between the 2 nodes before four time step
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,6,oldDistance4,err,error,*999)

                      !compute the new contraction velocity
                      !velocity=(distance-oldDistance)/timeStep
                      velocity=0.25_DP*((distance-oldDistance)/timeStep+(distance-oldDistance2)/(2.0_DP*timeStep)+ &
                        & (distance-oldDistance3)/(3.0_DP*timeStep)+(distance-oldDistance4)/(4.0_DP*timeStep))
                      IF(.NOT.calcClosestGaussPoint) THEN
                        !NOTE: maxVelocity is the max shortening velocity, and hence negative!!!
                        IF(velocity<maxVelocity) THEN
                          CALL FlagWarning('Exceeded maximum contraction velocity (shortening).',err,error,*999)
                          velocity=maxVelocity
                          !The max lengthening velocity is assumed to be maxVelocity/2.0
                        ELSEIF(velocity>(ABS(maxVelocity)/2.0_DP)) THEN
                          CALL FlagWarning('Exceeded maximum contraction velocity (lengthening).',err,error,*999)
                          velocity=-maxVelocity/2.0_DP
                        ENDIF
                      ENDIF

                      !store the relative contraction velocity in component 3 of the U2 variable of the monodomain independent field
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,3,velocity/ABS(maxVelocity),err,error,*999)                      
                      !store the node distance for contraction velocity calculation
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,distance,err,error,*999)                      
                      !store the node distance for contraction velocity calculation
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,4,oldDistance,err,error,*999)                      
                      !store the node distance for contraction velocity calculation
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,5,oldDistance2,err,error,*999)                      
                      !store the node distance for contraction velocity calculation
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,6,oldDistance3,err,error,*999)

                      !get the position in 1D of the previous node
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,1,valueLeft,err,error,*999)
                      !update the current 1D node position
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,valueLeft+distance,err,error,*999)

                      !get the initial sarcomere half length and initial node distance
                      CALL FieldVariable_ParameterSetGetConstant(monodomainIndependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 2,initialSarcomereLength,err,error,*999)
                      CALL FieldVariable_ParameterSetGetConstant(monodomainIndependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 3,initialDistance,err,error,*999)
                      !update the current sarcomere half length
                      VALUE=initialSarcomereLength*distance/initialDistance
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,VALUE,err,error,*999)

                      !update the first node to the same value as the second node (no better info available)
                      IF((n1==2).AND.(elementNumber==startElement)) THEN
                        !current sarcomere half length
                        CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable,FIELD_VALUES_SET_TYPE, &
                          & 1,1,nodeIdx-1,1,VALUE,err,error,*999)
                        !old node distance
                        CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                          & 1,1,nodeIdx-1,1,distance,err,error,*999)
                        !relative contraction velocity
                        CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                          & 1,1,nodeIdx-1,3,velocity/ABS(maxVelocity),err,error,*999)
                      ENDIF

                    ENDIF !((n1==1).AND.(elementNumber==startElement))

                    IF(calcClosestGaussPoint) THEN
                      !calculate the closest finite elasticity Gauss point of each bioelectrics node
                      distance=1000000.0_DP
                      gaussPoint=0
                      DO gaussPointIdx=1,numberOfGauss
                        !compute the distance between the bioelectrics node and the Gauss point
                        CALL BasisQuadratureScheme_GaussPositionGet(quadratureScheme,gaussPointIdx,gaussPosition,err,error,*999)
                        VALUE=SQRT((Xi(1)-gaussPosition(1))*(xi(1)-gaussPosition(1))+ &
                          & (xi(2)-gaussPosition(2))*(xi(2)-gaussPosition(2))+ &
                          & (xi(3)-gaussPosition(3))*(xi(3)-gaussPosition(3)))
                        IF(value<distance) THEN
                          distance=VALUE
                          gaussPoint=gaussPointIdx
                        ENDIF
                      ENDDO !gaussPointIdx
                      IF(gaussPoint==0) CALL FlagWarning("Closest Gauss Point not found",err,error,*999)
                      !store the nearest Gauss Point info and the inElement info (local element number!!!)
                      CALL Field_ParameterSetUpdateLocalNode(monodomainIndependentField,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4,gaussPoint,err,error,*999)
                      CALL Field_ParameterSetUpdateLocalNode(monodomainIndependentField,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5,elementNumber,err,error,*999)
                    ENDIF !calcClosestGaussPoint

                    IF(startElem==1) THEN
                      !fibres start in this element
                      xi(1)=xi(1)+1.0_DP/(REAL(nodesInXi1-1))
                    ELSEIF(startElem==0) THEN
                      !fibres don't start in this element
                      xi(1)=xi(1)+1.0_DP/(REAL(nodesInXi1))
                    ELSE
                      localError="The start element index is incorrect. The index is "// &
                        & TRIM(NumberToVString(startElem,"*",err,error))//" and should be zero or one." 
                      CALL FlagError(localError,err,error,*999)
                    ENDIF

                  ENDDO !n1

!!smooth of the velocity field
!!arithmetic mean of all rel_velo values within one FE element
                  !velocity=0.0_DP
                  !DO n1=1,nodesInXi1
                  !  nodeIdx2=nodeIdx2+1
                  !  CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !  CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !    & dofIdx,VALUE,err,error,*999)
                  !  velocity=velocity+VALUE
                  !ENDDO !n1
                  !velocity=velocity/nodesInXi1
                  !
                  !nodeIdx2=nodeIdx2-nodesInXi1
                  !DO n1=1,nodesInXi1
                  !  nodeIdx2=nodeIdx2+1
                  !  CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !    & 1,1,nodeIdx2,3,velocity,err,error,*999)
                  !ENDDO !n1
                  !
!!--------------------------------------------------------------------------
                  !
!!moving average
                  !offset=3                  
!!do the first three nodes of a fibre manually - arithmetic mean
                  !velocity=0.0_DP
                  !nodeIdx2=nodeIdx2+1                  
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE
                  !nodeIdx2=nodeIdx2+1
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE
                  ! 
                  !nodeIdx2=nodeIdx2+1
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE                  
                  !velocity=velocity/offset
                  !
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & 1,1,nodeIdx2-2,3,velocity,err,error,*999)
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & 1,1,nodeIdx2-1,3,velocity,err,error,*999)
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & 1,1,nodeIdx2,3,velocity,err,error,*999)
                  ! 
!!do the major part as moving average
                  !DO n1=1+offset,nodesInXi1-offset
                  !  nodeIdx2=nodeIdx2+1
                  !  velocity=0.0_DP
                  !  DO n4=nodeIdx2-offset,nodeIdx2+offset
                  !    CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,n4,3,dofIdx,err,error,*999)
                  !    CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !      & dofIdx,VALUE,err,error,*999)
                  !    velocity=velocity+VALUE
                  !  ENDDO !n4
                  !  velocity=velocity/(2*offset+1)
                  !  CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !    & 1,1,nodeIdx2,3,velocity,err,error,*999)
                  !ENDDO !n1
                  !
!!do the last three nodes of a fibre manually - arithmetic mean
                  !velocity=0.0_DP
                  !nodeIdx2=nodeIdx2+1
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE
                  !nodeIdx2=nodeIdx2+1
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE                  
                  !nodeIdx2=nodeIdx2+1
                  !CALL FieldVariable_LocalDOFIdxGet(monodomainIndependentVariable2,1,1,nodeIdx2,3,dofIdx,err,error,*999)
                  !CALL FieldVariable_ParameterSetGetLocalDOF(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & dofIdx,VALUE,err,error,*999)
                  !velocity=velocity+VALUE                  
                  !velocity=velocity/offset
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & 1,1,nodeIdx2-2,3,velocity,err,error,*999)                  
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_VALUES_SET_TYPE, &
                  !  & 1,1,nodeIdx2-1,3,velocity,err,error,*999)
                  !CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentVariable2,FIELD_U2_VARIABLE_TYPE, &
                  !  & FIELD_VALUES_SET_TYPE,1,1,nodeIdx2,3,velocity,err,error,*999)

                  !if there is not an adjacent element in positive XI_1 direction, go to the next FE element
                  CALL DecompositionElements_ElementNumberAdjacentGet(decompositionElements,ELEMENT_NORMAL_PLUS_XI1, &
                    & myElementIdx,numberOfAdjacentElements,err,error,*999)
                  IF(numberOfAdjacentElements==0) EXIT

                  !consider the adjacent element in positive XI_1 direction
                  CALL DecompositionElements_ElementAdjacentNumberGet(decompositionElements,1,ELEMENT_NORMAL_PLUS_XI1, &
                    &  myElementIdx,elementNumber,err,error,*999)

                  !if a fibre starts in the next element, go to the next FE elem
                  CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                    & elementNumber,4,startElem,err,error,*999)
                  !beginning of a fibre in this element: 1=yes, 0=no
                  IF (startElem==1) THEN
!!TODO: THIS IS JUST THE SAME AS ELEMENTIDX???
                    elementNumber=decompositionElements%ELEMENTS(elementIdx)%localNumber
                    EXIT
                  ENDIF

                  !find the elementIdx that corresponds to elementNumber
                  myElementIdx=0
                  DO elementIdx2=1,numberOfElements
                    CALL DecompositionElements_ElementAdjacentNumberGet(decompositionElements,1,0,elementIdx2, &
                      & elementNumber2,err,error,*999)
                    IF(elementNumber==elementNumber2) THEN
                      myElementIdx=elementIdx2
                      EXIT
                    ENDIF
                  ENDDO !elementIdx2
                  IF(myElementIdx==0) CALL FlagError("myElementIdx not found.",err,error,*999)                      

                  CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                    & elementNumber,1,nodesInXi1,err,error,*999)

                  startElem=0 !fibres don't start in this element

                  xi(1)=1.0_DP/(REAL(nodesInXi1))

                ENDDO !
                !for the beginning of the next fibre, go back to the element in which the last fibre started
                elementNumber=startElement

                CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                  & elementNumber,1,nodesInXi1,err,error,*999)

                myElementIdx=startElementIdx
                startElem=1 !fibres start in this element
                xi(1)=0.0_DP
                xi(2)=xi(2)+1.0_DP/(REAL(nodesInXi2))
              ENDDO !n2
              xi(1)=0.0_DP
              xi(2)=1.0_DP/(REAL(2*nodesInXi2))
              xi(3)=Xi(3)+1.0_DP/(REAL(nodesInXi3))
            ENDDO !n3

          ENDDO !elementIdx

        CASE(PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)

          !get the monodomain sub loop, solvers, solver, and finally geometric field and dependent field
          NULLIFY(monodomainControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,1,monodomainControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(monodomainControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(monodomainGeometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,monodomainGeometricField,err,error,*999)
          ! the Field_V_Variable_Type contains the 3D nodal positions
          NULLIFY(monodomainDependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,monodomainDependentField,err,error,*999)
          NULLIFY(monodomainIndependentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,monodomainIndependentField,err,error,*999)
          !get the finite elasticity sub loop, solvers, solver, and finally the dependent and independent fields
          NULLIFY(elasticityControlLoop)
          CALL ControlLoop_SubLoopGet(controlLoopParent,2,elasticityControlLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(elasticityControlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(elasticityGeometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,elasticityGeometricField,err,error,*999)
          NULLIFY(elasticityIndependentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,elasticityIndependentField,err,error,*999)
          !get the finite elasticity dependent field interpolation parameters of this element
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(equationsInterpolation)
          CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
          NULLIFY(interpolationParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
            & interpolationParameters,err,error,*999)
          NULLIFY(interpolatedPoint)
          CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,interpolatedPoint, &
            & err,error,*999)

          nodeIdx=0
          nodeIdx2=0
          fibreIdx=0
          CALL Field_VariableGet(monodomainDependentField,FIELD_V_VARIABLE_TYPE,monodomainDependentVariable,err,error,*999)
          CALL Field_VariableGet(monodomainGeometricField,FIELD_U_VARIABLE_TYPE,monodomainGeometricVariable,err,error,*999)
          CALL Field_VariableGet(elasticityIndependentField,FIELD_V_VARIABLE_TYPE,elasticityIndependentVariable,err,error,*999)

          NULLIFY(decomposition)
          CALL Field_DecompositionGet(elasticityGeometricField,decomposition,err,error,*999)
          NULLIFY(decompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
          NULLIFY(decompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
          NULLIFY(domain)
          CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)

          !loop over the elements of the finite elasticity mesh (internal and boundary elements)
          !no need to consider ghost elements here since only bioelectrical fields are changed
          DO elementIdx=1,numberOfElements

!!TODO: IS ELEMENTNUMBER NOT JUST THE SAME AS ELEMENTIDX???
            elementNumber=decompositionElements%ELEMENTS(elementIdx)%localNumber
            myElementIdx=elementIdx

            !the Field_V_Variable_Type of the FE independent field contains the number of nodes in each Xi-direction
            !of the bioelectrics grid
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
              & elementNumber,1,nodesInXi1,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
              & elementNumber,2,nodesInXi2,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
              & elementNumber,3,nodesInXi3,err,error,*999)
            !beginning of a fibre in this element: 1=yes, 0=no
            CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
              & elementNumber,4,startElem,err,error,*999)

            !if there is no bioelectrics grid in this finite elasticity element, or the fibres don't begin in this element,
            !jump to the next element
            IF((nodesInXi1==0).OR.(nodesInXi2==0).OR.(nodesInXi3==0).OR.(startElem==0)) CYCLE

            !get the finite elasticity dependent field interpolation parameters of this element
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,interpolationParameters, &
              & err,error,*999)
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
            NULLIFY(quadratureScheme)
            CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)

            startElement=elementNumber
            startElementIdx=myElementIdx

            !assume Xi(1) to be normal to the seed surface, i.e. the seed points have Xi(1)=0
            xi=[0.0_DP,1.0_DP/(REAL(2*nodesInXi2)),1.0_DP/(REAL(2*nodesInXi3))]

            !assume that the bioelectrics node numbers are increased in order Xi(1), Xi(2), Xi(3) 
            DO n3=1,nodesInXi3
              DO n2=1,nodesInXi2
                fibreIdx=fibreIdx+1

                !loop over the FE elements that contain nodes of the very same fibres
                DO                                    
                  DO n1=1,nodesInXi1
                    nodeIdx=nodeIdx+1

                    !store the fibre number this bioelectrics node belongs to.
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,3,fibreIdx,err,error,*999) 

                    !find the interpolated position of the bioelectric grid node from the FE dependent field
                    CALL Field_InterpolateXi(NO_PART_DERIV,XI,interpolatedPoint,err,error,*999)
                    !update the bioelectrics dependent field Field_V_Variable_Type
                    !the Field_V_Variable_Type of the monodomain dependent field contains the nodal positions in 3D
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,1,interpolatedPoint%values(1,NO_PART_DERIV),err,error,*999)
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,2,interpolatedPoint%values(2,NO_PART_DERIV),err,error,*999)
                    CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                      & 1,1,nodeIdx,3,interpolatedPoint%values(3,NO_PART_DERIV),err,error,*999)

                    IF((n1==1).AND.(elementNumber==startElement)) THEN
                      !a new line of bioelectrics grid nodes begins
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,0.0_DP,err,error,*999)
                    ELSE
                      !get the position in 3D of the previous node
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,1,previousNode(1),err,error,*999)
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,2,previousNode(2),err,error,*999)
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainDependentVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,3,previousNode(3),err,error,*999)

                      !compute the distance between the previous node and the actual node
                      distance=SQRT((interpolatedPoint%values(1,NO_PART_DERIV)- &
                        & previousNode(1))*(interpolatedPoint%values(1,NO_PART_DERIV)-previousNode(1))+ &
                        & (interpolatedPoint%values(2,NO_PART_DERIV)- &
                        & previousNode(2))*(interpolatedPoint%values(2,NO_PART_DERIV)-previousNode(2))+ &
                        & (interpolatedPoint%values(3,NO_PART_DERIV)- &
                        & previousNode(3))*(interpolatedPoint%values(3,NO_PART_DERIV)-previousNode(3)))

                      !get the position in 1D of the previous node
                      CALL FieldVariable_ParameterSetGetLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx-1,1,valueLeft,err,error,*999)
                      !update the current 1D node position
                      CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainGeometricVariable,FIELD_VALUES_SET_TYPE, &
                        & 1,1,nodeIdx,1,valueLeft+distance,err,error,*999)

                    ENDIF !((n1==1).AND.(elementNumber==startElement))

                    IF(calcClosestGaussPoint) THEN
                      !calculate the closest finite elasticity Gauss point of each bioelectrics node
                      distance=1000000.0_DP
                      gaussPoint=0
                      DO gaussPointIdx=1,numberOfGauss
                        !compute the distance between the bioelectrics node and the Gauss point
                        CALL BasisQuadratureScheme_GaussPositionGet(quadratureScheme,gaussPointIdx,gaussPosition,err,error,*999)
                        VALUE=SQRT((Xi(1)-gaussPosition(1))*(Xi(1)-gaussPosition(1))+ &
                          & (Xi(2)-gaussPosition(2))*(Xi(2)-gaussPosition(2))+ &
                          & (Xi(3)-gaussPosition(3))*(Xi(3)-gaussPosition(3)))
                        IF(value<distance) THEN
                          distance=VALUE
                          gaussPoint=gaussPointIdx
                        ENDIF
                      ENDDO !gaussPointIdx
                      IF(gaussPoint==0) CALL FlagWarning("Closest Gauss Point not found",err,error,*999)
                      !store the nearest Gauss Point info and the inElement info (local element number!!!)
                      CALL Field_ParameterSetUpdateLocalNode(monodomainIndependentField,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4,gaussPoint,err,error,*999)
                      CALL Field_ParameterSetUpdateLocalNode(monodomainIndependentField,FIELD_V_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5,elementNumber,err,error,*999)
                    ENDIF !calcClosestGaussPoint

                    IF(startElem==1) THEN
                      !fibres start in this element
                      xi(1)=xi(1)+1.0_DP/(REAL(nodesInXi1-1))
                    ELSEIF(startElem==0) THEN
                      !fibres don't start in this element
                      xi(1)=xi(1)+1.0_DP/(REAL(nodesInXi1))
                    ELSE
                      localError="The start element index is incorrect. The index is "// &
                        & TRIM(NumberToVString(startElem,"*",err,error))//" and should be zero or one." 
                      CALL FlagError(localError,err,error,*999)
                    ENDIF

                  ENDDO !n1

                  !if there is not an adjacent element in positive XI_1 direction, go to the next FE element
                  CALL DecompositionElements_ElementNumberAdjacentGet(decompositionElements,ELEMENT_NORMAL_PLUS_XI1, &
                    & myElementIdx,numberOfAdjacentElements,err,error,*999)
                  IF(numberOfAdjacentElements==0) EXIT

                  !consider the adjacent element in positive XI_1 direction
                  CALL DecompositionElements_ElementAdjacentNumberGet(decompositionElements,1,ELEMENT_NORMAL_PLUS_XI1, &
                    &  myElementIdx,elementNumber,err,error,*999)

                  !if a fibre starts in the next element, go to the next FE elem
                  CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                    & elementNumber,4,startElem,err,error,*999)
                  !beginning of a fibre in this element: 1=yes, 0=no
                  IF (startElem==1) THEN
!!TODO: THIS IS JUST THE SAME AS ELEMENTIDX???
                    elementNumber=decompositionElements%elements(elementIdx)%localNumber
                    EXIT
                  ENDIF

                  !find the elementIdx that corresponds to elementNumber
                  myElementIdx=0
                  CALL DecompositionElements_NumberOfElementsGet(decompositionElements,numberOfElements,err,error,*999)
                  DO elementIdx2=1,numberOfElements
                    CALL DecompositionElements_ElementAdjacentNumberGet(decompositionElements,1,0,elementIdx2, &
                      & elementNumber2,err,error,*999)
                    IF(elementNumber==elementNumber2) THEN
                      myElementIdx=elementIdx2
                      EXIT
                    ENDIF
                  ENDDO
                  IF(myElementIdx==0) CALL FlagError("myElementIdx not found.",err,error,*999)                      

                  CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                    & elementNumber,1,nodesInXi1,err,error,*999)

                  startElem=0 !fibres don't start in this element

                  xi(1)=1.0_DP/(REAL(nodesInXi1))

                ENDDO !
                !for the beginning of the next fibre, go back to the element in which the last fibre started
                elementNumber=startElement

                CALL FieldVariable_ParameterSetGetLocalElement(elasticityIndependentVariable,FIELD_VALUES_SET_TYPE, &
                  & elementNumber,1,nodesInXi1,err,error,*999)

                myElementIdx=startElementIdx
                startElem=1 !fibres start in this element
                xi(1)=0.0_DP
                xi(2)=xi(2)+1.0_DP/(REAL(nodesInXi2))
              ENDDO !n2
              xi(1)=0.0_DP
              xi(2)=1.0_DP/(REAL(2*nodesInXi2))
              xi(3)=Xi(3)+1.0_DP/(REAL(nodesInXi3))
            ENDDO !n3

          ENDDO !elementIdx

        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
            & " is not valid for a bioelectrics finite elasticity of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second problem specification of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for a multi physics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      !the main time loop - do nothing!
    ENDIF
    
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_UpdateGeometricField",err,error)
    EXITS("BioelectricFiniteElasticity_UpdateGeometricField")
    RETURN 1
    
  END SUBROUTINE BioelectricFiniteElasticity_UpdateGeometricField

  !
  !================================================================================================================================
  !

  !>Interpolates the finite elasticity independent field from the biolectrics independent field.
  !>NOTE: this is only temporary - will be replaced once embedded meshes are available
  SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryFinish,dofIdx,elementIdx,elementNumber,gaussPointIdx,inElement,internalStart,nearestGP,nodeIdx, &
      & numberOfGauss,numberOfLocal,numberOfSubLoops,pSpecification(3)
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_GAUSS_POINTS=64
    INTEGER(INTG) :: numberOfNodes(MAX_NUMBER_OF_GAUSS_POINTS)
    REAL(DP) :: a1,a1Values(MAX_NUMBER_OF_GAUSS_POINTS),a2,a2Values(MAX_NUMBER_OF_GAUSS_POINTS),activation, &
      & activationValues(MAX_NUMBER_OF_GAUSS_POINTS),activeStress,activeStressValues(MAX_NUMBER_OF_GAUSS_POINTS), &
      & titinStressBound,titinStressCrossFibreBound,titinStressCrossFibreUnbound,titinStressUnbound, &
      & titinStressValuesBound(MAX_NUMBER_OF_GAUSS_POINTS),titinStressValuesUnbound(MAX_NUMBER_OF_GAUSS_POINTS), &
      & titinStressValuesCrossFibreBound(MAX_NUMBER_OF_GAUSS_POINTS), &
      & titinStressValuesCrossFibreUnbound(MAX_NUMBER_OF_GAUSS_POINTS),x1,x1Values(MAX_NUMBER_OF_GAUSS_POINTS),x2, &
      & x2Values(MAX_NUMBER_OF_GAUSS_POINTS)
    TYPE(BasisType), POINTER :: basis
    TYPE(ControlLoopType), POINTER :: controlLoopRoot,controlLoopParent,elasticityControlLoop,monodomainControlLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: elementsMapping,nodesMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: monodomainIndependentField,elasticityIndependentField
    TYPE(FieldVariableType), POINTER :: fieldUVariable,fieldVVariable,fieldFEVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("BioelectricFiniteElasticity_IndependentFieldInterpolate",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
      & PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE,PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
      CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
      IF(numberOfSubLoops==0) THEN
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoopParent)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoopParent,err,error,*999)
        !--- MONODOMAIN ---
        NULLIFY(monodomainControlLoop)
        CALL ControlLoop_SubLoopGet(controlLoopParent,1,monodomainControlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(monodomainControlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(monodomainIndependentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,monodomainIndependentField,err,error,*999)

        !--- FINITE ELASTICITY ---
        NULLIFY(elasticityControlLoop)
        CALL ControlLoop_SubLoopGet(controlLoopParent,2,elasticityControlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(elasticityControlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(elasticityIndependentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,elasticityIndependentField,err,error,*999)

        !--- NOW INTERPOLATE ---
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(elasticityIndependentField,decomposition,err,error,*999)
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
        NULLIFY(domainMappings)
        CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
        NULLIFY(elementsMapping)
        CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(monodomainIndependentField,decomposition,err,error,*999)
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
        NULLIFY(domainMappings)
        CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
        NULLIFY(nodesMapping)
        CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)
        
        CALL Field_VariableGet(monodomainIndependentField,FIELD_U_VARIABLE_TYPE,fieldUVariable,err,error,*999)
        CALL Field_VariableGet(monodomainIndependentField,FIELD_V_VARIABLE_TYPE,fieldVVariable,err,error,*999)
        CALL Field_VariableGet(elasticityIndependentField,FIELD_U_VARIABLE_TYPE,fieldFEVariable,err,error,*999)
        
        !loop over the finite elasticity elements
        !first process the internal and boundary elements
        CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
        CALL DomainMapping_BoundaryFinishGet(elementsMapping,boundaryFinish,err,error,*999)
        DO elementIdx=internalStart,boundaryFinish
          CALL DomainMapping_NumberGet(elementsMapping,elementIdx,elementNumber,err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
          NULLIFY(quadratureScheme)
          CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
              
          IF(numberOfGauss>MAX_NUMBER_OF_GAUSS_POINTS) THEN
            localError="The number of Gauss points of "//TRIM(NumberToVString(numberOfGauss,"*",err,error))// &
              & " is greater than the allowed maximum number of Gauss points of "// &
              & TRIM(NumberToVString(MAX_NUMBER_OF_GAUSS_POINTS,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          
          numberOfNodes=0
          activeStressValues=0.0_DP
          titinStressValuesUnbound=0.0_DP
          titinStressValuesBound=0.0_DP
          titinStressValuesCrossFibreUnbound=0.0_DP
          titinStressValuesCrossFibreBound=0.0_DP
          activationValues=0.0_DP
          a1Values=0.0_DP
          a2Values=0.0_DP
          x1Values=0.0_DP
          x2Values=0.0_DP
          
          !loop over the bioelectrics nodes
          CALL DomainMapping_NumberOfLocalGet(nodesMapping,numberOfLocal,err,error,*999)
          DO nodeIdx=1,numberOfLocal
            !component 5 of variable V contains inElem info (LOCAL NUMBERING!!!)
            CALL FieldVariable_ParameterSetGetLocalNode(fieldVVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5,inElement, &
              & err,error,*999) 
            
            !check if the bioelectrics node is located within the finite elasticity element
            IF(inElement==elementNumber) THEN
              !component 4 of variable V contains Nearest Gauss Point info
              CALL FieldVariable_ParameterSetGetLocalNode(fieldVVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4,nearestGP, &
                & err,error,*999)
              IF(nearestGP>MAX_NUMBER_OF_GAUSS_POINTS) THEN
                localError="The nearest Gauss point number of "//TRIM(NumberToVString(nearestGP,"*",err,error))// &
                  & " is greater than the allowed maximum number of Gauss points of "// &
                  & TRIM(NumberToVString(MAX_NUMBER_OF_GAUSS_POINTS,"*",err,error))//"."
              ENDIF
              
              SELECT CASE(pSpecification(3))
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
                
                !component 1 of variable U contains the active stress
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1,activeStress, &
                  & err,error,*999)
                
                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                numberOfNodes(nearestGP)=numberOfNodes(nearestGP)+1
                !add up the active stress value
                activeStressValues(nearestGP)=activeStressValues(nearestGP)+activeStress
                
              CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
                
                !component 1 of variable U contains the active stress
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1,activeStress, &
                  & err,error,*999)
                    
                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                numberOfNodes(nearestGP)=numberOfNodes(nearestGP)+1
                !add up the active stress value
                activeStressValues(nearestGP)=activeStressValues(nearestGP)+activeStress
                
                !component 2 of variable U contains the titin stress unbound
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,2, &
                  & titinStressUnbound,err,error,*999)
                !component 3 of variable U contains the titin stress bound
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,3, &
                  & titinStressBound,err,error,*999)
                !component 4 of variable U contains the titin XF-stress (cross-fibre directions) unbound
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4, &
                  & titinStressCrossFibreUnbound,err,error,*999)
                !component 5 of variable U contains the titin XF-stress (cross-fibre directions) bound
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5, &
                  & titinStressCrossFibreBound,err,error,*999)
                !component 6 of variable U contains the titin activation
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,6, &
                  & activation,err,error,*999)
                
                titinStressValuesUnbound(nearestGP)=titinStressValuesUnbound(nearestGP)+titinStressUnbound
                titinStressValuesBound(nearestGP)=titinStressValuesBound(nearestGP)+titinStressBound
                titinStressValuesCrossFibreUnbound(nearestGP)=titinStressValuesCrossFibreUnbound(nearestGP) + &
                  & titinStressCrossFibreUnbound
                titinStressValuesCrossFibreBound(nearestGP)=titinStressValuesCrossFibreBound(nearestGP) + &
                  & titinStressCrossFibreBound
                activationValues(nearestGP)=activationValues(nearestGP)+activation
                
              CASE(PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
                
                !count the number of bioelectrics nodes that are closest to each finite elasticity Gauss point
                numberOfNodes(nearestGP)=numberOfNodes(nearestGP)+1
                
                !component 1 of variable U contains a1
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1,a1,err,error,*999)
                !component 2 of variable U contains a2
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,2,a2,err,error,*999)
                !component 3 of variable U contains x1
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,3,x1,err,error,*999)
                !component 4 of variable U contains x2
                CALL FieldVariable_ParameterSetGetLocalNode(fieldUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4,x2,err,error,*999)
                
                a1Values(nearestGP)=a1Values(nearestGP)+a1
                a2Values(nearestGP)=a2Values(nearestGP)+a2
                x1Values(nearestGP)=x1Values(nearestGP)+x1
                x2Values(nearestGP)=x2Values(nearestGP)+x2
              END SELECT
              
            ENDIF
          ENDDO !nodeIdx

          !loop over the finite elasticity Gauss points
          DO gaussPointIdx=1,numberOfGauss
            !make sure we don't divide by zero
            IF(numberOfNodes(gaussPointIdx)<=0) THEN
              activeStress=0.0_DP
              titinStressUnbound=0.0_DP
              titinStressBound=0.0_DP
              titinStressCrossFibreUnbound=0.0_DP
              titinStressCrossFibreBound=0.0_DP
              activation=0.0_DP
              a1=0.0_DP
              a2=0.0_DP
              x1=0.0_DP
              x2=0.0_DP
            ELSE
              activeStress=activeStressValues(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              titinStressUnbound=titinStressValuesUnbound(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              titinStressBound=titinStressValuesBound(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              titinStressCrossFibreUnbound=titinStressValuesCrossFibreUnbound(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              titinStressCrossFibreBound=titinStressValuesCrossFibreBound(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              activation=activationValues(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              a1=a1Values(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              a2=a2Values(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              x1=x1Values(gaussPointIdx)/numberOfNodes(gaussPointIdx)
              x2=x2Values(gaussPointIdx)/numberOfNodes(gaussPointIdx)
            ENDIF

            SELECT CASE(pSpecification(3))
            CASE(PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE,PROBLEM_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)

              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,1,activeStress,err,error,*999)
              
            CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
              
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,1,activeStress,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,2,titinStressUnbound,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,3,titinStressBound,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,4,titinStressCrossFibreUnbound,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,5,titinStressCrossFibreBound,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,6,activation,err,error,*999)
              
            CASE(PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
              
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,1,a1,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,2,a2,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,3,x1,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(fieldFEVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
                & elementNumber,4,x2,err,error,*999)
              
            END SELECT
            
          ENDDO !gaussPointIdx
        ENDDO !elementIdx

        !now the ghost elements -- get the relevant info from the other computational nodes
        CALL FieldVariable_ParameterSetUpdateStart(fieldFEVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fieldFEVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        
      ENDIF
    CASE DEFAULT
      localError="Independent field interpolation is not implemented for problem subtype " &
        & //TRIM(NumberToVString(pSpecification(3),"*",err,error))
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN
999 ERRORS("BioelectricFiniteElasticity_IndependentFieldInterpolate",err,error)
    EXITS("BioelectricFiniteElasticity_IndependentFieldInterpolate")
    RETURN 1

  END SUBROUTINE BioelectricFiniteElasticity_IndependentFieldInterpolate

  !
  !================================================================================================================================
  !

  !>Computes force enhancement based on the titin model of C Rode et al. (2009). Force depression is not yet implemented.
  !>Titin-induced force enhancement and force depression: A 'sticky-spring' mechanism in muscle contractions? C Rode, T Siebert, R Blickhan - Journal of theoretical biology, 2009.
  SUBROUTINE BioelectricFiniteElasticity_ComputeTitin(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !A pointer to the time control loop
    INTEGER(INTG), INTENT(OUT) :: err !The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,indexI,indexPseudo,indexRef,nodeIdx,numberOfLocal,numberOfSubLoops,pSpecification(3),switchModel
    INTEGER(INTG), PARAMETER :: DIM_DATA=250
    REAL(DP) :: actinMyosinDistance,d10,deltaF,diffQuot,elongation,elongationDistIG,elongationNew,elongationPEVK,f0,force, &
      & forceDistalIG,forcesDistIG(250),lengthDistIGF0,lengthDistIG,lengthsDistIG(250),lengthInitTitin,lengthTitin,sarcoLength, &
      & sarcoLengthAtActivation,slope,stiffnessDist,stiffnessPEVK,titinBound,titinUnbound,titinXFBound,titinXFUnbound
    REAL(DP), PARAMETER :: LENGTH_ACTIN=1.04_DP,LENGTH_MBAND=0.0625_DP,LENGTH_MYOSIN=0.7375_DP
    REAL(DP), PARAMETER :: LENGTH_ZERO=0.635_DP,LENGTH_ZDISC=0.05_DP
    REAL(DP), PARAMETER :: DX=0.001_DP
    REAL(DP), PARAMETER :: FORCE_INCREMENT=1.e-5_DP,TOL=1.e-5_DP
    REAL(DP), PARAMETER, DIMENSION(5) :: COEFF_MATRIX=[5.0239_DP,-0.6717_DP,-2.5841_DP,-5.0128_DP,-5.0239_DP]
    TYPE(ControlLoopType), POINTER :: controlLoopRoot,controlLoopParent,monodomainControlLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: nodesMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(FieldType), POINTER :: monodomainIndependentField
    TYPE(FieldVariableType), POINTER :: monodomainIndependentUVariable,monodomainIndependentU1Variable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    
    ENTERS("BioelectricFiniteElasticity_ComputeTitin",err,error,*999)

    !the realtion between length and force of the distal_Ig region is very nonlinear. Linear interpolation of Rode's data is used here.
    forcesDistIG= &
      & [0.0_DP,0.001_DP,0.002_DP,0.003_DP,0.004_DP,0.005_DP,0.006_DP,0.007_DP,0.008_DP,0.009_DP,0.01_DP,0.011_DP,0.012_DP, &
      & 0.013_DP,0.014_DP,0.015_DP,0.016_DP,0.017_DP,0.018_DP,0.019_DP,0.02_DP,0.021_DP,0.022_DP,0.023_DP,0.024_DP, &
      & 0.025_DP,0.026_DP,0.027_DP,0.028_DP,0.029_DP,0.03_DP,0.031_DP,0.032_DP,0.033_DP,0.034_DP,0.035_DP,0.036_DP, &
      & 0.037_DP,0.038_DP,0.039_DP,0.04_DP,0.041_DP,0.042_DP,0.043_DP,0.044_DP,0.045_DP,0.046_DP,0.047_DP,0.048_DP, &
      & 0.049_DP,0.05_DP,0.051_DP,0.052_DP,0.053_DP,0.054_DP,0.055_DP,0.056_DP,0.057_DP,0.058_DP,0.059_DP,0.06_DP, &
      & 0.061_DP,0.062_DP,0.063_DP,0.064_DP,0.065_DP,0.066_DP,0.067_DP,0.068_DP,0.069_DP,0.070_DP,0.071_DP,0.072_DP, &
      & 0.073_DP,0.074_DP,0.075_DP,0.076_DP,0.077_DP,0.078_DP,0.079_DP,0.080_DP,0.081_DP,0.082_DP,0.083_DP,0.084_DP, &
      & 0.085_DP,0.086_DP,0.087_DP,0.088_DP,0.089_DP,0.090_DP,0.091_DP,0.092_DP,0.093_DP,0.094_DP,0.095_DP,0.096_DP, &
      & 0.097_DP,0.098_DP,0.099_DP,0.1_DP,0.101_DP,0.102_DP,0.103_DP,0.104_DP,0.105_DP,0.106_DP,0.107_DP,0.108_DP, &
      & 0.109_DP,0.11_DP,0.111_DP,0.112_DP,0.113_DP,0.114_DP,0.115_DP,0.116_DP,0.117_DP,0.118_DP,0.119_DP,0.12_DP, &
      & 0.121_DP,0.122_DP,0.123_DP,0.124_DP,0.125_DP,0.126_DP,0.127_DP,0.128_DP,0.129_DP,0.13_DP,0.131_DP,0.132_DP, &
      & 0.133_DP,0.134_DP,0.135_DP,0.136_DP,0.137_DP,0.138_DP,0.139_DP,0.14_DP,0.141_DP,0.142_DP,0.143_DP,0.144_DP, &
      & 0.145_DP,0.146_DP,0.147_DP,0.148_DP,0.149_DP,0.15_DP,0.151_DP,0.152_DP,0.153_DP,0.154_DP,0.155_DP,0.156_DP, &
      & 0.157_DP,0.158_DP,0.159_DP,0.16_DP,0.161_DP,0.162_DP,0.163_DP,0.164_DP,0.165_DP,0.166_DP,0.167_DP, 0.168_DP, &
      & 0.169_DP,0.17_DP,0.171_DP,0.172_DP,0.173_DP,0.174_DP,0.175_DP,0.176_DP,0.177_DP,0.178_DP,0.179_DP,0.18_DP, &
      & 0.181_DP,0.182_DP,0.183_DP,0.184_DP,0.185_DP,0.186_DP,0.187_DP,0.188_DP,0.189_DP,0.19_DP,0.191_DP,0.192_DP, &
      & 0.193_DP,0.194_DP,0.195_DP,0.196_DP,0.197_DP,0.198_DP,0.199_DP,0.2_DP,0.201_DP,0.202_DP,0.203_DP,0.204_DP, &
      & 0.205_DP,0.206_DP,0.207_DP,0.208_DP,0.209_DP,0.21_DP,0.211_DP,0.212_DP,0.213_DP,0.214_DP,0.215_DP,0.216_DP, &
      & 0.217_DP,0.218_DP,0.219_DP,0.22_DP,0.221_DP,0.222_DP,0.223_DP,0.224_DP,0.225_DP,0.226_DP,0.227_DP,0.228_DP, &
      & 0.229_DP,0.23_DP,0.231_DP,0.232_DP,0.233_DP,0.234_DP,0.235_DP,0.236_DP,0.237_DP,0.238_DP,0.239_DP,0.24_DP, &
      & 0.241_DP,0.242_DP,0.243_DP,0.244_DP,0.245_DP,0.246_DP,0.247_DP,0.248_DP,0.249_DP]
    lengthsDistIG= &
      & [0.0_DP,0.03461753545561_DP,0.049729169766010_DP,0.058506390323323_DP,0.064606296848594_DP,0.06922519775133_DP, &
      & 0.0729080120998386_DP,0.0759458896446241_DP,0.0785230355395668_DP,0.0807314335143191_DP,0.0826674161660979_DP, &
      & 0.0843819721302_DP,0.0859161360822_DP,0.087300738288_DP,0.0885510536196_DP,0.08970061165_DP,0.090751366_DP, &
      & 0.0917285714001_DP,0.0926271710799_DP,0.093467026018_DP,0.094254010845_DP,0.094992014919_DP,0.09568255451_DP, &
      & 0.0963346932312_DP,0.0969518155718_DP,0.097537000419_DP,0.098093047899_DP,0.098622503780_DP,0.09912768169_DP, &
      & 0.0996103583026_DP,0.1000734170008_DP,0.100517614158_DP,0.100942907967_DP,0.101351270601_DP,0.101745244913_DP, &
      & 0.1021260518375_DP,0.1024947934706_DP,0.102851281365_DP,0.103194317381_DP,0.103528086731_DP,0.103853341524_DP, &
      & 0.1041686065999_DP,0.1044739635997_DP,0.104772842609_DP,0.105064758806_DP,0.105347078318_DP,0.105624524436_DP, &
      & 0.1058959740402_DP,0.1061596374590_DP,0.106419674124_DP,0.106673222977_DP,0.106921779223_DP,0.107167251003_DP, &
      & 0.1074055929211_DP,0.1076418938850_DP,0.107872798106_DP,0.108100580266_DP,0.108324935479_DP,0.108545154704_DP, &
      & 0.1087634145225_DP,0.1089769308376_DP,0.109189523632_DP,0.109397109133_DP,0.109604335645_DP,0.109806785604_DP, &
      & 0.1100090293896_DP,0.1102069598668_DP,0.110404790711_DP,0.110598542577_DP,0.110792294444_DP,0.110982362315_DP, &
      & 0.1111722575399_DP,0.1113591719379_DP,0.111545514634_DP,0.111729654458_DP,0.111912731875_DP,0.112094428489_DP, &
      & 0.1122745119339_DP,0.1124540532647_DP,0.112631398979_DP,0.112808744694_DP,0.112983883284_DP,0.113158733268_DP, &
      & 0.1133324054841_DP,0.1135049882670_DP,0.113677360515_DP,0.113847891889_DP,0.114018423263_DP,0.114187784974_DP, &
      & 0.1143564686814_DP,0.1145249703398_DP,0.114691998719_DP,0.114859027099_DP,0.115025270473_DP,0.115190825075_DP, &
      & 0.1153563796784_DP,0.1155207617063_DP,0.115685013868_DP,0.115849024045_DP,0.116012135436_DP,0.116175246828_DP, &
      & 0.1163378963393_DP,0.1165000194746_DP,0.116662142610_DP,0.116823704015_DP,0.116984982740_DP,0.117146261465_DP, &
      & 0.1173069696799_DP,0.1174675396300_DP,0.117628109580_DP,0.117788165045_DP,0.117948154078_DP,0.118108143111_DP, &
      & 0.1182677147299_DP,0.1184272433357_DP,0.118586771941_DP,0.118745999791_DP,0.118905181478_DP,0.119064363164_DP, &
      & 0.1192233610196_DP,0.1193823026790_DP,0.119541244338_DP,0.119700101992_DP,0.119858904246_DP,0.120017706500_DP, &
      & 0.1201764919257_DP,0.1203352494533_DP,0.120494006981_DP,0.120652768322_DP,0.120811570168_DP,0.120970372014_DP, &
      & 0.1211291738603_DP,0.1212880693042_DP,0.121446999172_DP,0.121605929041_DP,0.121764923098_DP,0.121924059629_DP, &
      & 0.1220831961613_DP,0.1222423326927_DP,0.122401700265_DP,0.122561117298_DP,0.122720534331_DP,0.122880045624_DP, &
      & 0.1230398124447_DP,0.1231995792652_DP,0.123359346085_DP,0.123519381317_DP,0.123679562894_DP,0.123839744470_DP, &
      & 0.1239999260467_DP,0.1241605622478_DP,0.124321219453_DP,0.124481876659_DP,0.124642637543_DP,0.124803827368_DP, &
      & 0.1249650171945_DP,0.1251262070202_DP,0.125287609380_DP,0.125449385133_DP,0.125611160886_DP,0.125772936638_DP, &
      & 0.1259350043157_DP,0.1260974158090_DP,0.126259827302_DP,0.126422238795_DP,0.126584979901_DP,0.126748073636_DP, &
      & 0.1269111673704_DP,0.1270742611047_DP,0.127237669513_DP,0.127401488846_DP,0.127565308178_DP,0.127729127511_DP, &
      & 0.1278931842348_DP,0.1280577695422_DP,0.128222354849_DP,0.128386940157_DP,0.128551614623_DP,0.128717003455_DP, &
      & 0.1288823922862_DP,0.1290477811173_DP,0.129213169948_DP,0.129379259581_DP,0.129545486803_DP,0.129711714025_DP, &
      & 0.1298779412472_DP,0.1300445897140_DP,0.130211687650_DP,0.130378785586_DP,0.130545883522_DP,0.130713029866_DP, &
      & 0.1308810284274_DP,0.1310490269884_DP,0.131217025549_DP,0.131385024110_DP,0.131553528414_DP,0.131722455223_DP, &
      & 0.1318913820322_DP,0.1320603088411_DP,0.132229235650_DP,0.132399074312_DP,0.132568954822_DP,0.132738835332_DP, &
      & 0.1329087158420_DP,0.1330788764544_DP,0.133249734060_DP,0.133420591667_DP,0.133591449273_DP,0.133762306879_DP, &
      & 0.1339336992411_DP,0.1341055553883_DP,0.134277411535_DP,0.134449267682_DP,0.134621123829_DP,0.134793694483_DP, &
      & 0.1349665687658_DP,0.1351394430487_DP,0.135312317331_DP,0.135485191614_DP,0.135658878561_DP,0.135832788820_DP, &
      & 0.1360066990800_DP,0.1361806093395_DP,0.136354519599_DP,0.136529253243_DP,0.136704215658_DP,0.136879178072_DP, &
      & 0.1370541404878_DP,0.1372291029027_DP,0.137404806924_DP,0.137580836098_DP,0.137756865271_DP,0.137932894445_DP, &
      & 0.1381089236188_DP,0.1382855157880_DP,0.138462624830_DP,0.138639733872_DP,0.138816842915_DP,0.138993951957_DP, &
      & 0.1391713448843_DP,0.1393495454909_DP,0.139527746097_DP,0.139705946704_DP,0.139884147310_DP,0.140062347917_DP, &
      & 0.1402415516727_DP,0.1404208541990_DP,0.140600156725270_DP,0.140779459251523_DP,0.140958761777775_DP]

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE)
      CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
      IF(numberOfSubLoops==0) THEN
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoopParent)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoopParent,err,error,*999)
        !The first control_loop is the one for monodomain
        NULLIFY(monodomainControlLoop)
        CALL ControlLoop_SubLoopGet(controlLoopParent,1,monodomainControlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(monodomainControlLoop,solvers,err,error,*999)
        !The second solver is associated with the diffusion part of the monodomain equation
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(monodomainIndependentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,monodomainIndependentField,err,error,*999)

        NULLIFY(monodomainIndependentUVariable)
        CALL Field_VariableGet(monodomainIndependentField,FIELD_U_VARIABLE_TYPE,monodomainIndependentUVariable,err,error,*999)
        NULLIFY(monodomainIndependentU1Variable)
        CALL Field_VariableGet(monodomainIndependentField,FIELD_U1_VARIABLE_TYPE,monodomainIndependentU1Variable,err,error,*999)

        NULLIFY(decomposition)
        CALL Field_DecompositionGet(monodomainIndependentField,decomposition,err,error,*999)
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
        NULLIFY(domainMappings)
        CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
        NULLIFY(nodesMapping)
        CALL DomainMappings_NodesMappingGet(domainMappings,nodesMapping,err,error,*999)
        
        ! Initialization
        indexRef=1

        CALL DomainMapping_NumberOfLocalGet(nodesMapping,numberOfLocal,err,error,*999)
        DO nodeIdx=1,numberOfLocal
          
          !the fourth component of the U1 variable contains the half sarcomere length at activation
          CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentU1Variable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4, &
            & sarcoLengthAtActivation,err,error,*999)
          
          !the first component of the U1 variable contains the actual half sarcomere length
          CALL FieldVariable_ParameterSetGetLocalNode(monodomainIndependentU1Variable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1, &
            & sarcoLength,err,error,*999)
          
          elongation=sarcoLength-sarcoLengthAtActivation
          
          IF(elongation<ZERO_TOLERANCE) THEN
            lengthTitin=sarcoLength-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
            !function to approximate the relation between the initial titin length and the initial passive force F0	
            titinUnbound=COEFF_MATRIX(1)*EXP(lengthTitin)+COEFF_MATRIX(2)*lengthTitin**3+COEFF_MATRIX(3)* &
              & lengthTitin**2+COEFF_MATRIX(4)*lengthTitin+COEFF_MATRIX(5) 
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,2, &
              & titinUnbound,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,3, &
              & titinUnbound,err,error,*999)
            
            ! calculate x-fibre titin stress (trigonometry):                      
            ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
            d10=-13.39_DP*sarcoLength+58.37_DP
            ! d10=-13.39_DP*1.0_DP+58.37_DP
            ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
            actinMyosinDistance=0.001_DP*(2.0_DP/3.0_DP*d10)  
            ! calculate x-fibre stress with tangens-function (unbound titin)
            titinXFUnbound=0.5_DP*titinUnbound*actinMyosinDistance/lengthTitin
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4, &
              & titinXFUnbound,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5, &
              & titinXFUnbound,err,error,*999)
            
          ELSE  ! Force enhancement  
            
            ! Calculate Titin-Force for unbound titin filaments
            lengthTitin=sarcoLength-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
            !function to approximate the relation between the initial titin length and the initial passive force F0	
            titinUnbound=COEFF_MATRIX(1)*EXP(lengthTitin)+COEFF_MATRIX(2)*lengthTitin**3+COEFF_MATRIX(3)* &
              & lengthTitin**2+COEFF_MATRIX(4)*lengthTitin+COEFF_MATRIX(5)
            
            ! Calculate Titin-Force for bound titin filaments
            ! Switch between different force-enhancent models/implementations 
            ! 1=linear approximation 2=square approximation 3=simple iterative solution of Rode's model
            ! 4=Newton's-method to solve Rode's model
            switchModel=4
            !linear approximation
            IF(switchModel==1) THEN
              lengthInitTitin=sarcoLengthAtActivation-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
              !function to approximate the relation between the initial titin length and the initial passive force F0	
              F0=COEFF_MATRIX(1)*EXP(lengthInitTitin)+COEFF_MATRIX(2)*lengthInitTitin**3+COEFF_MATRIX(3)* &
                & lengthInitTitin**2+COEFF_MATRIX(4)*lengthInitTitin+COEFF_MATRIX(5)
              !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
              stiffnessPEVK=1000.0_DP*(0.1880_DP*sarcoLengthAtActivation**4-0.8694_DP*sarcoLengthAtActivation**3+ &
                & 1.5084_DP*sarcoLengthAtActivation**2-1.1577_DP*sarcoLengthAtActivation+0.3345_DP)
              IF(elongation<0.02) THEN ! 0.02 -> offset value 
                force=0.0_DP
              ELSE 
                force=(elongation-0.02_DP)*stiffnessPEVK
              ENDIF
              titinBound=F0+force
              
              !square approximation
            ELSE IF(switchModel==2) THEN
              lengthInitTitin=sarcoLengthAtActivation-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
              !function to approximate the relation between the initial titin length and the initial passive force F0	
              F0=COEFF_MATRIX(1)*EXP(lengthInitTitin)+COEFF_MATRIX(2)*lengthInitTitin**3+COEFF_MATRIX(3)* &
                & lengthInitTitin**2+COEFF_MATRIX(4)*lengthInitTitin+COEFF_MATRIX(5)
              force=2.0_DP*elongation**2
              titinBound=F0+force
              
              !Rode's model -> Ekin's implementation
            ELSE IF(switchModel==3) THEN
              lengthInitTitin=sarcoLengthAtActivation-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
              !function to approximate the relation between the initial titin length and the initial passive force F0	
              F0=COEFF_MATRIX(1)*EXP(lengthInitTitin)+COEFF_MATRIX(2)*lengthInitTitin**3+COEFF_MATRIX(3)* &
                & lengthInitTitin**2+COEFF_MATRIX(4)*lengthInitTitin+COEFF_MATRIX(5)
              !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
              stiffnessPEVK=1000.0_DP*(0.1880_DP*sarcoLengthAtActivation**4-0.8694_DP*sarcoLengthAtActivation**3+ &
                & 1.5084_DP*sarcoLengthAtActivation**2-1.1577_DP*sarcoLengthAtActivation+0.3345_DP)
              
              IF(F0<=ZERO_TOLERANCE) THEN
                lengthDistIGF0=-35.63_DP+39.58889_DP*(F0+0.9_DP)
              ELSEIF(F0>=0.24_DP) THEN
                lengthDistIGF0=0.1411_DP+0.196576763485477_DP*(F0-0.2495_DP)
              ELSE                    
                indexPseudo=CEILING(F0/DX)
                indexI=indexRef+indexPseudo-1
                lengthDistIGF0=lengthsDistIG(indexI)-(lengthsDistIG(indexI+1)-lengthsDistIG(indexI))* &
                  & (forcesDistIG(indexI)-F0)/(forcesDistIG(indexI+1)-forcesDistIG(indexI))
              ENDIF
              
              elongationNew=-1.0_DP
              force=0.0_DP                  
              DO WHILE(elongationNew<elongation)
                force=force+FORCE_INCREMENT
                titinBound=force+F0
                elongationPEVK=force/stiffnessPEVK  
                IF (titinBound<=0.0_DP) THEN
                  lengthDistIG=-35.63_DP+39.58889_DP*(titinBound+0.9_DP)
                ELSE IF(titinBound>=0.24_DP) THEN
                  lengthDistIG=0.1411_DP+0.196576763485477_DP*(titinBound-0.2495_DP)
                ELSE                  
                  indexPseudo=CEILING(titinBound/DX)
                  indexI=indexRef+indexPseudo-1
                  lengthDistIG=lengthsDistIG(indexI)-(lengthsDistIG(indexI+1)-lengthsDistIG( &
                    & indexI))*(forcesDistIG(indexI)-titinBound)/(forcesDistIG(indexI+1)-forcesDistIG(indexI))
                END IF
                elongationDistIG=lengthDistIG-lengthDistIGF0                  
                elongationNew=elongationPEVK+elongationDistIG
              ENDDO
              
              !Rode's titin model -> solve with Newton's method
            ELSE IF(switchModel==4) THEN
              lengthInitTitin=sarcoLengthAtActivation-LENGTH_MYOSIN-LENGTH_MBAND-LENGTH_ZDISC
              !function to approximate the relation between the initial titin length and the initial passive force F0	
              F0=COEFF_MATRIX(1)*EXP(lengthInitTitin)+COEFF_MATRIX(2)*lengthInitTitin**3+COEFF_MATRIX(3)* &
                & lengthInitTitin**2+COEFF_MATRIX(4)*lengthInitTitin+COEFF_MATRIX(5)
              !function to approximate the relation between the initial sarcomere length and the stiffness of the PEVK region.
              stiffnessPEVK=1000.0_DP*(0.1880_DP*sarcoLengthAtActivation**4-0.8694_DP*sarcoLengthAtActivation**3+ &
                & 1.5084_DP*sarcoLengthAtActivation**2-1.1577_DP*sarcoLengthAtActivation+0.3345_DP)
              
              !calculate lengthDistIGF0 with linear inter- or extrapolation
              indexI=2
              slope=0.0_DP
              lengthDistIGF0=0.0_DP
              IF(F0>=forcesDistIG(DIM_DATA)) THEN
                indexI=DIM_DATA
              ELSE
                DO WHILE(F0>forcesDistIG(indexI))
                  indexI=indexI+1
                ENDDO
              ENDIF
              slope=(forcesDistIG(indexI)-forcesDistIG(indexI-1))/ &
                & (lengthsDistIG(indexI)-lengthsDistIG(indexI-1))
              lengthDistIGF0=lengthsDistIG(indexI-1)+slope*(F0-forcesDistIG(indexI-1))

              !initialize Newton-method to calculate the titin force
              indexI=2                
              stiffnessDist=1.0_DP   
              forceDistalIG=F0      
              force=0.0_DP                  ! delta P
              titinBound=0.0_DP            ! total titin force (= P_PEVK)
              deltaF=10.0_DP               ! start value for iteration
              diffQuot=1.0_DP              ! numerical derivative              
              elongationPEVK=elongation/2.0_DP                               
              lengthDistIG=lengthDistIGF0+elongation-elongationPEVK     
              
              DO WHILE(ABS(deltaF)>TOL) !Newton-method (solve forceDistalIG-force_0=0)
                
                IF(lengthDistIG>=lengthsDistIG(DIM_DATA)) THEN  !Extrapolation if lengthsDistIG(end)>lengthDistIG
                  indexI=DIM_DATA
                ELSE 
                  DO WHILE(lengthDistIG>lengthsDistIG(indexI))
                    indexI=indexI+1
                  ENDDO
                ENDIF
                stiffnessDist=(forcesDistIG(indexI)-forcesDistIG(indexI-1))/ &
                  & (lengthsDistIG(indexI)-lengthsDistIG(indexI-1))
                forceDistalIG=forcesDistIG(indexI-1)+stiffnessDist* &
                  & (lengthDistIG-lengthsDistIG(indexI-1))                
                
                force=stiffnessPEVK*elongationPEVK
                titinBound=force+F0
                
                deltaF=titinBound-forceDistalIG !new iteration for forceDistalIG-titinBound
                diffQuot=stiffnessPEVK+stiffnessDist                 
                elongationPEVK=elongationPEVK-deltaF/diffQuot                 
                lengthDistIG=lengthDistIGF0+elongation-elongationPEVK       
              ENDDO
            ENDIF ! switch model
            
            ! calculate x-fibre titin stress                       
            ! fitting function to calculate the filament-lattice parameter d10 in nanometer (from Elliot 1963 -mammalian muscle)
            d10=-13.39_DP*sarcoLength+58.37_DP
            ! d10=-13.39_DP*1.0_DP+58.37_DP
            ! calculate the distance between actin and myosin filament in micro meter (geometrical relationship)
            actinMyosinDistance=0.001_DP*(2.0_DP/3.0_DP*d10) 
            ! calculate x-fibre stress with tangent (unbound titin)
            titinXFUnbound=0.5_DP*titinUnbound*actinMyosinDistance/lengthTitin 
            ! calculate x-fibre stress with tangent (for bound titin filaments)
            titinXFBound=0.5_DP*titinBound*actinMyosinDistance/lengthDistIG  
            
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,2, &
              & titinUnbound,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,3, &
              & titinBound,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,4, &
              & titinXFUnbound,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalNode(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,5, &
              & titinXFBound,err,error,*999)
            
          ENDIF ! Check if elongation is positive or not
        ENDDO ! Over the nodes

        !now the ghost elements -- get the relevant info from the other computational nodes
        CALL FieldVariable_ParameterSetUpdateStart(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(monodomainIndependentUVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        
      ENDIF
    CASE DEFAULT
      CALL FlagError("Problem subtype not implemented for titin",err,error,*999)
    END SELECT

    EXITS("BioelectricFiniteElasticity_ComputeTitin")
    RETURN
999 ERRORSEXITS("BioelectricFiniteElasticity_ComputeTitin",err,error)
    RETURN 1

  END SUBROUTINE BioelectricFiniteElasticity_ComputeTitin

  !
  !================================================================================================================================
  !

END MODULE BioelectricFiniteElasticityRoutines

