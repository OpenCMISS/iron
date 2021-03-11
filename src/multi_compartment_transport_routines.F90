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
!> Contributor(s): Andrew Cookson, Chris Bradley
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

!>TThis module handles all routines pertaining to (advection-)diffusion coupled to (advection-)diffusion.


MODULE MultiCompartmentTransportRoutines

  USE AdvectionDiffusionEquationsRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
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
  USE FieldAccessRoutines
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

  PUBLIC MultiCompartmentTransport_EquationsSetSetup
  
  PUBLIC MultiCompartmentTransport_EquationsSetSolutionMethodSet

  PUBLIC MultiCompartmentTransport_ProblemSetup
  
  PUBLIC MultiCompartmentTransport_ProblemSpecificationSet
  
  PUBLIC MultiCompartmentTransport_FiniteElementCalculate

  PUBLIC MultiCompartmentTransport_PreSolve
  
  PUBLIC MultiCompartmentTransport_PostSolve
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a coupled diffusion & advection-diffusion equation type of a multi physics equations set class.
  SUBROUTINE MultiCompartmentTransport_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("MultiCompartmentTransport_EquationsSetSolutionMethodSet",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)
       
    EXITS("MultiCompartmentTransport_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("MultiCompartmentTransport_EquationsSetSolutionMethodSet",err,error)
    EXITS("MultiCompartmentTransport_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the multi-compartment coupled advection-diffusion & diffusion transport equation.
  SUBROUTINE MultiCompartmentTransport_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MultiCompartmentTransport_EquationsSetSetup",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)
             
    EXITS("MultiCompartmentTransport_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("MultiCompartmentTransport_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE MultiCompartmentTransport_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a multi-compartment coupled advection-diffusion & diffusion transport equation.
  SUBROUTINE MultiCompartmentTransport_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("MultiCompartmentTransport_FiniteElementCalculate",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)
      
    EXITS("MultiCompartmentTransport_FiniteElementCalculate")
    RETURN
999 ERRORS("MultiCompartmentTransport_FiniteElementCalculate",err,error)
    EXITS("MultiCompartmentTransport_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the problem subtype for a coupled diffusion & advection-diffusion equation type .
  SUBROUTINE MultiCompartmentTransport_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiCompartmentTransport_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a multi-compartment coupled transport equation type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_MULTI_PHYSICS_CLASS,PROBLEM_MULTI_COMPARTMENT_TRANSPORT_TYPE,problemSubtype]

    EXITS("MultiCompartmentTransport_ProblemSpecificationSet")
    RETURN
999 ERRORS("MultiCompartmentTransport_ProblemSpecificationSet",err,error)
    EXITS("MultiCompartmentTransport_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the coupled diffusion-diffusion equations problem.
  SUBROUTINE MultiCompartmentTransport_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solverDiffusion, solverAdvectionDiffusion
    TYPE(SolverEquationsType), POINTER :: solverEquationsDiffusion, solverEquationsAdvectionDiffusion
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiCompartmentTransport_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " does not equal a standard multi-component transport equation subtype."
      CALL FlagError(localError,err,error,*999)      
    END SELECT
 
    !--------------------------------------------------------------------
    !   monolithic coupled source diffusion-diffusion problem
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
          & " is invalid for a multi-compartment transport equation."
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
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a multi-compartment transport equation."
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
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !Set the solver to be a linear solver for the diffusion problem
        NULLIFY(solverDiffusion)
        CALL Solvers_SolverGet(solvers,1,solverDiffusion,err,error,*999)
        CALL Solver_TypeSet(solverDiffusion,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_DynamicOrderSet(solverDiffusion,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solverDiffusion,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solverDiffusion,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_LibraryTypeSet(solverDiffusion,SOLVER_CMISS_LIBRARY,err,error,*999)
        !
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a multi-compartment transport equation."
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
        !Get the diffusion solver and create the diffusion solver equations
        NULLIFY(solverDiffusion)
        CALL Solvers_SolverGet(solvers,1,solverDiffusion,err,error,*999)
        NULLIFY(solverEquationsDiffusion)
        CALL SolverEquations_CreateStart(solverDiffusion,solverEquationsDiffusion,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquationsDiffusion,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquationsDiffusion,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquationsDiffusion,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the creation of the diffusion solver equations
        NULLIFY(solverDiffusion)
        CALL Solvers_SolverGet(solvers,1,solverDiffusion,err,error,*999)
        NULLIFY(solverEquationsDiffusion)
        CALL Solver_SolverEquationsGet(solverDiffusion,solverEquationsDiffusion,err,error,*999)
        CALL SolverEquations_CreateFinish(solverEquationsDiffusion,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a multi-compartment transport equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a multi-compartment transport equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("MultiCompartmentTransport_ProblemSetup")
    RETURN
999 ERRORSEXITS("MultiCompartmentTransport_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_ProblemSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the multi-compartment coupled transport problem pre-solve.
  SUBROUTINE MultiCompartmentTransport_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations 
    TYPE(SolverMappingType), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiCompartmentTransport_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
      IF(ASSOCIATED(equationsAnalytic)) THEN
        CALL MultiCompartmentTransport_PreSolveUpdateAnalyticValues(solver,err,error,*999)
        !CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
        !IF(solverNumber==1) THEN
        !  !copy current value of concentration_one to another variable
        !  CALL ADVEC_DIFFUSION_EQUATION_PRE_SOLVE_STORE_CURRENT_SOLN(solver,err,error,*999)
        !  !Set source term to be updated value of concentration_two
        !  CALL ADVECTION_DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(solver,err,error,*999)
        !ELSE IF(solverNumber==2) THEN
        !  !compute value of constant source term - evaluated from lamdba*(0.5*(c_1^{t+1}+c_1^{t}) - c_2^{t})
        !  CALL DIFFUSION_EQUATION_PRE_SOLVE_GET_SOURCE_VALUE(solver,err,error,*999)
        !ENDIF
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a multi-compartment transport type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("MultiCompartmentTransport_PreSolve")
    RETURN
999 ERRORSEXITS("MultiCompartmentTransport_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_PreSolve
      
  !   
  !================================================================================================================================
  !
  !updates the boundary conditions and source term to the required analytic values
  SUBROUTINE MultiCompartmentTransport_PreSolveUpdateAnalyticValues(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionCheckVariable,componentIdx,derivativeIdx,dimensionIdx, &
      & equationsSetIdx,globalDerivativeIndex,localDOFIdx,nodeIdx,numberOfComponents,numberOfDimensions,numberOfEquationsSets, &
      & numberOfNodes,numberOfNodeDerivatives,pSpecification(3),variableIndex,variableType
    REAL(DP) :: analyticValue,currentTime,normal(3),tangents(3,3),timeIncrement,X(3),valueSource
!     REAL(DP) :: k_xx, k_yy, k_zz
    REAL(DP) :: A1,A2,A3,A4,D1,D2,D3,D4,lambda12,lambda13,lambda23
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet 
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField,sourceField
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,geometricVariable,materialsVariable,sourceVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
!    REAL(DP), POINTER :: boundaryValues(:)

    ENTERS("MultiCompartmentTransport_PreSolveUpdateAnalyticValues",err,error,*999)

    A1=0.4_DP
    A2=0.3_DP
    A3=0.2_DP
    A4=0.1_DP
    D1=1.0_DP
    D2=1.0_DP
    D3=1.0_DP
    D4=1.0_DP
    lambda12=0.1_DP
    lambda13=0.1_DP
    lambda23=0.1_DP

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_Specification(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !loop over all the equation sets and set the appropriate field variable type BCs and
      !the source field associated with each equation set
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
          NULLIFY(geometricParameters)
          CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(analyticField)
          CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
          NULLIFY(analyticVariable)
          NULLIFY(analyticParameters)
          IF(ASSOCIATED(analyticField)) THEN
            CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
            CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
          ENDIF
          NULLIFY(materialsField)
          CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
          NULLIFY(materialsVariable)
          NULLIFY(materialsParameters)
          IF(ASSOCIATED(materialsField)) THEN
            CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
            CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
          ENDIF
          NULLIFY(sourceField)
          CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
          NULLIFY(sourceVariable)
          IF(ASSOCIATED(sourceField)) CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
          !CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          !DO variableIdx=1,numberOfVariables
          variableIndex=2*equationsSetIdx-1
          NULLIFY(dependentVariable)
          CALL Field_VariableIndexGet(dependentField,variableIndex,dependentVariable,variableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
              & err,error,*999)
            NULLIFY(domain)
            CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodes(domainTopology,domainNodes,err,error,*999)
            NULLIFY(boundaryConditions)
            CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
            NULLIFY(boundaryConditionsVariable)
            CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
              DO dimensionIdx=1,numberOfDimensions
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDofGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
                x(dimensionIdx)=geometricParameters(localDOFIdx)
              ENDDO !dimensionIdx
              CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
              !Loop over the derivatives
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
                CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,currentTime, &
                  & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,analyticValue, &
                  & err,error,*999)
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                  & err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                  & analyticValue,err,error,*999)
                CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOFIdx, &
                  & boundaryConditionCheckVariable,err,error,*999)
                IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) THEN
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                    & analyticValue,err,error,*999)
                ENDIF                
                !IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                !  IF(boundaryNode) THEN
                !    !If we are a boundary node then set the analytic value on the boundary
                !    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,variable_type,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
                !      & analyticValue,err,error,*999)
                !  ENDIF
                !ENDIF
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Sources
          IF(ASSOCIATED(sourceField)) THEN
            CALL FieldVariable_NumberOfComponentsGet(sourceVariable,numberOfComponents,err,error,*999)
            DO componentIdx=1,numberOfComponents
              CALL FieldVariable_ComponentInterpolationCheck(sourceVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
                & err,error,*999)
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(sourceVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodes(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
                DO dimensionIdx=1,numberOfDimensions
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDofGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
                  x(dimensionIdx)=geometricParameters(localDOFIdx)
                ENDDO !dimensionIdx
                !Loop over the derivatives
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  SELECT CASE(analyticFunctionType)
                  CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
                    SELECT CASE(equationsSetIdx)
                    CASE(1)
                      valueSource=EXP(-1.0_DP*currentTime)*(-1.0_DP*A1*(X(1)*X(1)+X(2)*X(2))-4.0_DP*D1*A1+lambda12*(A1-A2)* &
                        & (X(1)*X(1)+X(2)*X(2)))
                    CASE(2)
                      valueSource=EXP(-1.0_DP*currentTime)*(-1.0_DP*A2*(X(1)*X(1)+X(2)*X(2))-4.0_DP*D2*A2+lambda12*(A2-A1)* &
                        & (X(1)*X(1)+X(2)*X(2)))
                    END SELECT
                  CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
                    SELECT CASE(equationsSetIdx)
                    CASE(1)
                      valueSource=EXP(-1.0_DP*currentTime)*(-1.0_DP*A1*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))- &
                        & 6.0_DP*D1*A1+lambda13*(A1-A3)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+ &
                        & lambda12*(A1-A2)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                    CASE(2)
                      valueSource=EXP(-1.0_DP*currentTime)*(-1.0_DP*A2*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))- &
                        & 6.0_DP*D2*A2+lambda12*(A2-A1)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+ &
                        & lambda23*(A2-A3)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                    CASE(3)
                      valueSource=EXP(-1.0_DP*currentTime)*(-1.0_DP*A3*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))- &
                        & 6.0_DP*D3*A3+lambda13*(A3-A1)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))+ &
                        & lambda23*(A3-A2)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)))
                    END SELECT
                  CASE DEFAULT
                    localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))//&
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(sourceVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                    & err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(sourceVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,valueSource, &
                    & err,error,*999)
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
            CALL FieldVariable_ParameterSetUpdateStart(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          ENDIF
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(ASSOCIATED(sourceField)) &
            & CALL FieldVariable_ParameterSetUpdateFinish(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)          
        !ENDDO !variable_idx
          CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
          CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateStart(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a multi-physics coupled diffusion equation type of a multi-physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("MultiCompartmentTransport_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORS("MultiCompartmentTransport_PreSolveUpdateAnalyticValues",err,error)
    EXITS("MultiCompartmentTransport_PreSolveUpdateAnalyticValues")
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_PreSolveUpdateAnalyticValues
  !   
  !================================================================================================================================
  !
  !>Sets up the multi-compartment transport problem post solve.
  SUBROUTINE MultiCompartmentTransport_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3),solverNumber
    !REAL(DP), POINTER :: outputData1(:),outputData2(:),outputData3(:),outputData4(:),outputData5(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("MultiCompartmentTransport_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      IF(solverNumber==1) THEN
        !CALL AdvectionDiffusion_PostSolve(solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)

        !NULLIFY(outputData1)
        !CALL Field_ParameterSetDataGet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, & 
        !  & FIELD_VALUES_SET_TYPE,outputData1,err,error,*999)
        !NULLIFY(outputData2)
        !CALL Field_ParameterSetDataGet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, & 
        !  & FIELD_VALUES_SET_TYPE,outputData2,err,error,*999)
        !NULLIFY(outputData3)
        !CALL Field_ParameterSetDataGet(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, & 
        !  & FIELD_VALUES_SET_TYPE,outputData3,err,error,*999)
        !NULLIFY(outputData4)
        !CALL Field_ParameterSetDataGet(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, & 
        !  & FIELD_VALUES_SET_TYPE,outputData4,err,error,*999)
        !NULLIFY(outputData5)
        !CALL Field_ParameterSetDataGet(equationsSet%dependent%dependentField,FIELD_U3_VARIABLE_TYPE, & 
        !  & FIELD_VALUES_SET_TYPE,outputData5,err,error,*999)
        
      ELSE IF(solverNumber==2) THEN
        !CALL Diffusion_PostSolve(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(controlLoop%problem%SPECIFICATION(3),"*",err,error))// &
        & " is not valid for a multi-compartment type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("MultiCompartmentTransport_PostSolve")
    RETURN
999 ERRORSEXITS("MultiCompartmentTransport_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_PostSolve

  !
  !================================================================================================================================
  !

  !>Sets up the diffuion-diffusion problem post solve output data.
  SUBROUTINE MultiCompartmentTransport_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("MultiCompartmentTransport_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_MULTI_COMPARTMENT_TRANSPORT_SUBTYPE)
      !CALL AdvectionDiffusion_PostSolveOutputData(solver,err,error,*999)
      !CALL Diffusion_PostSolveOutputData(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a multi-compartment transport type of a multi physics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("MultiCompartmentTransport_PostSolveOutputData")
    RETURN
999 ERRORS("MultiCompartmentTransport_PostSolveOutputData",err,error)
    EXITS("MultiCompartmentTransport_PostSolveOutputData")
    RETURN 1
    
  END SUBROUTINE MultiCompartmentTransport_PostSolveOutputData
      
  !   
  !================================================================================================================================
  !

END MODULE MultiCompartmentTransportRoutines
