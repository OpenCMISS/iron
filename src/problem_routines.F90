!> \file
!> \author Chris Bradley
!> \brief This module handles all problem routines.
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

!> This module handles all problem routines.
MODULE ProblemRoutines

  USE BaseRoutines
  USE BIOELECTRIC_ROUTINES
  USE ClassicalFieldRoutines
  USE ComputationAccessRoutines
  USE ContextAccessRoutines
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE ElasticityRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsSetConstants
  USE EquationsSetRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FIELD_IO_ROUTINES
  USE FINITE_ELASTICITY_ROUTINES
  USE FittingRoutines
  USE FluidMechanicsRoutines
  USE InputOutput
  USE INTERFACE_CONDITIONS_CONSTANTS
  USE INTERFACE_CONDITIONS_ROUTINES
  USE InterfaceRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MULTI_PHYSICS_ROUTINES
  USE PROBLEM_CONSTANTS
  USE ProblemAccessRoutines
  USE REACTION_DIFFUSION_EQUATION_ROUTINES
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SOLVER_MATRICES_ROUTINES
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Problems_Initialise,Problems_Finalise
  
  PUBLIC Problem_CellMLEquationsCreateStart,Problem_CellMLEquationsCreateFinish
  
  PUBLIC Problem_CreateStart,Problem_CreateFinish,Problem_Destroy
  
  PUBLIC Problem_SpecificationGet,Problem_SpecificationSizeGet
  
  PUBLIC Problem_ControlLoopCreateStart,Problem_ControlLoopCreateFinish
  
  PUBLIC Problem_ControlLoopDestroy
  
  PUBLIC Problem_SolverDAECellMLRHSEvaluate
  
  PUBLIC Problem_SolverEquationsBoundaryConditionsAnalytic

  PUBLIC Problem_SolverEquationsCreateStart,Problem_SolverEquationsCreateFinish
  
  PUBLIC Problem_SolverEquationsDestroy
  
  PUBLIC Problem_SolverJacobianEvaluate,Problem_SolverResidualEvaluate
  
  PUBLIC Problem_SolverNonlinearMonitor

  PUBLIC Problem_SolverOptimiserMonitor
  
  PUBLIC Problem_Solve
  
  PUBLIC Problem_SolversCreateStart,Problem_SolversCreateFinish
  
  PUBLIC Problem_SolversDestroy

  PUBLIC Problem_WorkGroupSet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the creation of the CellML equations for the problem solver. \see OpenCMISS::Iron::cmfe_Problem_SolverCellMLEquationsCreateFinish
  SUBROUTINE Problem_CellMLEquationsCreateFinish(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finish the CellML equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: PROBLEM_SETUP_INFO

    ENTERS("Problem_CellMLEquationsCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(PROBLEM_SETUP_INFO,err,error,*999)
    PROBLEM_SETUP_INFO%setupType=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
    PROBLEM_SETUP_INFO%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish problem specific startup
    CALL Problem_Setup(problem,PROBLEM_SETUP_INFO,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(PROBLEM_SETUP_INFO,err,error,*999)
      
    EXITS("Problem_CellMLEquationsCreateFinish")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsCreateFinish
  
  !
  !================================================================================================================================
  !

  !>Start the creation of CellML equations for a problem solver. \see OpenCMISS::Iron::cmfe_Problem_SolverCellMLEquationsCreateStart
  SUBROUTINE Problem_CellMLEquationsCreateStart(problem,err,error,*)

    !Argument variablesg
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to start the creation of the CellML equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: PROBLEM_SETUP_INFO

    ENTERS("Problem_CellMLEquationsCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(PROBLEM_SETUP_INFO,err,error,*999)
    PROBLEM_SETUP_INFO%setupType=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
    PROBLEM_SETUP_INFO%actionType=PROBLEM_SETUP_START_ACTION
    !Start the problem specific control setup
    CALL Problem_Setup(problem,PROBLEM_SETUP_INFO,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(PROBLEM_SETUP_INFO,err,error,*999)
       
    EXITS("Problem_CellMLEquationsCreateStart")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsCreateStart

  !
  !================================================================================================================================
  !

  !>Solves CellML equations for a problem.
  SUBROUTINE Problem_CellMLEquationsSolve(cellMLEquations,err,error,*)

   !Argument variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: cellMLEquations !<A pointer to the CellML equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_CellMLEquationsSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
    IF(.NOT.cellMLEquations%CELLML_EQUATIONS_FINISHED) CALL FlagError("CellML equations have not been finished.",err,error,*999)

    NULLIFY(solver)
    CALL CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*999)
    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"CellML equations solve: ",solver%label,err,error,*999)
    ENDIF

    SELECT CASE(cellMLEquations%timeDependence)
    CASE(CELLML_EQUATIONS_STATIC)
      !Do nothing
    CASE(CELLML_EQUATIONS_QUASISTATIC,CELLML_EQUATIONS_DYNAMIC)
      NULLIFY(controlLoop)
      CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      CALL CellMLEquations_TimeSet(cellMLEquations,currentTime,err,error,*999)
    CASE DEFAULT
      localError="The CellML equations time dependence type of "// &
        & TRIM(NumberToVString(cellMLEquations%timeDependence,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    CALL Solver_Solve(solver,err,error,*999)
      
    EXITS("Problem_CellMLEquationsSolve")
    RETURN
999 ERRORSEXITS("Problem_CellMLEquationsSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CellMLEquationsSolve

  !
  !================================================================================================================================
  !

  !>Solves CellML equations for a problem.
  SUBROUTINE Problem_SolverDAECellMLRHSEvaluate(cellML,time,dofIdx,stateData,rateData,err,error,*)

   !Argument variables
    TYPE(CELLML_TYPE), POINTER :: cellML !<A pointer to the CellML to evaluate
    REAL(DP), INTENT(IN) :: time !<The time to evaluate the CellML model at
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The index of the DOF to evaluate
    REAL(DP), POINTER :: stateData(:) !<The states data to evaluate the model at
    REAL(DP), POINTER :: rateData(:) !<On exit, the evaluated rates
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofOrderType,intermediateDataOffset,maxNumberOfIntermediates,maxNumberOfParameters,maxNumberOfStates, &
      modelIdx,parameterDataOffset
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP), POINTER :: intermediateData(:),parameterData(:)
    TYPE(CELLML_MODEL_TYPE), POINTER :: model
    TYPE(FieldType), POINTER :: intermediateField,modelsField,parametersField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    
    ENTERS("Problem_SolverDAECellMLRHSEvaluate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(cellML)) CALL FlagError("CellML is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(cellML%MODELS_FIELD)) CALL FlagError("CellML models field is not associated.",err,error,*999)
    modelsField=>cellML%MODELS_FIELD%MODELS_FIELD
    IF(.NOT.ASSOCIATED(modelsField)) CALL FlagError("Models field not associated.",err,error,*999)
   
    maxNumberOfStates=cellML%MAXIMUM_NUMBER_OF_STATE
    maxNumberOfIntermediates=cellML%MAXIMUM_NUMBER_OF_INTERMEDIATE
    maxNumberOfParameters=cellML%maximumNumberOfParameters
    !Make sure CellML fields have been updated to the current value of any mapped fields
    NULLIFY(modelsVariable)
    CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
    CALL Field_DOFOrderTypeGet(modelsField,FIELD_U_VARIABLE_TYPE,dofOrderType,err,error,*999)
    CALL Field_ParameterSetDataGet(modelsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
    modelIdx=modelsData(dofIdx)
    model=>cellML%models(modelIdx)%ptr
    IF(.NOT.ASSOCIATED(model)) CALL FlagError("Model is not associated.",err,error,*999)
    IF(dofOrderType==FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
      parameterDataOffset=modelsVariable%totalNumberOfDofs
      intermediateDataOffset=modelsVariable%totalNumberOfDofs
    ELSE
      parameterDataOffset=maxNumberOfParameters
      intermediateDataOffset=maxNumberOfIntermediates
    ENDIF
    NULLIFY(parameterData)
    !Get the parameters information if this environment has any.
    IF(ASSOCIATED(cellML%PARAMETERS_FIELD)) THEN
      parametersField=>cellML%PARAMETERS_FIELD%PARAMETERS_FIELD
      IF(ASSOCIATED(parametersField)) THEN
        CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parameterData, &
          & err,error,*999)
      ENDIF
    ENDIF
    !Get the intermediate information if this environment has any.
    NULLIFY(intermediateData)
    IF(ASSOCIATED(cellML%INTERMEDIATE_FIELD)) THEN
      intermediateField=>cellml%INTERMEDIATE_FIELD%INTERMEDIATE_FIELD
      IF(ASSOCIATED(intermediateField)) THEN
        CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData, &
          & err,error,*999)
      ENDIF
    ENDIF!associated intermediate
    
    !Evaluate the CellML RHS
    CALL Solver_DAECellMLRHSEvaluate(model,time,1,1,stateData,dofIdx,parameterDataOffset,parameterData,dofIdx, &
      intermediateDataOffset,intermediateData,1,1,rateData,err,error,*999)
    
    EXITS("Problem_SolverDAECellMLRHSEvaluate")
    RETURN
999 ERRORSEXITS("Problem_SolverDAECellMLRHSEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverDAECellMLRHSEvaluate

  !
  !================================================================================================================================
  !

  !>Solves a problem control loop.
  RECURSIVE SUBROUTINE Problem_ControlLoopSolve(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iterationIdx,loopIdx,solverIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop
    TYPE(ControlLoopSimpleType), POINTER :: simpleLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_ControlLoopSolve",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)

    !Solve this control loop
    IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Control loop: ",controlLoop%label,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Control loop level = ",controlLoop%controlLoopLevel,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Sub loop index     = ",controlLoop%subLoopIndex,err,error,*999)
    ENDIF
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Control loop: ",controlLoop%label,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Control loop level = ",controlLoop%controlLoopLevel,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sub loop index     = ",controlLoop%subLoopIndex,err,error,*999)
    ENDIF
    SELECT CASE(controlLoop%loopType)
    CASE(CONTROL_SIMPLE_TYPE)
      NULLIFY(simpleLoop)
      CALL ControlLoop_SimpleLoopGet(controlLoop,simpleLoop,err,error,*999)
      IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Simple control loop: ",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Simple control loop: ",err,error,*999)
      ENDIF
      CALL Problem_ControlLoopPreLoop(controlLoop,err,error,*999)
      IF(controlLoop%numberOfSubLoops==0) THEN
        !If there are no sub loops then solve.
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
          
          CALL Problem_SolverSolve(solver,err,error,*999)
          
        ENDDO !solverIdx
      ELSE
        !If there are sub loops the recursively solve those control loops
        DO loopIdx=1,controlLoop%numberOfSubLoops
          NULLIFY(controlLoop2)
          CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
          CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
        ENDDO !loopIdx
      ENDIF
      CALL Problem_ControlLoopPostLoop(controlLoop,err,error,*999)
    CASE(CONTROL_FIXED_LOOP_TYPE)
      NULLIFY(fixedLoop)
      CALL ControlLoop_FixedLoopGet(controlLoop,fixedLoop,err,error,*999)
      DO iterationIdx=fixedLoop%startIteration,fixedLoop%stopIteration,fixedLoop%iterationIncrement
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Fixed control loop iteration: ",iterationIdx,err,error,*999)
        ENDIF
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Fixed control loop iteration: ",iterationIdx,err,error,*999)
        ENDIF
        fixedLoop%iterationNumber=iterationIdx
        CALL Problem_ControlLoopPreLoop(controlLoop,err,error,*999)
        IF(controlLoop%numberOfSubLoops==0) THEN
          !If there are no sub loops then solve
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,controlLoop%numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        CALL Problem_ControlLoopPostLoop(controlLoop,err,error,*999)
      ENDDO !iterationIdx
    CASE(CONTROL_TIME_LOOP_TYPE)
      NULLIFY(timeLoop)
      CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
      !Set the current time to be the start time. Solvers should use the first time step to do any initialisation.
      timeLoop%currentTime=timeLoop%startTime
      
      !Precompute the number of iterations from total time span and time increment if it was not specified explicitely 
      IF(timeLoop%numberOfIterations==0) THEN
        timeLoop%numberOfIterations=CEILING((timeLoop%stopTime-timeLoop%startTime)/timeLoop%timeIncrement)
        !If number of iterations was specified but does not match TIME_INCREMENT, e.g. TIME_INCREMENT is still at the default value, compute correct TIME_INCREMENT
      ELSE IF(CEILING((timeLoop%stopTime-timeLoop%startTime)/timeLoop%timeIncrement) /= timeLoop%numberOfIterations) THEN
        timeLoop%timeIncrement = (timeLoop%stopTime-timeLoop%startTime)/timeLoop%numberOfIterations
      ENDIF
      
      timeLoop%iterationNumber=0
      
      DO WHILE(timeLoop%iterationNumber<timeLoop%numberOfIterations)
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Time control loop iteration: ",timeLoop%iterationNumber, &
            & err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Total number of iterations: ",timeLoop%numberOfIterations, &
            & err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Current time   = ",timeLoop%currentTime,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Stop time      = ",timeLoop%stopTime,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Time increment = ",timeLoop%timeIncrement,err,error,*999)
        ENDIF
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Time control loop iteration: ",timeLoop%iterationNumber, &
            & err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Total number of iterations: ",timeLoop%numberOfIterations, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Current time   = ",timeLoop%currentTime,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Stop time      = ",timeLoop%stopTime,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Time increment = ",timeLoop%timeIncrement,err,error,*999)
        ENDIF
        !Perform any pre-loop actions.
        CALL Problem_ControlLoopPreLoop(controlLoop,err,error,*999)
        IF(controlLoop%numberOfSubLoops==0) THEN
          !If there are no sub loops then solve.
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,controlLoop%numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        !Perform any post loop actions.
        CALL Problem_ControlLoopPostLoop(controlLoop,err,error,*999)
        !Increment loop counter and time
        timeLoop%iterationNumber=timeLoop%iterationNumber+1
        timeLoop%globalIterationNumber=timeLoop%globalIterationNumber+1
        timeLoop%currentTime=timeLoop%currentTime+timeLoop%timeIncrement
      ENDDO !time loop
    CASE(CONTROL_WHILE_LOOP_TYPE)
      NULLIFY(whileLoop)
      CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,ERR,ERROR,*999)
      whileLoop%iterationNumber=0
      whileLoop%continueLoop=.TRUE.
      DO WHILE(whileLoop%continueLoop.AND.whileLoop%iterationNumber<whileLoop%maximumNumberOfIterations)
        whileLoop%iterationNumber=whileLoop%iterationNumber+1
        IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"While control loop iteration: ",whileLoop%iterationNumber,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of iterations = ", &
            & whileLoop%maximumNumberOfIterations,err,error,*999)
        ENDIF
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"While control loop iteration: ",whileLoop%iterationNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum number of iterations = ", &
            & whileLoop%maximumNumberOfIterations,err,error,*999)
        ENDIF
        CALL Problem_ControlLoopPreLoop(controlLoop,err,error,*999)
        IF(controlLoop%numberOfSubLoops==0) THEN
          !If there are no sub loops then solve
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            IF(ASSOCIATED(solver%SOLVER_EQUATIONS)) &
              & CALL Problem_SolverLoadIncrementApply(solver%SOLVER_EQUATIONS,1,1,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,controlLoop%numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        CALL Problem_ControlLoopPostLoop(controlLoop,err,error,*999)
      ENDDO !while loop
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
      NULLIFY(loadIncrementLoop)
      CALL ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*999)
      loadIncrementLoop%iterationNumber=0
      IF (loadIncrementLoop%maximumNumberOfIterations<1) THEN
        !Automatic stepping
        CALL FlagError("Automatic load incrementing is not implemented yet.",err,error,*999)
      ELSE
        !Fixed number of steps
        DO WHILE(loadIncrementLoop%iterationNumber<loadIncrementLoop%maximumNumberOfIterations)
          loadIncrementLoop%iterationNumber=loadIncrementLoop%iterationNumber+1
          IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
            CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Load increment control loop iteration: ", &
              & loadIncrementLoop%iterationNumber,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of iterations = ", &
              & loadIncrementLoop%maximumNumberOfIterations,err,error,*999)
          ENDIF
          IF(diagnostics1) THEN
            CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Load increment control loop iteration: ", &
              & loadIncrementLoop%iterationNumber,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Maximum number of iterations = ", &
              & loadIncrementLoop%maximumNumberOfIterations,err,error,*999)
          ENDIF
          CALL Problem_ControlLoopPreLoop(controlLoop,err,error,*999)
          IF(controlLoop%numberOfSubLoops==0) THEN
            !If there are no sub loops then solve
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
            DO solverIdx=1,solvers%NUMBER_OF_SOLVERS
              NULLIFY(solver)
              CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
              !Apply incremented boundary conditions here => 
              IF(ASSOCIATED(solver%SOLVER_EQUATIONS)) &
                & CALL Problem_SolverLoadIncrementApply(solver%SOLVER_EQUATIONS,loadIncrementLoop%iterationNumber, &
                & loadIncrementLoop%maximumNumberOfIterations,err,error,*999)
              
              CALL Problem_SolverSolve(solver,err,error,*999)
                
            ENDDO !solverIdx
          ELSE
            !If there are sub loops the recursively solve those control loops
            DO loopIdx=1,controlLoop%numberOfSubLoops
              NULLIFY(controlLoop2)
              CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
              CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
            ENDDO !loopIdx
          ENDIF
          CALL Problem_ControlLoopPostLoop(controlLoop,err,error,*999)
        ENDDO !while loop
      ENDIF
    CASE DEFAULT
      localError="The control loop loop type of "//TRIM(NumberToVString(controlLoop%loopType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Problem_ControlLoopSolve")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopSolve

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a problem. \see OpenCMISS::Iron::cmfe_Problem_CreateFinish
  SUBROUTINE Problem_CreateFinish(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finish creating.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemIdx
    TYPE(ProblemSetupType) :: problemSetupInfo
    TYPE(ProblemType), POINTER :: problem2
    TYPE(ProblemsType), POINTER :: problems

    ENTERS("Problem_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_INITIAL_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish the problem specific setup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
    !Finish the problem creation
    problem%problemFinished=.TRUE.
    
    IF(diagnostics1) THEN
      NULLIFY(problems)
      CALL Problem_ProblemsGet(problem,problems,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Number of problems = ",problems%numberOfProblems,err,error,*999)
      DO problemIdx=1,problems%numberOfProblems
        NULLIFY(problem2)
        CALL Problems_ProblemGet(problems,problemIdx,problem2,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Problem number  = ",problemIdx,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  User number     = ",problem2%userNumber,err,error,*999)
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Global number   = ",problem2%globalNumber,err,error,*999)
        CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,SIZE(problem2%specification,1),8,8, &
          & problem2%specification,'("  Problem specification = ",8(X,I3))','(16X,8(X,I3))', &
          & err,error,*999)
      ENDDO !problemIdx    
    ENDIF
    
    EXITS("Problem_CreateFinish")
    RETURN
999 ERRORSEXITS("Problem_CreateFinish",err,error)    
    RETURN 1
   
  END SUBROUTINE Problem_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating a problem defined by userNumber. \see OpenCMISS::Iron::cmfe_Problem_CreateStart
  SUBROUTINE Problem_CreateStart(userNumber,problems,problemSpecification,problem,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the problem to create
    TYPE(ProblemsType), POINTER :: problems !<The problems to create the problem for. 
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification array, eg. [problem_class, problem_type, problem_subtype]
    TYPE(ProblemType), POINTER :: problem !<On return, a pointer to the created problem. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemIdx
    TYPE(ComputationEnvironmentType), POINTER :: computationEnvironment
    TYPE(ContextType), POINTER :: context
    TYPE(ProblemType), POINTER :: newProblem
    TYPE(ProblemPtrType), POINTER :: newProblems(:)
    TYPE(ProblemSetupType) :: problemSetupInfo
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: worldWorkGroup
 
    NULLIFY(newProblem)
    NULLIFY(newProblems)

    ENTERS("Problem_CreateStart",err,error,*999)

    IF(ASSOCIATED(problem)) CALL FlagError("Problem is already associated.",err,error,*999)   
    NULLIFY(problem)
    CALL Problem_UserNumberFind(problems,userNumber,problem,err,error,*999)
    IF(ASSOCIATED(problem)) THEN
      localError="Problem number "//TRIM(NumberToVString(userNumber,"*",err,error))//" has already been created."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Allocate the new problem
    ALLOCATE(newProblem,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new problem.",err,error,*999)
    !Initalise problem
    CALL Problem_Initialise(newProblem,err,error,*999)
    !Set default problem values
    newProblem%userNumber=userNumber
    newProblem%globalNumber=problems%numberOfProblems+1
    newProblem%problems=>problems
    NULLIFY(context)
    CALL Problems_ContextGet(problems,context,err,error,*999)
    NULLIFY(computationEnvironment)
    CALL Context_ComputationEnvironmentGet(context,computationEnvironment,err,error,*999)
    NULLIFY(worldWorkGroup)
    CALL ComputationEnvironment_WorldWorkGroupGet(computationEnvironment,worldWorkGroup,err,error,*999)
    newProblem%workGroup=>worldWorkGroup
    !Set problem specification
    CALL Problem_SpecificationSet(newProblem,problemSpecification,err,error,*999)
    !For compatibility with old code, set class, type and subtype
    newProblem%problemFinished=.FALSE.
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_INITIAL_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_START_ACTION
    !Start problem specific setup
    CALL Problem_Setup(newProblem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
    !Add new problem into list of problems
    ALLOCATE(newProblems(problems%numberOfProblems+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new problems.",err,error,*999)
    DO problemIdx=1,problems%numberOfProblems
      newProblems(problemIdx)%ptr=>problems%problems(problemIdx)%ptr
    ENDDO !problemIdx
    newProblems(problems%numberOfProblems+1)%ptr=>newProblem
    IF(ASSOCIATED(problems%problems)) DEALLOCATE(problems%problems)
    problems%problems=>newProblems
    problems%numberOfProblems=problems%numberOfProblems+1
    problem=>newProblem
    
    EXITS("Problem_CreateStart")
    RETURN
999 ERRORSEXITS("Problem_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Destroys a problem. \see OpenCMISS::Iron::cmfe_Problem_Destroy
  SUBROUTINE Problem_Destroy(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to destroy 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemIdx,problemPosition
    TYPE(ProblemPtrType), POINTER :: newProblems(:)
    TYPE(ProblemsType), POINTER :: problems

    NULLIFY(newProblems)

    ENTERS("Problem_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*998)
    
    NULLIFY(problems)
    CALL Problem_ProblemsGet(problem,problems,err,error,*999)
    IF(.NOT.ASSOCIATED(problems%problems)) CALL FlagError("Problem problems is not associated.",err,error,*999)
        
    problemPosition=problem%globalNumber
      
    !Destroy all the problem components
    CALL Problem_Finalise(problem,err,error,*999)
    
    !Remove the problem from the list of problems
    IF(problems%numberOfProblems>1) THEN
      ALLOCATE(newProblems(problems%numberOfProblems-1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new problems.",err,error,*999)
      DO problemIdx=1,problems%numberOfProblems
        IF(problemIdx<problemPosition) THEN
          newProblems(problemIdx)%ptr=>problems%problems(problemIdx)%ptr
        ELSE IF(problemIdx>problemPosition) THEN
          problems%problems(problemIdx)%ptr%globalNumber=problems%problems(problemIdx)%ptr%globalNumber-1
          newProblems(problemIdx-1)%ptr=>problems%problems(problemIdx)%ptr
        ENDIF
      ENDDO !problemIdx
      DEALLOCATE(problems%problems)
      problems%problems=>newProblems
      problems%numberOfProblems=problems%numberOfProblems-1
    ELSE
      DEALLOCATE(problems%problems)
      problems%numberOfProblems=0
    ENDIF
        
    EXITS("Problem_Destroy")
    RETURN
999 IF(ASSOCIATED(newProblems)) DEALLOCATE(newProblems)
998 ERRORSEXITS("Problem_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Destroy
  
  !
  !================================================================================================================================
  !

  !>Finalise the problem setup and deallocate all memory.
  SUBROUTINE Problem_SetupFinalise(problemSetupInfo,err,error,*)

    !Argument variables
    TYPE(ProblemSetupType), INTENT(OUT) :: problemSetupInfo !<The problem setup to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_SetupFinalise",err,error,*999)

    problemSetupInfo%setupType=0
    problemSetupInfo%actionType=0
       
    EXITS("Problem_SetupFinalise")
    RETURN
999 ERRORSEXITS("Problem_SetupFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SetupFinalise

 !
  !================================================================================================================================
  !

  !>Initialise the problem setup.
  SUBROUTINE Problem_SetupInitialise(problemSetupInfo,err,error,*)

    !Argument variables
    TYPE(ProblemSetupType), INTENT(OUT) :: problemSetupInfo !<The problem setup to intialise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_SetupInitialise",err,error,*999)

    problemSetupInfo%setupType=0
    problemSetupInfo%actionType=0
        
    EXITS("Problem_SetupInitialise")
    RETURN
999 ERRORSEXITS("Problem_SetupInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SetupInitialise

  !
  !================================================================================================================================
  !

  !>Finalise the problem and deallocate all memory.
  SUBROUTINE Problem_Finalise(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_Finalise",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(ASSOCIATED(problem%controlLoop)) CALL ControlLoop_Destroy(problem%controlLoop,err,error,*999)
      IF(ALLOCATED(problem%specification)) DEALLOCATE(problem%specification)
      DEALLOCATE(problem)
    ENDIF
       
    EXITS("Problem_Finalise")
    RETURN
999 ERRORSEXITS("Problem_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Finalise

  !
  !================================================================================================================================
  !

  !>Initialises a problem.
  SUBROUTINE Problem_Initialise(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<The pointer to the problem
    INTEGER(INTG), INTENT(OUT) :: err !<The error code !<The errror code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string !<The error string
    !Local Variables
 
    ENTERS("Problem_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    problem%userNumber=0
    problem%globalNumber=0
    problem%problemFinished=.FALSE.
    NULLIFY(problem%problems)
    NULLIFY(problem%workGroup)
    NULLIFY(problem%controlLoop)
       
    EXITS("Problem_Initialise")
    RETURN
999 ERRORSEXITS("Problem_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Initialise

  !
  !================================================================================================================================
  !

  !>Finish the creation of the control for the problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopCreateFinish
  SUBROUTINE Problem_ControlLoopCreateFinish(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finish the control for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_ControlLoopCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_CONTROL_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish problem specific startup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
    !Finish problem control creation
    controlLoop%controlLoopFinished=.TRUE.
     
    EXITS("Problem_ControlLoopCreateFinish")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopCreateFinish
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a control loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopCreateStart
  !>The default values of the PROBLEM CONTROL LOOP attributes are:
  !>- LOOP_TYPE: PROBLEM_CONTROL_SIMPLE_TYPE
  !>- CONTROL_LOOP_LEVEL: 1
  !>- NUMBER_OF_SUB_LOOPS: 0
  SUBROUTINE Problem_ControlLoopCreateStart(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to start the creation of a control for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_ControlLoopCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ASSOCIATED(problem%controlLoop)) CALL FlagError("The problem control loop is already associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_CONTROL_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_START_ACTION
    !Start the problem specific control setup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
    
    EXITS("Problem_ControlLoopCreateStart")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the control loop for a problem. \see OpenCMISS::Iron::cmfe_Problem_ControlLoopDestroy
  SUBROUTINE Problem_ControlLoopDestroy(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to destroy the control for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop

    ENTERS("Problem_ControlLoopDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
    
    CALL ControlLoop_Destroy(problem%controlLoop,err,error,*999)
       
    EXITS("Problem_ControlLoopDestroy")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopDestroy

  !
  !================================================================================================================================
  !

  !>Sets up the specifices for a problem.
  SUBROUTINE Problem_Setup(problem,problemSetupInfo,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetupInfo !<The problem setup information.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_Setup",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      CALL Elasticity_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      CALL BIOELECTRIC_PROBLEM_SETUP(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_FITTING_CLASS)
      CALL Fitting_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_PROBLEM_SETUP(problem,problemSetupInfo,err,error,*999)
    CASE DEFAULT
      localError="The first problem specification of "//TRIM(NumberToVString(problem%specification(1),"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Problem_Setup")
    RETURN
999 ERRORSEXITS("Problem_Setup",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Setup

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear problem solver.
  SUBROUTINE Problem_SolverJacobianEvaluate(solver,err,error,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,solverMatrixIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SOLVER_TYPE), POINTER :: cellMLSolver,linkingSolver
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverJacobianEvaluate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    
    IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) THEN
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver vector values:",err,error,*999)
      DO solverMatrixIdx=1,solverMatrices%NUMBER_OF_MATRICES
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,solverMatrix%SOLVER_VECTOR,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
    !Check if the nonlinear solver is linked to a dynamic solver 
    linkingSolver=>solver%LINKING_SOLVER
    IF(ASSOCIATED(linkingSolver)) THEN
      IF(linkingSolver%SOLVE_TYPE/=SOLVER_DYNAMIC_TYPE) &
        & CALL FlagError("Solver equations linking solver mapping is not dynamic.",err,error,*999)
      !Update the field values from the dynamic factor * current solver values AND add in mean predicted displacements/
      CALL Solver_VariablesDynamicNonlinearUpdate(solver,err,error,*999)
      !check for a linked CellML solver 
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        cellMLSolver=>newtonSolver%CELLML_EVALUATOR_SOLVER
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        cellMLSolver=>quasiNewtonSolver%CELLML_EVALUATOR_SOLVER
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Calculate the Jacobian
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        !Assemble the equations for dynamic problems
        CALL EquationsSet_JacobianEvaluate(equationsSet,err,error,*999)
      ENDDO !equationsSetIdx
      !Assemble the dynamic nonlinear solver matrices
      CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_JACOBIAN_ONLY,err,error,*999)
    ELSE
      !Otherwise perform as steady nonlinear
      !Copy the current solution vector to the dependent field
      CALL SOLVER_VARIABLES_FIELD_UPDATE(solver,err,error,*999)
      !check for a linked CellML solver 
!!TODO: This should be generalised for nonlinear solvers in general and not just Newton solvers.
      newtonSolver=>solver%NONLINEAR_SOLVER%NEWTON_SOLVER
      IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
      cellMLSolver=>newtonSolver%CELLML_EVALUATOR_SOLVER
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Calculate the Jacobian
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        !Assemble the equations for linear problems
        CALL EquationsSet_JacobianEvaluate(equationsSet,err,error,*999)
      ENDDO !equationsSetIdx
      !Update interface matrices
      !DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      !  interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
      !  !Assemble the interface condition for the Jacobian LHS
      !  CALL WriteString(GENERAL_OUTPUT_TYPE,"********************Jacobian evaluation******************",err,error,*999)
      !  CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
      !ENDDO
      !Assemble the static nonlinear solver matrices
      CALL Solver_StaticAssemble(solver,SOLVER_MATRICES_JACOBIAN_ONLY,err,error,*999)
    ENDIF
    
    EXITS("Problem_SolverJacobianEvaluate")
    RETURN
999 ERRORSEXITS("Problem_SolverJacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverJacobianEvaluate
  
  !
  !================================================================================================================================
  ! 

  !>Evaluates the residual for a nonlinear problem solver.
  SUBROUTINE Problem_SolverResidualEvaluate(solver,err,error,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,solverMatrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
    TYPE(SOLVER_TYPE), POINTER :: cellMLSolver,linkingSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix    
    TYPE(VARYING_STRING) :: localError
    
    NULLIFY(cellMLSolver)
    NULLIFY(linkingSolver)

    ENTERS("Problem_SolverResidualEvaluate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    
    IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) THEN
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver vector values:",err,error,*999)
      DO solverMatrixIdx=1,solverMatrices%NUMBER_OF_MATRICES
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,solverMatrix%SOLVER_VECTOR,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
    !Check if the nonlinear solver is linked to a dynamic solver 
    linkingSolver=>solver%LINKING_SOLVER
    IF(ASSOCIATED(linkingSolver)) THEN
      IF(linkingSolver%SOLVE_TYPE/=SOLVER_DYNAMIC_TYPE) &
        & CALL FlagError("Solver equations linking solver mapping is not dynamic.",err,error,*999)
      !Update the field values from the dynamic factor*current solver values AND add in predicted displacements
      CALL Solver_VariablesDynamicNonlinearUpdate(solver,err,error,*999)
      !Caculate the strain field for an CellML evaluator solver
      CALL Problem_PreResidualEvaluate(solver,err,error,*999)
      !check for a linked CellML solver
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        cellMLSolver=>newtonSolver%CELLML_EVALUATOR_SOLVER
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        cellMLSolver=>quasiNewtonSolver%CELLML_EVALUATOR_SOLVER
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Calculate the residual for each element (M, C, K and g)
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        SELECT CASE(equations%linearity)
        CASE(EQUATIONS_LINEAR)
          !Assemble the equations for linear equations
          CALL EquationsSet_Assemble(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR)
          !Evaluate the residual for nonlinear equations
          CALL EquationsSet_ResidualEvaluate(equationsSet,err,error,*999)
        CASE DEFAULT
          localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
      !Assemble the final solver residual.
      CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,err,error,*999)
    ELSE
      !Perform as normal nonlinear solver
      !Copy the current solution vector to the dependent field
      CALL SOLVER_VARIABLES_FIELD_UPDATE(solver,err,error,*999)
      !Caculate the strain field for an CellML evaluator solver
      CALL Problem_PreResidualEvaluate(solver,err,error,*999)
      !check for a linked CellML solver
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        cellMLSolver=>newtonSolver%CELLML_EVALUATOR_SOLVER
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        cellMLSolver=>quasiNewtonSolver%CELLML_EVALUATOR_SOLVER
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Make sure the equations sets are up to date
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        SELECT CASE(equations%linearity)
        CASE(EQUATIONS_LINEAR)
          !Assemble the equations for linear equations
          CALL EquationsSet_Assemble(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR)
          !Evaluate the residual for nonlinear equations
          CALL EquationsSet_ResidualEvaluate(equationsSet,err,error,*999)
        CASE DEFAULT
          localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
      !Note that the linear interface matrices are not required to be updated since these matrices do not change
      !Update interface matrices
      !DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      !  interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
      !  !Assemble the interface condition for the Jacobian LHS
      !  CALL WriteString(GENERAL_OUTPUT_TYPE,"********************Residual evaluation******************",err,error,*999)
      !  CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
      !ENDDO
      !Assemble the solver matrices
      CALL Solver_StaticAssemble(solver,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,err,error,*999)
    END IF
    CALL Problem_PostResidualEvaluate(solver,err,error,*999)
     
    EXITS("Problem_SolverResidualEvaluate")
    RETURN
999 ERRORSEXITS("Problem_SolverResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Pre-evaluates the residual for the solver
  SUBROUTINE Problem_PreResidualEvaluate(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer the solver to pre-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_PreResidualEvaluate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)

    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver pre-residual: ",solver%label,err,error,*999)
    ENDIF
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL Equations_AssertIsFinished(equations,err,error,*999)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)            
        CALL FlagError("Can not pre-evaluate a residual for linear equations.",err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        SELECT CASE(equations%timeDependence)
        CASE(EQUATIONS_STATIC, &
          & EQUATIONS_QUASISTATIC, &
          & EQUATIONS_FIRST_ORDER_DYNAMIC, &
          & EQUATIONS_SECOND_ORDER_DYNAMIC)                        
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            IF(.NOT.ALLOCATED(equationsSet%specification)) &
              & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
            IF(SIZE(equationsSet%specification,1)<1) &
              & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
            SELECT CASE(equationsSet%specification(1))
            CASE(EQUATIONS_SET_ELASTICITY_CLASS)
              CALL Elasticity_FiniteElementPreResidualEvaluate(equationsSet,err,error,*999)
            CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
              CALL FluidMechanics_FiniteElementPreResidualEvaluate(equationsSet,err,error,*999)
            CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
              !Pre residual evaluate not used
            CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
              !Pre residual evaluate not used
            CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
              !Pre residual evaluate not used
            CASE(EQUATIONS_SET_MODAL_CLASS)
              !Pre residual evaluate not used
            CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
              !Pre residual evaluate not used
            CASE DEFAULT
              localError="The first equations set specification of "// &
                & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
              CALL FlagError(localError,err,error,*999)
            END SELECT !equationsSet%specification(1)
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            SELECT CASE(equationsSet%specification(1))
            CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
              !Pre residual evaluate not used
            CASE DEFAULT
              localError="The first equations set specification of "// &
                & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
              CALL FlagError(localError,err,error,*999)
            END SELECT !equationsSet%specification(1)
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
            localError="The equations set solution method  of "// &
              & TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT !equationsSet%solutionMethod
        CASE(EQUATIONS_TIME_STEPPING)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The equations set time dependence type of "// &
            & TRIM(NumberToVString(equations%timeDependence,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !equationsSetIdx
       
    EXITS("Problem_PreResidualEvaluate")
    RETURN
999 ERRORSEXITS("Problem_PreResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_PreResidualEvaluate
     
  !
  !================================================================================================================================
  !

  !>Post-evaluates the residual for the solver
  SUBROUTINE Problem_PostResidualEvaluate(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer the solver to post-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_PostResidualEvaluate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    
    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver post-residual: ",solver%label,err,error,*999)
    ENDIF
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL Equations_AssertIsFinished(equations,err,error,*999)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)            
        CALL FlagError("Can not post-evaluate a residual for linear equations.",err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        SELECT CASE(equations%timeDependence)
        CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC,EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            IF(.NOT.ALLOCATED(equationsSet%specification)) &
              & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
            IF(SIZE(equationsSet%specification,1)<1) &
              & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
            SELECT CASE(equationsSet%specification(1))
            CASE(EQUATIONS_SET_ELASTICITY_CLASS)
              CALL Elasticity_FiniteElementPostResidualEvaluate(equationsSet,err,error,*999)
            CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
              !Post residual evaluate not used
            CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
              !Post residual evaluate not used
            CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
              !Post residual evaluate not used
            CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
              !Post residual evaluate not used
            CASE(EQUATIONS_SET_MODAL_CLASS)
              !Post residual evaluate not used
            CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
              !Post residual evaluate not used
            CASE DEFAULT
              localError="The first equations set specification of "// &
                & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
              CALL FlagError(localError,err,error,*999)
            END SELECT !equationsSet%specification(1)
          CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
            SELECT CASE(equationsSet%specification(1))
            CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
              !Post residual evaluate not used
            CASE DEFAULT
              localError="The first equations set specification of "// &
                & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))// &
                & " is not valid with the nodal solution method."
              CALL FlagError(localError,err,error,*999)
            END SELECT !equationsSet%specification(1)
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
            localError="The equations set solution method  of "// &
              & TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT !equationsSet%solutionMethod
        CASE(EQUATIONS_TIME_STEPPING)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The equations set time dependence type of "// &
            & TRIM(NumberToVString(equations%timeDependence,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity of "// &
          & TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !equationsSetIdx
       
    EXITS("Problem_PostResidualEvaluate")
    RETURN
999 ERRORSEXITS("Problem_PostResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_PostResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Finish the creation of solvers for a problem. \see OpenCMISS::Iron::cmfe_Problem_SolversCreateFinish
  SUBROUTINE Problem_SolversCreateFinish(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finish the creation of the solvers for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: problemSetupInfo
     
    ENTERS("Problem_SolversCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_SOLVERS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish the problem specific solvers setup.
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
       
    EXITS("Problem_SolversCreateFinish")
    RETURN
999 ERRORSEXITS("Problem_SolversCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolversCreateFinish
  
  !
  !================================================================================================================================
  !

  !>Start the creation of a solvers for the problem. \see OpenCMISS::Iron::cmfe_Problem_SolversCreateStart
  SUBROUTINE Problem_SolversCreateStart(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to create the solvers for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_SolversCreateStart",err,error,*999)
    
    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_SOLVERS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_START_ACTION
    !Start the problem specific solvers setup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
    
    EXITS("Problem_SolversCreateStart")
    RETURN
999 ERRORSEXITS("Problem_SolversCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolversCreateStart
  
  !
  !================================================================================================================================
  !

  !>Solves a problem. \see OpenCMISS::Iron::cmfe_Problem_Solve
  SUBROUTINE Problem_Solve(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    
    ENTERS("Problem_Solve",err,error,*999)

    CALL Problem_AssertIsFinished(problem,err,error,*999)

    NULLIFY(controlLoop)
    CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
    CALL ControlLoop_FieldVariablesCalculate(controlLoop,err,error,*999)
    CALL Problem_ControlLoopSolve(controlLoop,err,error,*999)
       
    EXITS("Problem_Solve")
    RETURN
999 ERRORSEXITS("Problem_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_Solve

  !
  !================================================================================================================================
  !

  !> Apply the load increment for each equations_set associated with solver.
  SUBROUTINE Problem_SolverLoadIncrementApply(solverEquations,iterationNumber,maximumNumberOfIterations,err,error,*)
    
    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    INTEGER(INTG) :: equationsSetIdx

    ENTERS("Problem_SolverLoadIncrementApply",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Make sure the equations sets are up to date
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_LoadIncrementApply(equationsSet,solverEquations%BOUNDARY_CONDITIONS,iterationNumber, &
        & maximumNumberOfIterations,err,error,*999)
    ENDDO !equationsSetIdx
    
    EXITS("Problem_SolverLoadIncrementApply")
    RETURN
999 ERRORSEXITS("Problem_SolverLoadIncrementApply",err,error)
    RETURN 1

  END SUBROUTINE Problem_SolverLoadIncrementApply

  !
  !================================================================================================================================
  !

  !>Executes before each loop of a control loop, ie before each time step for a time loop
  SUBROUTINE Problem_ControlLoopPreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_ControlLoopPreLoop",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Control pre-loop: ",controlLoop%label,err,error,*999)
    ENDIF
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    !For all time loops, update the previous values from the current values
    IF(controlLoop%loopType==CONTROL_TIME_LOOP_TYPE) CALL ControlLoop_PreviousValuesUpdate(controlLoop,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      CALL Elasticity_ControlLoopPreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      !do nothing
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_ControlLoopPreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      !do nothing
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      !do nothing
    CASE(PROBLEM_FITTING_CLASS)
      !do nothing
    CASE(PROBLEM_MODAL_CLASS)
      !do nothing
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_CONTROL_LOOP_PRE_LOOP(controlLoop,err,error,*999)
    CASE DEFAULT
      localError="Problem class "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" &
        & is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Problem_ControlLoopPreLoop")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopPreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopPreLoop

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE Problem_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Problem_ControlLoopPostLoop",err,error,*999)

   IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
 
   IF(controlLoop%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
     CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
     CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Control post-loop: ",controlLoop%label,err,error,*999)
   ENDIF
   NULLIFY(problem)
   CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
   IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
   IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
   SELECT CASE(problem%specification(1))
   CASE(PROBLEM_ELASTICITY_CLASS)
     CALL Elasticity_ControlLoopPostLoop(controlLoop,err,error,*999)
   CASE(PROBLEM_BIOELECTRICS_CLASS)
     CALL BIOELECTRIC_CONTROL_LOOP_POST_LOOP(controlLoop,err,error,*999)
   CASE(PROBLEM_FLUID_MECHANICS_CLASS)
     CALL FluidMechanics_ControlLoopPostLoop(controlLoop,err,error,*999)
   CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
     !Do nothing
   CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
     IF(SIZE(problem%specification,1)<2) CALL FlagError("Problem specification must have at least two entries.",err,error,*999)
     CALL ClassicalField_ControlLoopPostLoop(controlLoop,err,error,*999)        
     SELECT CASE(problem%specification(2))
     CASE(PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE)
       CALL REACTION_DIFFUSION_CONTROL_LOOP_POST_LOOP(controlLoop,err,error,*999)
     CASE DEFAULT
       !do nothing
     END SELECT
   CASE(PROBLEM_FITTING_CLASS)
     !Do nothing
   CASE(PROBLEM_MODAL_CLASS)
     !Do nothing
   CASE(PROBLEM_MULTI_PHYSICS_CLASS)
     CALL MULTI_PHYSICS_CONTROL_LOOP_POST_LOOP(controlLoop,err,error,*999)
   CASE DEFAULT
     localError="The first problem specification of "// &
       & TRIM(NumberToVString(problem%specification(1),"*",err,error))//" is not valid."
     CALL FlagError(localError,err,error,*999)
   END SELECT
    
    EXITS("Problem_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopPostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopPostLoop

  !
  !================================================================================================================================
  !

  !>Executes pre solver routines for a problem.
  SUBROUTINE Problem_SolverPreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_SolverPreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) &
      & CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    
    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver pre-solve: ",solver%label,err,error,*999)
    ENDIF
    
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      CALL Elasticity_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      CALL BIOELECTRIC_PRE_SOLVE(solver,err,error,*999)
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      !Do nothing???
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_FITTING_CLASS)
      CALL Fitting_PreSolve(solver,err,error,*999)
    CASE(PROBLEM_MODAL_CLASS)
      !Do nothing???
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_PRE_SOLVE(controlLoop,solver,err,error,*999)
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" &
        & is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Problem_SolverPreSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverPreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverPreSolve

  !
  !================================================================================================================================
  !

  !>Executes post solver routines for a problem.
  SUBROUTINE Problem_SolverPostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverPostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) &
      & CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    
    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver post-solve: ",solver%label,err,error,*999)
    ENDIF
    
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      CALL Elasticity_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      CALL BIOELECTRIC_POST_SOLVE(solver,err,error,*999)
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      !Do nothing???
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)                
      CALL ClassicalField_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_FITTING_CLASS)
      CALL Fitting_PostSolve(solver,err,error,*999)
    CASE(PROBLEM_MODAL_CLASS)
      !Do nothing???
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_POST_SOLVE(controlLoop,solver,err,error,*999)
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Problem_SolverPostSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverPostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverPostSolve

  !
  !================================================================================================================================
  !

  !>Solves solver equations for a problem.
  SUBROUTINE Problem_SolverEquationsSolve(solverEquations,err,error,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverEquationsSolve",err,error,*999)

    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver equations solve: ",solver%label,err,error,*999)
    ENDIF
    SELECT CASE(solverEquations%timeDependence)
    CASE(SOLVER_EQUATIONS_STATIC)
      SELECT CASE(solverEquations%linearity)
      CASE(SOLVER_EQUATIONS_LINEAR)
        CALL Problem_SolverEquationsStaticLinearSolve(solverEquations,err,error,*999)
      CASE(SOLVER_EQUATIONS_NONLINEAR)
        CALL Problem_SolverEquationsStaticNonlinearSolve(solverEquations,err,error,*999)
      CASE DEFAULT
        localError="The solver equations linearity of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_EQUATIONS_QUASISTATIC)
      SELECT CASE(solverEquations%linearity)
      CASE(SOLVER_EQUATIONS_LINEAR)
        CALL Problem_SolverEquationsQuasistaticLinearSolve(solverEquations,err,error,*999)
      CASE(SOLVER_EQUATIONS_NONLINEAR)
        CALL Problem_SolverEquationsQuasistaticNonlinearSolve(solverEquations,err,error,*999)
      CASE DEFAULT
        localError="The solver equations linearity of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(solverEquations%linearity)
      CASE(SOLVER_EQUATIONS_LINEAR)
        CALL Problem_SolverEquationsDynamicLinearSolve(solverEquations,err,error,*999)
      CASE(SOLVER_EQUATIONS_NONLINEAR)
        CALL Problem_SolverEquationsDynamicNonlinearSolve(solverEquations,err,error,*999)
      CASE DEFAULT
        localError="The solver equations linearity of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The solver equations time dependence type of "// &
        & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Problem_SolverEquationsSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsSolve

  !
  !================================================================================================================================
  !

  !>Solves dynamic linear solver equations.
  SUBROUTINE Problem_SolverEquationsDynamicLinearSolve(solverEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    
    ENTERS("Problem_SolverEquationsDynamicLinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)

    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    !Make sure the equations sets are up to date
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Set the solver time
    CALL SOLVER_DYNAMIC_TIMES_SET(solver,currentTime,timeIncrement,err,error,*999)
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values for linear problems
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
    ENDDO !equationsSetIdx
   
    EXITS("Problem_SolverEquationsDynamicLinearSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsDynamicLinearSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsDynamicLinearSolve

  !
  !================================================================================================================================
  !

  !>Solves dynamic nonlinear solver equations.
  SUBROUTINE Problem_SolverEquationsDynamicNonlinearSolve(solverEquations,err,error,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: dynamicSolver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverEquationsDynamicNonlinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      IF(dynamicSolver%restart.OR..NOT.dynamicSolver%SOLVER_INITIALISED) THEN!.OR.dynamicSolver%FSI) THEN
        !If we need to restart or we haven't initialised yet or we have an FSI scheme, make sure the equations sets are up to date
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        SELECT CASE(equations%linearity)
        CASE(EQUATIONS_LINEAR)
          !Assemble the equations
          CALL EquationsSet_Assemble(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR)
          !Evaluate the residuals
          CALL EquationsSet_ResidualEvaluate(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR_BCS)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The equations linearity type of "// &
            & TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Set the solver time
    CALL SOLVER_DYNAMIC_TIMES_SET(solver,currentTime,timeIncrement,err,error,*999)
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    
    EXITS("Problem_SolverEquationsDynamicNonlinearSolve")
    RETURN
999 ERRORS("Problem_SolverEquationsDynamicNonlinearSolve",err,error)
    EXITS("Problem_SolverEquationsDynamicNonlinearSolve")
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsDynamicNonlinearSolve

  !
  !================================================================================================================================
  !

  !>Solves quasistatic linear solver equations.
  SUBROUTINE Problem_SolverEquationsQuasistaticLinearSolve(solverEquations,err,error,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
     
    ENTERS("Problem_SolverEquationsQuasistaticLinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    !Make sure the equations sets are up to date
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the current times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)    
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Solve for the current time
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values for linear problems
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
    ENDDO !equationsSetIdx
   
    EXITS("Problem_SolverEquationsQuasistaticLinearSolve")
    RETURN
999 ERRORS("Problem_SolverEquationsQuasistaticLinearSolve",err,error)
    EXITS("Problem_SolverEquationsQuasistaticLinearSolve")
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsQuasistaticLinearSolve

  !
  !================================================================================================================================
  !

  !>Solves quasistatic nonlinear solver equations.
  SUBROUTINE Problem_SolverEquationsQuasistaticNonlinearSolve(solverEquations,err,error,*)
    
    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Problem_SolverEquationsQuasistaticNonlinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    !Make sure the equations sets are up to date
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    !Update the rhs field variable with residuals or backsubstitute for any linear
    !equations sets
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        CALL EquationsSet_Backsubstitute(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        CALL EquationsSet_NonlinearRHSUpdate(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
      CASE DEFAULT
        localError="The equations linearity type of "// &
          & TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !equationsSetIdx
    
    EXITS("Problem_SolverEquationsQuasistaticNonlinearSolve")
    RETURN
999 ERRORS("Problem_SolverEquationsQuasistaticNonlinearSolve",err,error)
    EXITS("Problem_SolverEquationsQuasistaticNonlinearSolve")
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsQuasistaticNonlinearSolve

  !
  !================================================================================================================================
  !

  !>Solves static linear solver equations.
  SUBROUTINE Problem_SolverEquationsStaticLinearSolve(solverEquations,err,error,*)

   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    
#ifdef TAUPROF
    CHARACTER(12) :: CVAR
    INTEGER :: PHASE(2) = [ 0, 0 ]
    SAVE PHASE
#endif

    ENTERS("Problem_SolverEquationsStaticLinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Make sure the equations sets are up to date
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
#ifdef TAUPROF
      WRITE (CVAR,'(a8,i2)') 'Assemble',equationsSetIdx
      CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
      CALL TAU_PHASE_START(PHASE)
#endif
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
#ifdef TAUPROF
      CALL TAU_PHASE_STOP(PHASE)
#endif
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
#ifdef TAUPROF
      WRITE (CVAR,'(a8,i2)') 'Interface',interfaceConditionIdx
      CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
      CALL TAU_PHASE_START(PHASE)
#endif
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
#ifdef TAUPROF
      CALL TAU_PHASE_STOP(PHASE)
#endif
    ENDDO !interfaceConditionIdx

    !Solve
    CALL Solver_Solve(solver,err,error,*999)

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START('EquationsSet_Backsubstitute()')
#endif
    !Back-substitute to find flux values for linear problems
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
    ENDDO !equationsSetIdx
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP('EquationsSet_Backsubstitute()')
#endif
    
    EXITS("Problem_SolverEquationsStaticLinearSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsStaticLinearSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsStaticLinearSolve
  
  !
  !================================================================================================================================
  !

  !>Solves static nonlinear solver equations.
  SUBROUTINE Problem_SolverEquationsStaticNonlinearSolve(solverEquations,err,error,*)
    
   !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
#ifdef TAUPROF
    CHARACTER(12) :: CVAR
    INTEGER :: PHASE(2) = [ 0, 0 ]
    SAVE PHASE
#endif
    ENTERS("Problem_SolverEquationsStaticNonlinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Apply boundary conditition
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Assemble the equations set
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
#ifdef TAUPROF
      WRITE (CVAR,'(a8,i2)') 'Interface',interfaceConditionIdx
      CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
      CALL TAU_PHASE_START(PHASE)
#endif
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
#ifdef TAUPROF
      CALL TAU_PHASE_STOP(PHASE)
#endif
    ENDDO !interfaceConditionIdx
    !Solve
    CALL Solver_Solve(solver,err,error,*999)
    !Update the rhs field variable with residuals or backsubstitute for any linear
    !equations sets
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR_BCS)
        CALL EquationsSet_Backsubstitute(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
      CASE(EQUATIONS_NONLINEAR)
        CALL EquationsSet_NonlinearRHSUpdate(equationsSet,solverEquations%BOUNDARY_CONDITIONS,err,error,*999)
      CASE DEFAULT
        localError="The equations linearity type of "// &
          & TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !equationsSetIdx
    
    EXITS("Problem_SolverEquationsStaticNonlinearSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsStaticNonlinearSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsStaticNonlinearSolve

  !
  !================================================================================================================================
  !


  !>Solves a solver for a problem.
  SUBROUTINE Problem_SolverSolve(solver,err,error,*)

   !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Problem_SolverSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    IF(solver%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver: ",solver%label,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Solver index = ",solver%globalNumber,err,error,*999)
    ENDIF
      
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START('Pre solve')
#endif
      
    CALL Problem_SolverPreSolve(solver,err,error,*999)
      
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP('Pre solve')
    
    CALL TAU_STATIC_PHASE_START('Solve')
#endif
    
    IF(ASSOCIATED(solver%SOLVER_EQUATIONS)) THEN
      !A solver with solver equations.
      CALL Problem_SolverEquationsSolve(solver%SOLVER_EQUATIONS,err,error,*999)
    ELSE
      !Check for other equations.
      IF(ASSOCIATED(solver%CELLML_EQUATIONS)) THEN
        !A solver with CellML equations.
        CALL Problem_CellMLEquationsSolve(solver%CELLML_EQUATIONS,err,error,*999)
      ELSEIF(solver%SOLVE_TYPE==SOLVER_GEOMETRIC_TRANSFORMATION_TYPE) THEN
        CALL Problem_SolverGeometricTransformationSolve(solver%geometricTransformationSolver,err,error,*999)
      ELSE
        !Do nothing now. 
        !CALL FlagError("Solver does not have any equations associated.",err,error,*999)
      ENDIF
    ENDIF

#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP('Solve')
      
    CALL TAU_STATIC_PHASE_START('Post solve')
#endif
    
    CALL Problem_SolverPostSolve(solver,err,error,*999)
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP('Post solve')
#endif
      
    EXITS("Problem_SolverSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverSolve

  !
  !================================================================================================================================
  !

  !>Destroy the solvers for a problem. \see OpenCMISS::Iron::cmfe_Problem_SolversDestroy
  SUBROUTINE Problem_SolversDestroy(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to destroy the solvers for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop

    ENTERS("Problem_SolversDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
    CALL ControlLoop_SolversDestroy(controlLoop,err,error,*999)
       
    EXITS("Problem_SolversDestroy")
    RETURN
999 ERRORSEXITS("Problem_SolversDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolversDestroy

  !
  !================================================================================================================================
  !

  !>Set boundary conditions for solver equations according to the analytic equations. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsBoundaryConditionsAnalytic
  SUBROUTINE Problem_SolverEquationsBoundaryConditionsAnalytic(solverEquations,err,error,*)

    !Argument variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations to get the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(EquationsSetType), POINTER :: equationsSet

    ENTERS("Problem_SolverEquationsBoundaryConditionsAnalytic",err,error,*999)

    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)

    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_BoundaryConditionsAnalytic(equationsSet,boundaryConditions,err,error,*999)
    ENDDO !equationsSetIdx

    EXITS("Problem_SolverEquationsBoundaryConditionsAnalytic")
    RETURN
999 ERRORS("Problem_SolverEquationsBoundaryConditionsAnalytic",err,error)
    EXITS("Problem_SolverEquationsBoundaryConditionsAnalytic")
    RETURN 1

  END SUBROUTINE Problem_SolverEquationsBoundaryConditionsAnalytic

  !
  !================================================================================================================================
  !

  !>Finish the creation of the solver equations for the problem. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsCreateFinish
  SUBROUTINE Problem_SolverEquationsCreateFinish(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to finish the solver equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_SolverEquationsCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish problem specific startup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
      
    EXITS("Problem_SolverEquationsCreateFinish")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsCreateFinish
  
  !
  !================================================================================================================================
  !

  !>Start the creation of solver equations for a problem. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsCreateStart
  !>The default values of the solver attributes are:
  !>- SOLVE_TYPE: 1 (SOLVER_LINEAR_TYPE)
  !>- OUTPUT_TYPE: 0 (SOLVER_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (SOLVER_SPARSE_MATRICES)
  SUBROUTINE Problem_SolverEquationsCreateStart(problem,err,error,*)

    !Argument variablesg
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to start the creation of the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_SolverEquationsCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_START_ACTION
    !Start the problem specific control setup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
       
    EXITS("Problem_SolverEquationsCreateStart")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsCreateStart

  !
  !================================================================================================================================
  !

  !!TODO: this should be removed - just call the solver equations destroy directly???
  
  !>Destroy the solver equations for a problem. \see OpenCMISS::Iron::cmfe_Problem_SolverEquationsDestroy
  SUBROUTINE Problem_SolverEquationsDestroy(problem,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to destroy the solver equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop

    ENTERS("Problem_SolverEquationsDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Problem_ControlLoopRootGet(problem,controlLoop,err,error,*999)
    CALL ControlLoop_SolverEquationsDestroy(controlLoop,err,error,*999)
       
    EXITS("Problem_SolverEquationsDestroy")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsDestroy
  
  !
  !================================================================================================================================
  !

  !>Solves geometric transformation for a field 
  SUBROUTINE Problem_SolverGeometricTransformationSolve(geometricTransformationSolver,err,error,*) !\todo: Add rotation operations.
    
   !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer to the geometric transformation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldParameterSetType), POINTER :: boundaryConditionsParameterSet
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
    TYPE(ControlLoopSimpleType), POINTER :: simpleLoop
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    INTEGER(INTG) :: componentIdx,versionIdx,derivativeIdx,nodeIdx,noGeomComp
    INTEGER(INTG) :: localNodeNumber,userNodeNumber,incrementIdx,iterationNumber
    REAL(DP) :: nodalParameters(3),nodalParametersTrans(3),transformationMatrix(4,4)
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    LOGICAL :: transformBC=.FALSE.,sameBases=.TRUE.
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverGeometricTransformationSolve",err,error,*999) 
    
    IF(.NOT.ASSOCIATED(geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(geometricTransformationSolver%field)) &
      & CALL FlagError("The field of geometric transformation solver is not associated.",err,error,*999)
    
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(geometricTransformationSolver%field,geometricTransformationSolver%fieldVariableType, &
      & fieldVariable,err,error,*999)
    NULLIFY(boundaryConditionsParameterSet)
    CALL FieldVariable_ParameterSetCheck(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,boundaryConditionsParameterSet, &
      & err,error,*999)
    IF(ASSOCIATED(boundaryConditionsParameterSet)) transformBC=.TRUE. !if the BC is defined on the field variable to be transformed
    
    noGeomComp=SIZE(geometricTransformationSolver%transformationMatrices,1)-1 ! Number of geometric components
    !**********************************************************************************************************************
    !Determine iteration/load increment number 
    IF(geometricTransformationSolver%numberOfIncrements>1) THEN
      solver=>geometricTransformationSolver%solver
      IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
      NULLIFY(controlLoop)
      CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
      SELECT CASE(controlLoop%loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        NULLIFY(simpleLoop)
        CALL ControlLoop_SimpleLoopGet(controlLoop,simpleLoop,err,error,*999)
        iterationNumber=1
      CASE(CONTROL_FIXED_LOOP_TYPE)
        NULLIFY(fixedLoop)
        CALL ControlLoop_FixedLoopGet(controlLoop,fixedLoop,err,error,*999)
        iterationNumber=fixedLoop%iterationNumber
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL FlagError("Geometric transformation for time loop is not implemented.",err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        NULLIFY(whileLoop)
        CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,err,error,*999)
        iterationNumber=whileLoop%iterationNumber
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        NULLIFY(loadIncrementLoop)
        CALL ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*999)
        iterationNumber=loadIncrementLoop%iterationNumber
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(controlLoop%loopType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(iterationNumber>geometricTransformationSolver%numberOfIncrements) THEN
        !If load increment is not specified for that iteration, loop around
        incrementIdx=MOD(iterationNumber-1,geometricTransformationSolver%numberOfIncrements)+1
      ELSE
        incrementIdx=iterationNumber !If load increment is specified for that iteration, use that load increment
      ENDIF
    ELSE
      incrementIdx=1
    ENDIF
    !Determine the transformation matrix to use
    IF(geometricTransformationSolver%arbitraryPath .OR. geometricTransformationSolver%numberOfIncrements==1) THEN
      transformationMatrix(1:noGeomComp+1,1:noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
        & (1:noGeomComp+1,1:noGeomComp+1,incrementIdx)
    ELSE !If need to scale transformation matrix (i.e. transformation applied through several load increment.)
      IF(incrementIdx==1) THEN ! 1st load increment, rotation is applied
        transformationMatrix(1:noGeomComp,1:noGeomComp)=geometricTransformationSolver%transformationMatrices &
          & (1:noGeomComp,1:noGeomComp,1)
      ELSE !No rotation operation in any other load increments
        DO componentIdx=1,noGeomComp
          transformationMatrix(componentIdx,componentIdx)=1.0_DP
        ENDDO !componentIdx
      ENDIF
      !Translation is scaled for every load increment 
      IF(ALLOCATED(geometricTransformationSolver%scalings)) THEN
        transformationMatrix(1:noGeomComp,noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
          & (1:noGeomComp,noGeomComp+1,1)*geometricTransformationSolver%scalings(incrementIdx)
      ELSE !if no scaling just take 1/numberOfIncrements as scaling
        transformationMatrix(1:noGeomComp,noGeomComp+1)=geometricTransformationSolver%transformationMatrices &
          & (1:noGeomComp,noGeomComp+1,1)/geometricTransformationSolver%numberOfIncrements
      ENDIF
    ENDIF
    !**********************************************************************************************************************
    ! Transform the field
    ! Determine if the all components have the same mesh components/ bases
    DO componentIdx=1,noGeomComp-1
      IF(fieldVariable%COMPONENTS(componentIdx)%meshComponentNumber/= &
        & fieldVariable%COMPONENTS(componentIdx+1)%meshComponentNumber) sameBases=.FALSE.
    ENDDO
    IF(sameBases) THEN
      domain=>fieldVariable%components(1)%domain !Use the 1st component domain since they are the same for all components
      IF(ASSOCIATED(domain)) THEN
        domainNodes=>domain%topology%nodes
        DO nodeIdx=1,domainNodes%numberOfNodes
          localNodeNumber=domainNodes%nodes(nodeIdx)%localNumber
          userNodeNumber=domainNodes%nodes(nodeIdx)%userNumber
          DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%numberOfDerivatives
            DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
              DO componentIdx=1,noGeomComp !Get all component for a nodal derivative
                CALL Field_ParameterSetGetNode(geometricTransformationSolver%field,geometricTransformationSolver% &
                  & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                  & nodalParameters(componentIdx),err,error,*999)
              ENDDO !componentIdx
              !Rotate the nodal parameters
              userNodeNumber=domainNodes%nodes(nodeIdx)%userNumber
              nodalParametersTrans(1:noGeomComp)=MATMUL(transformationMatrix(1:noGeomComp,1:noGeomComp), &
                & nodalParameters(1:noGeomComp))
              DO componentIdx=1,noGeomComp !Update all component for a nodal derivative
                CALL Field_ParameterSetUpdateNode(geometricTransformationSolver%field,geometricTransformationSolver% &
                  & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                  & nodalParametersTrans(componentIdx),err,error,*999)
                IF(derivativeIdx==1) THEN ! Translate nodal coordinate
                  CALL Field_ParameterSetAddNode(geometricTransformationSolver%field,geometricTransformationSolver% &
                    & fieldVariableType,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,componentIdx, &
                    & transformationMatrix(componentIdx,1+noGeomComp),err,error,*999)
                ENDIF !derivativeIdx==1
                IF(transformBC) THEN
                  CALL Field_ParameterSetUpdateNode(geometricTransformationSolver%field,geometricTransformationSolver% &
                    & fieldVariableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber, &
                    & componentIdx,nodalParametersTrans(componentIdx),err,error,*999)
                  IF(derivativeIdx==1) THEN ! Translate nodal coordinate for BC
                    CALL Field_ParameterSetAddNode(geometricTransformationSolver%field,geometricTransformationSolver% &
                      & fieldVariableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber, &
                      & componentIdx,transformationMatrix(componentIdx,1+noGeomComp),err,error,*999)
                  ENDIF !derivativeIdx==1
                ENDIF !transformBC
              ENDDO !componentIdx
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ELSE
        CALL FlagError("Domain is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Transformation for different component bases not implemented.",err,error,*999)
    ENDIF
      
    EXITS("Problem_SolverGeometricTransformationSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverGeometricTransformationSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverGeometricTransformationSolve

  !
  !================================================================================================================================
  !

  !>Monitors the problem nonlinear solve
  SUBROUTINE Problem_SolverNonlinearMonitor(solver,iterationNumber,residualNorm,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to monitor
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The number of iterations
    REAL(DP), INTENT(IN) :: residualNorm !<The residual norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceType), POINTER :: interface
    LOGICAL :: reproject
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverNonlinearMonitor",err,error,*998)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      IF(SIZE(problem%specification,1)/=3) &
        & CALL FlagError("Problem specification must have three entries for an elasticity problem.",err,error,*999)
      SELECT CASE(problem%specification(2))
      CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE)
        !Output meshes at iterations
        IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NULLIFY(nonlinearSolver)
          CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
          CALL Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*999)
        ENDIF
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        SELECT CASE(problem%specification(3))
        CASE(PROBLEM_LE_CONTACT_TRANSFORM_SUBTYPE,PROBLEM_FE_CONTACT_TRANSFORM_SUBTYPE) !Reproject at iteration 0 before the nonlinear solve to update xi location since the field is transformed.
          IF(iterationNumber==0) THEN
            reproject=.TRUE.
          ELSE
            reproject=.FALSE.
          ENDIF
        CASE(PROBLEM_LE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_LE_CONTACT_REPROJECT_SUBTYPE, &
          & PROBLEM_FE_CONTACT_TRANSFORM_REPROJECT_SUBTYPE,PROBLEM_FE_CONTACT_REPROJECT_SUBTYPE)
          reproject=.TRUE.
        CASE DEFAULT
          localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))//" &
            & is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(reproject) THEN
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            NULLIFY(interfaceCondition)
            CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
            IF(interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_REPROJECT_OPERATOR .OR. &
              & interfaceCondition%OPERATOR==INTERFACE_CONDITION_FLS_CONTACT_OPERATOR) THEN !Only reproject for contact operator
              IF(interfaceCondition%integrationType==INTERFACE_CONDITION_DATA_POINTS_INTEGRATION) THEN !Only reproject for data point interpolated field
                NULLIFY(INTERFACE)
                CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
                CALL InterfacePointsConnectivity_DataReprojection(INTERFACE,interfaceCondition,err,error,*999)
                CALL INTERFACE_CONDITION_ASSEMBLE(interfaceCondition,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !interfaceConditionIdx
        ENDIF !Reproject
        !Output meshes at iterations
        IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
          NULLIFY(nonlinearSolver)
          CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
          CALL Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The problem type of "//TRIM(NumberToVString(problem%specification(2),"*",err,error))//" &
          & is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
      & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
      !Do nothing???
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" &
        & is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Nonlinear solve monitor--progress output if required
    IF(solver%SOLVE_TYPE==SOLVER_NONLINEAR_TYPE) THEN
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      CALL SOLVER_NONLINEAR_MONITOR(nonlinearSolver,iterationNumber,residualNorm,err,error,*999)
    ELSE
      localError="Invalid solve type. The solve type of "//TRIM(NumberToVString(solver%SOLVE_TYPE,"*",err,error))// &
        & " does not correspond to a nonlinear solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    EXITS("Problem_SolverNonlinearMonitor")
    RETURN
999 NULLIFY(solver)
998 ERRORSEXITS("Problem_SolverNonlinearMonitor",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverNonlinearMonitor
  
  !
  !================================================================================================================================
  !

  !> Output fields at Newton iterations. This is in temporarily for debug output. It may be removed at a later date.
  SUBROUTINE Problem_SolverNewtonFieldsOutput(solver,iterationNumber,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to solver to output the fields for
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<Iteration number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: load_step
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: directory
    
    INTEGER(INTG) :: interfaceConditionIdx, interfaceElementNumber, dataPointIdx, globalDataPointNumber, coupledElementNumber, &
      & coupledMeshFaceLineNumber, coupledMeshIdx,component
    TYPE(InterfaceType), POINTER :: interface !<A pointer to the interface 
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity !<A pointer to the interface points connectivity
    TYPE(FieldType), POINTER :: coupledMeshDependentField
    TYPE(FieldInterpolationParametersPtrType), POINTER :: interpolationParameters(:)
    TYPE(FieldInterpolatedPointPtrType), POINTER :: interpolatedPoints(:)
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData !<A pointer to the decomposition data point topology
    TYPE(DataPointsType), POINTER :: interfaceDatapoints
    TYPE(DataProjectionType), POINTER :: dataProjection

    TYPE(ProblemType), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations

    INTEGER(INTG) :: IUNIT
    CHARACTER(LEN=100) :: filenameOutput,groupname

    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: solve_call

    ENTERS("Problem_SolverNewtonFieldsOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)
    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      IF(SIZE(problem%specification,1)/=3) &
        & CALL FlagError("Problem specification must have three entries for an elasticity problem.",err,error,*999)
      SELECT CASE(problem%SPECIFICATION(2))
      CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE,PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE, &
        & PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)

        !This is not how diagnostics should be used
        ! IF(diagnostics1) THEN          
        !   directory="results_iter/"
        !   INQUIRE(FILE=CHAR(directory),EXIST=dirExists)
        !   IF(.NOT.dirExists) THEN
        !     CALL SYSTEM(CHAR("mkdir "//directory))
        !   ENDIF
        
        !   ! Find how many times the problem solve command has been issued.
        !   max_solve_calls=100
        !   coupledMeshIdx=1
        !   load_step=1
        !   firstIterationNumber=0
        !   DO solve_call=1,max_solve_calls
        !     fileToCheck=directory// &
        !       & "mesh"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))// &
        !       & "_solveCall"//TRIM(NumberToVString(solve_call,"*",err,error))// &
        !       & "_load"//TRIM(NumberToVString(load_step,"*",err,error))// &
        !       & "_iter"//TRIM(NumberToVString(firstIterationNumber,"*",err,error))//".part0.exnode"
        !     INQUIRE(FILE=CHAR(fileToCheck),EXIST=fileExists)
        !     IF(.NOT.fileExists) THEN
        !       EXIT
        !     ENDIF
        !   ENDDO
        
        !   load_step=solver%SOLVERS%controlLoop%loadIncrementLoop%iterationNumber
        
        !   IF((iterationNumber > 0).OR.(load_step > 1))THEN
        !     solve_call = solve_call - 1
        !   ENDIF
        
        !   WRITE(*,'(1X,''SolveCall: '',I4)') solve_call
        !   WRITE(*,'(1X,''  LoadStep: '',I4)') load_step
        !   WRITE(*,'(1X,''    Iteration: '',I4)') iterationNumber
        
        !   DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        !     region=>solverMapping%equationsSets(equationsSetIdx)%ptr%REGION
        !     IF(ASSOCIATED(region))THEN
        !       NULLIFY(fields)
        !       fields=>region%FIELDS
        !       fileName=directory//"mesh"//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        !         & "_solveCall"//TRIM(NumberToVString(solve_call,"*",err,error))// &
        !         & "_load"//TRIM(NumberToVString(load_step,"*",err,error))// &
        !         & "_iter"//TRIM(NumberToVString(iterationNumber,"*",err,error))
        !       method="FORTRAN"
        !       CALL FIELD_IO_ELEMENTS_EXPORT(fields,fileName,method,err,error,*999)
        !       CALL FIELD_IO_NODES_EXPORT(fields,fileName,method,err,error,*999)
        !     ELSE
        !       CALL FlagError("Region is not associated.",err,error,*999)
        !     ENDIF
        !   ENDDO
        ! ENDIF

      CASE DEFAULT
        localError="The problem type of "//TRIM(NumberToVString(problem%specification(2),"*",err,error))//" &
          & is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
      & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
      !Do nothing???
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(problem%SPECIFICATION(1),"*",err,error))//" &
        & is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    SELECT CASE(problem%SPECIFICATION(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      SELECT CASE(problem%SPECIFICATION(2))
      CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE)
        ! Pass
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        
        IF(diagnostics1) THEN
          IUNIT = 300
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
            interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
            INTERFACE=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr%INTERFACE
            pointsConnectivity=>interface%pointsConnectivity
            interfaceDatapoints=>pointsConnectivity%dataPoints
            IF(ASSOCIATED(pointsConnectivity)) THEN
              DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
                filenameOutput=directory//"PointsConnectivity"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))// &
                  & "_solveCall"//TRIM(NumberToVString(solve_call,"*",err,error))// &
                  & "_load"//TRIM(NumberToVString(load_step,"*",err,error))// &
                  & "_iter"//TRIM(NumberToVString(iterationNumber,"*",err,error))//".exdata"
                OPEN(UNIT=IUNIT,FILE=filenameOutput,STATUS="UNKNOWN",ACTION="WRITE",IOSTAT=err)
                groupname="PointsConnectivity"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))
                WRITE(IUNIT,'( '' Group name: '',A)') groupname
                WRITE(IUNIT,'(1X,''#Fields=4'')')
                WRITE(IUNIT,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
                WRITE(IUNIT,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''2) error, field, rectangular cartesian, #Components=3'')')
                WRITE(IUNIT,'(1X,''  x.  Value index= 4, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  y.  Value index= 5, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  z.  Value index= 6, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''3) projectedCoordinate, field, rectangular cartesian, #Components=3'')')
                WRITE(IUNIT,'(1X,''  x.  Value index= 7, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  y.  Value index= 8, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''  z.  Value index= 9, #Derivatives=0'')')
                WRITE(IUNIT,'(1X,''4) exitTag, field, rectangular cartesian, #Components=1'')')
                WRITE(IUNIT,'(1X,''  tag.  Value index= 10, #Derivatives=0'')')
                coupledMeshDependentField=>interfaceCondition%dependent%equationsSets(coupledMeshIdx)%ptr% &
                  & dependent%dependentField
                NULLIFY(interpolationParameters)
                CALL Field_InterpolationParametersInitialise(coupledMeshDependentField,interpolationParameters,err,error, &
                  & *999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                NULLIFY(interpolatedPoints)
                CALL Field_InterpolatedPointsInitialise(interpolationParameters,interpolatedPoints,err,error,*999, &
                  & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                interpolatedPoint=>interpolatedPoints(FIELD_U_VARIABLE_TYPE)%ptr
                dataProjection=>interfaceDatapoints%dataProjections%dataProjections(coupledMeshIdx+1)%ptr
                DO interfaceElementNumber=1,SIZE(pointsConnectivity%coupledElements,1)
                  decompositionElementData=>interfaceCondition%LAGRANGE%lagrangeField%DECOMPOSITION%TOPOLOGY%dataPoints% &
                    & elementDataPoints(interfaceElementNumber)
                  DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                    globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
                    WRITE(IUNIT,'(1X,''Node:'',I4)') globalDataPointNumber
                    DO component=1,3
                      WRITE(IUNIT,'(1X,3E25.15)') interfaceDatapoints%dataPoints(globalDataPointNumber)%position(component)
                    ENDDO !component
                    coupledElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                      & coupledElementNumber
                    coupledMeshFaceLineNumber=coupledMeshDependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS% &
                      & ELEMENTS(coupledElementNumber)% &
                      & elementFaces(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                      & elementLineFaceNumber)
                    CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
                      & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                    CALL Field_InterpolateXi(NO_PART_DERIV,pointsConnectivity%pointsConnectivity(globalDataPointNumber, &
                      & coupledMeshIdx)%reducedXi(:),interpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                    DO component=1,3
                      WRITE(IUNIT,'(1X,3E25.15)') interpolatedPoint%VALUES(component,NO_PART_DERIV) - &
                        & interfaceDatapoints%dataPoints(globalDataPointNumber)%position(component)
                    ENDDO !component
                    DO component=1,3
                      WRITE(IUNIT,'(1X,3E25.15)') interpolatedPoint%VALUES(component,NO_PART_DERIV)
                    ENDDO !component
                    WRITE(IUNIT,'(1X,I2)') dataProjection%dataProjectionResults(globalDataPointNumber)%exitTag
                  ENDDO !dataPointIdx
                ENDDO !interfaceElementNumber
                CALL Field_InterpolationParametersFinalise(interpolationParameters,err,error,*999)
                CALL Field_InterpolatedPointsFinalise(interpolatedPoints,err,error,*999)
                OPEN(UNIT=IUNIT)
              ENDDO !coupledMeshIdx
            ENDIF
          ENDDO !interfaceConditionIdx
        ENDIF
        
      CASE DEFAULT
        localError="The problem type of "//TRIM(NumberToVString(problem%SPECIFICATION(2),"*",err,error))//" &
          & is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
      & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
      !Do nothing???
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(problem%SPECIFICATION(1),"*",err,error))//" &
        & is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Problem_SolverNewtonFieldsOutput")
    RETURN
999 ERRORSEXITS("Problem_SolverNewtonFieldsOutput",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverNewtonFieldsOutput
  
  !
  !================================================================================================================================
  !

  !>Monitors the problem optimiser solve
  SUBROUTINE Problem_SolverOptimiserMonitor(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver to monitor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver
    
    ENTERS("Problem_SolverOptimiserMonitor",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    
    !IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    !IF(SIZE(problem%specification,1)<1) CALL FlagError("Problem specification must have at least one entry.",err,error,*999)

    !SELECT CASE(problem%specification(1))
    !CASE DEFAULT
    !  localError="The problem class of "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" is invalid."
    !  CALL FlagError(localError,err,error,*999)
    !END SELECT

    !Optimiser solve monitor--progress output if required
    CALL Solver_AssertIsOptimiser(solver,err,error,*999)
    NULLIFY(optimiserSolver)
    CALL Solver_OptimiserSolverGet(solver,optimiserSolver,err,error,*999)
    CALL Solver_OptimiserMonitor(optimiserSolver,err,error,*999)
    
    EXITS("Problem_SolverOptimiserMonitor")
    RETURN
999 ERRORSEXITS("Problem_SolverOptimiserMonitor",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverOptimiserMonitor
  
  !
  !================================================================================================================================
  !

  !>Gets the problem specification array for a problem identified by a pointer. \see OpenCMISS::Iron::cmfe_Problem_SpecificationGet
  SUBROUTINE Problem_SpecificationGet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the specification for.
    INTEGER(INTG), INTENT(INOUT) :: problemSpecification(:) !<On return, The problem specifcation array. Must be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: specificationLength
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_SpecificationGet",err,error,*999)

    CALL Problem_AssertIsFinished(problem,err,error,*999)

    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    specificationLength=SIZE(problem%specification,1)
    IF(SIZE(problemSpecification,1)<specificationLength) THEN
      localError="The problem specification size is "//TRIM(NumberToVstring(specificationLength,"*",err,error))// &
        & " but the input array only has size "//TRIM(NumberToVstring(SIZE(problemSpecification,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    problemSpecification(1:specificationLength)=problem%specification(1:specificationLength)

    EXITS("Problem_SpecificationGet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationGet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification
  SUBROUTINE Problem_SpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for.
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification array to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemClass

    ENTERS("Problem_SpecificationSet",err,error,*999)

    CALL Problem_AssertNotFinished(problem,err,error,*999)
    IF(SIZE(problemSpecification,1)<1) CALL FlagError("Problem specification array must have one or more entries.",err,error,*999)
    
    problemClass=problemSpecification(1)
    SELECT CASE(problemClass)
    CASE(PROBLEM_ELASTICITY_CLASS)
      CALL Elasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      CALL Bioelectric_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_FITTING_CLASS)
      CALL Fitting_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE(PROBLEM_OPTIMISATION_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_ProblemSpecificationSet(problem,problemSpecification,err,error,*999)
    CASE DEFAULT
      localError="The first problems specification of "//TRIM(NumberToVstring(problemClass,"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Problem_SpecificationSet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationSet

  !
  !================================================================================================================================
  !

  !>Gets the size of the problem specification array for a problem identified by a pointer. \see OpenCMISS::Iron::cmfe_Problem_SpecificationSizeGet
  SUBROUTINE Problem_SpecificationSizeGet(problem,specificationSize,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: specificationSize !<On return, the size of the problem specifcation array.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_SpecificationSizeGet",err,error,*999)

    specificationSize=0
    CALL Problem_AssertIsFinished(problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    
    specificationSize=SIZE(problem%specification,1)

    EXITS("Problem_SpecificationSizeGet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationSizeGet

  !
  !================================================================================================================================
  !

  !>Sets the work group for a problem. \see OpenCMISS::Iron::cmfe_Problem_WorkGroupGet
  SUBROUTINE Problem_WorkGroupSet(problem,workGroup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER, INTENT(INOUT) :: problem !<A pointer to the problem to set the work group for.
    TYPE(WorkGroupType), POINTER, INTENT(IN) :: workGroup !<A pointer to the workgroup to set for the problem.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problem_WorkGroupSet",err,error,*999)

    CALL Problem_AssertNotFinished(problem,err,error,*999)
    CALL WorkGroup_AssertIsFinished(workGroup,err,error,*999)
 
    problem%workGroup=>workGroup

    EXITS("Problem_WorkGroupSet")
    RETURN
999 ERRORSEXITS("Problem_WorkGroupSet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_WorkGroupSet

  !
  !================================================================================================================================
  !

  !>Finalises all problems and deallocates all memory.
  SUBROUTINE Problems_Finalise(problems,err,error,*)

    !Argument variables
    TYPE(ProblemsType), POINTER :: problems !<A pointer to the problems to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Problems_Finalise",err,error,*999)

    IF(ASSOCIATED(problems)) THEN
      DO WHILE(problems%numberOfProblems>0)
        CALL Problem_Destroy(problems%problems(1)%ptr,err,error,*999)
      ENDDO !problemIdx
      DEALLOCATE(problems)
    ENDIF
    
    EXITS("Problems_Finalise")
    RETURN
999 ERRORSEXITS("Problems_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Problems_Finalise

  !
  !================================================================================================================================
  !

  !>Intialises all problems for a context.
  SUBROUTINE Problems_Initialise(context,err,error,*)

    !Argument variables
    TYPE(ContextType), POINTER :: context !<The context to intialise the problems for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Problems_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(context)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(ASSOCIATED(context%problems)) CALL FlagError("Context problems is already associated.",err,error,*998)
    
    ALLOCATE(context%problems,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate problems.",err,error,*999)
    !Initialise
    context%problems%context=>context
    context%problems%numberOfProblems=0
    NULLIFY(context%problems%problems)
    
    EXITS("Problems_Initialise")
    RETURN
999 CALL Problems_Finalise(context%problems,dummyErr,dummyError,*998)
998 ERRORSEXITS("Problems_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Problems_Initialise
  
  !
  !================================================================================================================================
  !

  
END MODULE ProblemRoutines

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the Jacobian for a Newton like nonlinear solver
SUBROUTINE Problem_SolverJacobianEvaluatePetsc(snes,x,A,B,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE DistributedMatrixVector
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SOLVER_MATRICES_ROUTINES
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
 
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc snes
  TYPE(PetscVecType), INTENT(INOUT) :: X !<The PETSc x Vec
  TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(DistributedVectorType), POINTER :: solverVector
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*998)

  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*999)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
  IF(solverMatrices%NUMBER_OF_MATRICES/=1) THEN
    localError="The number of solver matrices of "// &
      & TRIM(NumberToVString(solverMatrices%NUMBER_OF_MATRICES,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*998)
  ENDIF
  
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverVectorGet(solverMatrix,solverVector,err,error,*999)
  
  CALL DistributedVector_OverrideSetOn(solverVector,x,err,error,*999)
  
  CALL Problem_SolverJacobianEvaluate(ctx,err,error,*999)
  
  CALL DistributedVector_OverrideSetOff(solverVector,err,error,*999)
  
!!TODO: move this to Problem_SolverJacobianEvaluate or elsewhere?
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)
  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*999)
  SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
  CASE(SOLVER_NONLINEAR_NEWTON)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    newtonSolver%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS=newtonSolver%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS+1
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    quasiNewtonSolver%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS=quasiNewtonSolver%TOTAL_NUMBER_OF_JACOBIAN_EVALUATIONS+1
  CASE DEFAULT
    localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))// &
      & " is invalid."
    CALL FlagError(localError,err,error,*999)
  END SELECT
  
  RETURN
999 CALL DistributedVector_OverrideSetOff(solverVector,dummyErr,dummyError,*998)
998 CALL WriteError(err,error,*997)
997 CALL FlagWarning("Error evaluating nonlinear Jacobian.",err,error,*996)
996 RETURN
  
END SUBROUTINE Problem_SolverJacobianEvaluatePetsc

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the Jacobian for a Newton like nonlinear solver using PETSc's FD Jacobian
!>calculation.
SUBROUTINE Problem_SolverJacobianFDCalculatePetsc(snes,x,A,B,ctx,err)

  USE BaseRoutines
  USE CmissPetsc
  USE CmissPetscTypes
  USE DistributedMatrixVector
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SOLVER_MATRICES_ROUTINES
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE

  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec
  TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: linesearchSolver
  TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
  TYPE(QUASI_NEWTON_LINESEARCH_SOLVER_TYPE), POINTER :: quasiNewtonLinesearchSolver
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
  TYPE(PetscMatFDColoringType), POINTER :: jacobianMatFDColoring
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*998)
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)
  
  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*999)
  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*999)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
  IF(solverMatrices%NUMBER_OF_MATRICES/=1) THEN
    localError="The number of solver matrices of "// &
      & TRIM(NumberToVString(solverMatrices%NUMBER_OF_MATRICES,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*998)
  ENDIF
  
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)  
  SELECT CASE(solverEquations%sparsityType)
  CASE(SOLVER_SPARSE_MATRICES)
    SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      linesearchSolver=>newtonSolver%LINESEARCH_SOLVER
      IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Newton solver linesearch solver is not associated.",err,error,*999)
      jacobianMatFDColoring=>linesearchSolver%jacobianMatFDColoring
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      quasiNewtonLinesearchSolver=>quasiNewtonSolver%LINESEARCH_SOLVER
      IF(.NOT.ASSOCIATED(quasiNewtonLinesearchSolver)) &
        & CALL FlagError("Quasi-Newton solver linesearch solver is not associated.",err,error,*999)
      jacobianMatFDColoring=>quasiNewtonLinesearchSolver%jacobianMatFDColoring
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(.NOT.ASSOCIATED(jacobianMatFDColoring)) CALL FlagError("Linesearch solver FD colouring is not associated.",err,error,*998)
    CALL Petsc_SnesComputeJacobianDefaultColor(snes,x,A,B,jacobianMatFDColoring,err,error,*999)
  CASE(SOLVER_FULL_MATRICES)
    CALL Petsc_SnesComputeJacobianDefault(snes,x,A,B,ctx,err,error,*999)
  CASE DEFAULT
    localError="The specified solver equations sparsity type of "// &
      & TRIM(NumberToVString(solverEquations%sparsityType,"*",err,error))//" is invalid."
    CALL FlagError(localError,err,error,*999)
  END SELECT
  IF(ctx%outputType>=SOLVER_MATRIX_OUTPUT) THEN
    CALL DistributedMatrix_OverrideSetOn(solverMatrices%matrices(1)%ptr%matrix,A,err,error,*999)
    CALL SOLVER_MATRICES_OUTPUT(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES_JACOBIAN_ONLY,solverMatrices,err,error,*998)
    CALL DistributedMatrix_OverrideSetOff(solverMatrices%matrices(1)%ptr%matrix,err,error,*999)
  ENDIF

  RETURN
999 CALL DistributedMatrix_OverrideSetOff(solverMatrix%matrix,dummyErr,dummyError,*998)
998 CALL WriteError(err,error,*997)
997 CALL FlagWarning("Error evaluating nonlinear Jacobian.",err,error,*996)
996 RETURN
  
END SUBROUTINE Problem_SolverJacobianFDCalculatePetsc

!
!================================================================================================================================
!

!>Called from the PETSc TAO solvers to evaluate the objective for an optimiser solver
SUBROUTINE Problem_SolverObjectiveEvaluatePetsc(tao,x,f,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE DistributedMatrixVector
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The PETSc tao type
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec type
  REAL(DP), INTENT(OUT) :: f !<On exit, the evaluated objective
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(DistributedVectorType), POINTER :: solverVector
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*997)
  CALL Solver_AssertIsOptimiser(ctx,err,error,*999)
  
  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*999)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)  
  IF(solverMatrices%NUMBER_OF_MATRICES/=1) THEN
    localError="The number of solver matrices of "//TRIM(NumberToVString(solverMatrices%NUMBER_OF_MATRICES,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*997)          
  ENDIF
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverVectorGet(solverMatrix,solverVector,err,error,*999)

  CALL DistributedVector_OverrideSetOn(solverVector,x,err,error,*999)
    
  !CALL Problem_SolverObjectiveEvaluate(ctx,err,error,*999)

  f=0.0_DP
                    
  CALL DistributedVector_OverrideSetOff(solverVector,err,error,*999)

!!TODO: move this to Problem_SolverResidualEvaluate or elsewhere?
  !optimiserSolver%totalNumberOfObjectiveEvaluations=optimiserSolver%totalNumberOfFunctionEvaluations+1
  
  RETURN
999 CALL DistributedVector_OverrideSetOff(solverVector,dummyErr,dummyError,*998)  
998 CALL WriteError(err,error,*997)
997 CALL FlagWarning("Error evaluating optimiser objective.",err,error,*996)
996 RETURN

END SUBROUTINE Problem_SolverObjectiveEvaluatePetsc

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to evaluate the residual for a Newton like nonlinear solver
SUBROUTINE Problem_SolverResidualEvaluatePetsc(snes,x,f,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE DistributedMatrixVector
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc snes type
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec type
  TYPE(PetscVecType), INTENT(INOUT) :: f !<The PETSc f Vec type
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(DistributedVectorType), POINTER :: residualVector,solverVector
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
  TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
  TYPE(SOLVER_MATRICES_TYPE), POINTER :: solverMatrices
  TYPE(SOLVER_MATRIX_TYPE), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*997)
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)

  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*997)
  SELECT CASE(nonLinearSolver%NONLINEAR_SOLVE_TYPE)
  CASE(SOLVER_NONLINEAR_NEWTON)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*997)
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*997)
  CASE DEFAULT
    localError="The nonlinear solver type of "//TRIM(NumberToVString(nonLinearSolver%NONLINEAR_SOLVE_TYPE,"*",err,error))// &
      & " is invalid."
    CALL FlagError(localError,err,error,*997)
  END SELECT
  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*997)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*997)
  IF(solverMatrices%NUMBER_OF_MATRICES/=1) THEN
    localError="The number of solver matrices of "// &
      & TRIM(NumberToVString(solverMatrices%NUMBER_OF_MATRICES,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*997)          
  ENDIF
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*997)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverVectorGet(solverMatrix,solverVector,err,error,*997)
  NULLIFY(residualVector)
  CALL SolverMatrices_ResidualVectorGet(solverMatrices,residualVector,err,error,*997)
  
  CALL DistributedVector_OverrideSetOn(solverVector,X,err,error,*998)
  CALL DistributedVector_OverrideSetOn(residualVector,F,err,error,*999)                
                    
  CALL Problem_SolverResidualEvaluate(ctx,err,error,*999)
  
  CALL DistributedVector_OverrideSetOff(residualVector,err,error,*999)
  CALL DistributedVector_OverrideSetOff(solverVector,err,error,*998)
  
!!TODO: move this to Problem_SolverResidualEvaluate or elsewhere?
  SELECT CASE(nonLinearSolver%NONLINEAR_SOLVE_TYPE)
  CASE(SOLVER_NONLINEAR_NEWTON)
    newtonSolver%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS=newtonSolver%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS+1
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    quasiNewtonSolver%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS=quasiNewtonSolver%TOTAL_NUMBER_OF_FUNCTION_EVALUATIONS+1
  CASE DEFAULT
    !Do nothing?
  END SELECT
  
  RETURN
999 CALL DistributedVector_OverrideSetOff(residualVector,dummyErr,dummyError,*998)
998 CALL DistributedVector_OverrideSetOff(solverVector,dummyErr,dummyError,*997)  
997 CALL WriteError(err,error,*996)
996 CALL FlagWarning("Error evaluating nonlinear residual.",err,error,*995)
995 RETURN

END SUBROUTINE Problem_SolverResidualEvaluatePetsc

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to test convergence for a Newton like nonlinear solver
SUBROUTINE Problem_SolverConvergenceTestPetsc(snes,iterationNumber,xnorm,gnorm,fnorm,reason,ctx,err)

  USE BaseRoutines
  USE CmissPetsc
  USE CmissPetscTypes
  USE Constants
  USE DistributedMatrixVector
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE Strings
  USE Types
 
  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES type
  INTEGER(INTG), INTENT(INOUT) :: iterationNumber !< The current iteration (1 is the first and is before any Newton step)
  REAL(DP), INTENT(INOUT) :: xnorm !<The 2-norm of current iterate
  REAL(DP), INTENT(INOUT) :: gnorm !<The 2-norm of current step
  REAL(DP), INTENT(INOUT) :: fnorm !<The 2-norm of function
  INTEGER(INTG), INTENT(INOUT) :: reason !<The reason for convergence/divergence
  TYPE(SOLVER_TYPE), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(PetscVecType) :: x,f,y,w,g
  TYPE(NEWTON_SOLVER_TYPE), POINTER :: newtonSolver
  TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
  TYPE(QUASI_NEWTON_SOLVER_TYPE), POINTER :: quasiNewtonSolver
  TYPE(PetscSnesLinesearchType) :: lineSearch
  REAL(DP) :: energy,normalisedEnergy
  TYPE(VARYING_STRING) :: error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*999)
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)

  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*999)
  
  SELECT CASE(nonlinearSolver%NONLINEAR_SOLVE_TYPE)
  CASE(SOLVER_NONLINEAR_NEWTON)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    reason=PETSC_SNES_CONVERGED_ITERATING
    SELECT CASE(newtonSolver%convergenceTestType)
    CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
      IF(iterationNumber>0) THEN
        CALL Petsc_SnesLineSearchInitialise(lineSearch,err,error,*999)
        CALL Petsc_SnesGetLineSearch(snes,lineSearch,err,error,*999)
        CALL Petsc_VecInitialise(x,err,error,*999)
        CALL Petsc_VecInitialise(f,err,error,*999)
        CALL Petsc_VecInitialise(y,err,error,*999)
        CALL Petsc_VecInitialise(w,err,error,*999)
        CALL Petsc_VecInitialise(g,err,error,*999)
        CALL Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*999)
        CALL Petsc_VecDot(y,g,energy,err,error,*999)
        IF(iterationNumber==1) THEN
          IF(ABS(energy)<ZERO_TOLERANCE) THEN
            reason=PETSC_SNES_CONVERGED_FNORM_ABS
          ELSE
            newtonSolver%convergenceTest%energyFirstIter=energy
            newtonSolver%convergenceTest%normalisedEnergy=1.0
          ENDIF
        ELSE
          normalisedEnergy=energy/newtonSolver%convergenceTest%energyFirstIter
          newtonSolver%convergenceTest%normalisedEnergy=normalisedEnergy
          IF(ABS(normalisedEnergy)<newtonSolver%ABSOLUTE_TOLERANCE) THEN
            reason=PETSC_SNES_CONVERGED_FNORM_ABS
            newtonSolver%convergenceTest%energyFirstIter=0.0_DP
            newtonSolver%convergenceTest%normalisedEnergy=0.0_DP
          ENDIF
          CALL WriteString(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Normalised energy = ",normalisedEnergy,err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
        ENDIF
        CALL Petsc_SnesLineSearchFinalise(lineSearch,err,error,*999)
      ENDIF
    CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
      CALL FlagError("Differentiated ratio convergence test not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The specified convergence test type of "//TRIM(NumberToVString( &
        & newtonSolver%convergenceTestType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    reason=PETSC_SNES_CONVERGED_ITERATING
    SELECT CASE(quasiNewtonSolver%convergenceTestType)
    CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
      IF(iterationNumber>0) THEN
        CALL Petsc_SnesLineSearchInitialise(lineSearch,err,error,*999)
        CALL Petsc_SnesGetLineSearch(snes,lineSearch,err,error,*999)
        CALL Petsc_VecInitialise(x,err,error,*999)
        CALL Petsc_VecInitialise(f,err,error,*999)
        CALL Petsc_VecInitialise(y,err,error,*999)
        CALL Petsc_VecInitialise(w,err,error,*999)
        CALL Petsc_VecInitialise(g,err,error,*999)
        CALL Petsc_SnesLineSearchGetVecs(lineSearch,x,f,y,w,g,err,error,*999)
        CALL Petsc_VecDot(y,g,energy,err,error,*999)
        IF(iterationNumber==1) THEN
          IF(ABS(energy)<ZERO_TOLERANCE) THEN
            reason=PETSC_SNES_CONVERGED_FNORM_ABS
          ELSE
            quasiNewtonSolver%convergenceTest%energyFirstIter=energy
            quasiNewtonSolver%convergenceTest%normalisedEnergy=1.0
          ENDIF
        ELSE
          normalisedEnergy=energy/quasiNewtonSolver%convergenceTest%energyFirstIter
          quasiNewtonSolver%convergenceTest%normalisedEnergy=normalisedEnergy
          IF(ABS(normalisedEnergy)<quasiNewtonSolver%ABSOLUTE_TOLERANCE) THEN
            reason=PETSC_SNES_CONVERGED_FNORM_ABS
            quasiNewtonSolver%convergenceTest%energyFirstIter=0.0_DP
            quasiNewtonSolver%convergenceTest%normalisedEnergy=0.0_DP
          ENDIF
          CALL WriteString(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Normalised energy = ",normalisedEnergy,err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"*********************************************",err,error,*999)
        ENDIF
        CALL Petsc_SnesLineSearchFinalise(lineSearch,err,error,*999)
      ELSE
        quasiNewtonSolver%convergenceTest%energyFirstIter=0.0_DP
        quasiNewtonSolver%convergenceTest%normalisedEnergy=0.0_DP
      ENDIF
    CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
      CALL FlagError("Differentiated ratio convergence test not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The specified convergence test type of "//TRIM(NumberToVString( &
        & quasiNewtonSolver%convergenceTestType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
  CASE DEFAULT
    !Do nothing?
  END SELECT
  
  RETURN
999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error in convergence test.",err,error,*997)
997 RETURN    

END SUBROUTINE Problem_SolverConvergenceTestPetsc

!
!================================================================================================================================
!


!>Called from the PETSc TS solvers to solve cellml DAE
SUBROUTINE Problem_SolverDAECellMLRHSPetsc(ts,time,states,rates,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE CmissPetsc
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE Types

  IMPLICIT NONE

  !Argument variables
  TYPE(PetscTSType), INTENT(INOUT) :: ts !<The PETSc TS type
  REAL(DP), INTENT(INOUT) :: time !<The current time
  TYPE(PetscVecType), INTENT(INOUT) :: states !<current states
  TYPE(PetscVecType), INTENT(INOUT) :: rates !<returned rates
  TYPE(CellMLPETScContextType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(CELLML_TYPE), POINTER :: cellML
  TYPE(SOLVER_TYPE), POINTER :: solver
  TYPE(VARYING_STRING) :: error
  INTEGER(INTG) :: dofIdx
  REAL(DP), POINTER :: stateData(:)

  NULLIFY(stateData)

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Context is not associated.",err,error,*999)
  solver=>ctx%solver
  IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Context solver is not associated.",err,error,*999)
  cellML=>ctx%cellml
  IF(.NOT.ASSOCIATED(cellml)) CALL FlagError("Context cellml is not associated.",err,error,*999)  
  dofIdx=ctx%dofIdx
  
  !Get the state data
  NULLIFY(stateData)
  CALL Petsc_VecGetArrayReadF90(states,stateData,err,error,*999)
  !Evaluate the CellML model
  CALL Problem_SolverDAECellMLRHSEvaluate(cellML,time,dofIdx,stateData,ctx%rates,err,error,*999)
  !Restore the state data
  CALL Petsc_VecRestoreArrayReadF90(states,stateData,err,error,*999)
  !Set the PETSc rates vector
  CALL Petsc_VecSetValues(rates,SIZE(stateData,1),ctx%ratesIndices,ctx%rates,PETSC_INSERT_VALUES,err,error,*999)
  CALL VecAssemblyBegin(rates,err,error,*999)
  CALL VecAssemblyEnd(rates,err,error,*999)
  
  RETURN
999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error calling Problem_SolverDAECellMLRHSPetsc routine from PETSc.",err,error,*997)
997 RETURN    

END SUBROUTINE Problem_SolverDAECellMLRHSPetsc


!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to monitor a nonlinear solver
SUBROUTINE Problem_SolverNonlinearMonitorPETSC(snes,iterationNumber,residualNorm,context,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE DistributedMatrixVector
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SolverAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc snes type
  INTEGER(INTG), INTENT(INOUT) :: iterationNumber !<The iteration number
  REAL(DP), INTENT(INOUT) :: residualNorm !<The residual norm
  TYPE(SOLVER_TYPE), POINTER :: context !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(VARYING_STRING) :: error

  IF(.NOT.ASSOCIATED(context)) CALL FlagError("Solver context is not associated.",err,error,*999)
  CALL Solver_AssertIsNonlinear(context,err,error,*999)
  
  CALL Problem_SolverNonlinearMonitor(context,iterationNumber,residualNorm,err,error,*999)
  
  RETURN
999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error evaluating nonlinear residual.",err,error,*997)
997 RETURN    

END SUBROUTINE Problem_SolverNonlinearMonitorPETSC

!
!================================================================================================================================
!

!>Called from the PETSc TAO solvers to monitor an optimiser solver
SUBROUTINE Problem_SolverOptimiserMonitorPETSC(tao,context,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemRoutines
  USE SolverAccessRoutines
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The PETSc tao type
  TYPE(SOLVER_TYPE), POINTER :: context !<The passed through context (solver)
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(VARYING_STRING) :: error

  IF(.NOT.ASSOCIATED(context)) CALL FlagError("Solver context is not associated.",err,error,*999)
  CALL Solver_AssertIsOptimiser(context,err,error,*999)

  CALL Problem_SolverOptimiserMonitor(context,err,error,*999)
  
  RETURN

999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error monitoring optimiser.",err,error,*997)
997 RETURN    

END SUBROUTINE Problem_SolverOptimiserMonitorPETSC
