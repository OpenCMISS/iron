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
  USE BioelectricRoutines
  USE CellMLAccessRoutines
  USE ClassicalFieldRoutines
  USE ComputationAccessRoutines
  USE ContextAccessRoutines
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DataPointAccessRoutines
  USE DataProjectionAccessRoutines
  USE DistributedMatrixVector
  USE ElasticityRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsSetRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FIELD_IO_ROUTINES
  USE FiniteElasticityRoutines
  USE FittingRoutines
  USE FluidMechanicsRoutines
  USE InputOutput
  USE InterfaceConditionRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MultiPhysicsRoutines
  USE ProblemAccessRoutines
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesRoutines
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
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_CellMLEquationsCreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_FINISH_ACTION
    !Finish problem specific startup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
      
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
    TYPE(ProblemSetupType) :: problemSetupInfo

    ENTERS("Problem_CellMLEquationsCreateStart",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    
    !Initialise the problem setup information
    CALL Problem_SetupInitialise(problemSetupInfo,err,error,*999)
    problemSetupInfo%setupType=PROBLEM_SETUP_CELLML_EQUATIONS_TYPE
    problemSetupInfo%actionType=PROBLEM_SETUP_START_ACTION
    !Start the problem specific control setup
    CALL Problem_Setup(problem,problemSetupInfo,err,error,*999)
    !Finalise the problem setup information
    CALL Problem_SetupFinalise(problemSetupInfo,err,error,*999)
       
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
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer to the CellML equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: timeDependence,solverOutputType
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_CellMLEquationsSolve",err,error,*999)

    CALL CellMLEquations_AssertIsFinished(cellMLEquations,err,error,*999)
 
    NULLIFY(solver)
    CALL CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    IF(solverOutputType>=SOLVER_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"CellML equations solve: ",solver%label,err,error,*999)
    ENDIF

    CALL CellMLEquations_TimeDependenceTypeGet(cellMLEquations,timeDependence,err,error,*999)
    SELECT CASE(timeDependence)
    CASE(CELLML_EQUATIONS_STATIC)
      !Do nothing
    CASE(CELLML_EQUATIONS_QUASISTATIC,CELLML_EQUATIONS_DYNAMIC)
      NULLIFY(controlLoop)
      CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      CALL CellMLEquations_TimeSet(cellMLEquations,currentTime,err,error,*999)
    CASE DEFAULT
      localError="The CellML equations time dependence type of "//TRIM(NumberToVString(timeDependence,"*",err,error))// &
        & " is invalid."
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
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML to evaluate
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
    TYPE(CellMLModelType), POINTER :: model
    TYPE(FieldType), POINTER :: intermediateField,modelsField,parametersField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    
    ENTERS("Problem_SolverDAECellMLRHSEvaluate",err,error,*999)

    NULLIFY(modelsField)
    CALL CellML_ModelsFieldGet(cellML,modelsField,err,error,*999)

    CALL CellML_MaximumNumberOfStateGet(cellML,maxNumberOfStates,err,error,*999)
    CALL CellML_MaximumNumberOfIntermediateGet(cellML,maxNumberOfIntermediates,err,error,*999)
    CALL CellML_MaximumNumberOfParametersGet(cellML,maxNumberOfParameters,err,error,*999)
    
    !Make sure CellML fields have been updated to the current value of any mapped fields
    NULLIFY(modelsVariable)
    CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
    CALL Field_DOFOrderTypeGet(modelsField,FIELD_U_VARIABLE_TYPE,dofOrderType,err,error,*999)
    CALL Field_ParameterSetDataGet(modelsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
    modelIdx=modelsData(dofIdx)
    NULLIFY(model)
    CALL CellML_CellMLModelGet(cellML,modelIdx,model,err,error,*999)
    IF(dofOrderType==FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,parameterDataOffset,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,intermediateDataOffset,err,error,*999)
    ELSE
      parameterDataOffset=maxNumberOfParameters
      intermediateDataOffset=maxNumberOfIntermediates
    ENDIF
    !Get the parameters information if this environment has any.
    NULLIFY(parametersField)
    NULLIFY(parameterData)
    CALL CellML_ParametersFieldExists(cellML,parametersField,err,error,*999)
    IF(ASSOCIATED(parametersField)) THEN
      CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parameterData, &
        & err,error,*999)
    ENDIF
    !Get the intermediate information if this environment has any.
    NULLIFY(intermediateField)
    NULLIFY(intermediateData)
    CALL CellML_IntermediateFieldExists(cellML,intermediateField,err,error,*999)
    IF(ASSOCIATED(intermediateField)) THEN
      CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData, &
        & err,error,*999)
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

  !>Sets up a problem control loop
  SUBROUTINE Problem_RootControlLoopSetup(rootControlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: rootControlLoop !<A pointer to the root control loop to setup.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("Problem_RootControlLoopSetup",err,error,*999)

    CALL ControlLoop_AssertIsFinished(rootControlLoop,err,error,*999)
    CALL ControlLoop_AssertIsRootLoop(rootControlLoop,err,error,*999)

    !Setup any variables involved in the control loop (and any subloops)   
    CALL ControlLoop_FieldVariablesCalculate(rootControlLoop,err,error,*999)
    
    !Setup the times of any time control loops.
    CALL ControlLoop_TimesSetup(rootControlLoop,0.0_DP,err,error,*999)
    
    !Initialise any solvers in the control loop
    CALL Problem_ControlLoopSolversSetup(rootControlLoop,err,error,*999)
       
    EXITS("Problem_RootControlLoopSetup")
    RETURN
999 ERRORSEXITS("Problem_RootControlLoopSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_RootControlLoopSetup

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
    INTEGER(INTG) :: controlOutputType,currentIteration,inputFrequency,iterationIdx,iterationIncrement,loopIdx,loopType, &
      & numberOfSolvers,numberOfSubLoops,outputFrequency,solverIdx,startIteration,stopIteration
    TYPE(ControlLoopType), POINTER :: controlLoop2
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop
    TYPE(ControlLoopSimpleType), POINTER :: simpleLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_ControlLoopSolve",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)

    !Solve this control loop
    CALL ControlLoop_OutputTypeGet(controlLoop,controlOutputType,err,error,*999)
    IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
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
    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    SELECT CASE(loopType)
    CASE(CONTROL_SIMPLE_TYPE)
      NULLIFY(simpleLoop)
      CALL ControlLoop_SimpleLoopGet(controlLoop,simpleLoop,err,error,*999)
      IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Simple control loop: ",err,error,*999)
      ENDIF
      IF(diagnostics1) THEN
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Simple control loop: ",err,error,*999)
      ENDIF
      CALL Problem_PreLoop(controlLoop,err,error,*999)
      IF(numberOfSubLoops==0) THEN
        !If there are no sub loops then solve.
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
        DO solverIdx=1,numberOfSolvers
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
          
          CALL Problem_SolverSolve(solver,err,error,*999)
          
        ENDDO !solverIdx
      ELSE
        !If there are sub loops the recursively solve those control loops
        DO loopIdx=1,numberOfSubLoops
          NULLIFY(controlLoop2)
          CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
          CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
        ENDDO !loopIdx
      ENDIF
      CALL Problem_PostLoop(controlLoop,err,error,*999)
    CASE(CONTROL_FIXED_LOOP_TYPE)
      NULLIFY(fixedLoop)
      CALL ControlLoop_FixedLoopGet(controlLoop,fixedLoop,err,error,*999)
      CALL ControlLoop_CurrentFixedInformationGet(controlLoop,currentIteration,startIteration,stopIteration,iterationIncrement, &
        & outputFrequency,inputFrequency,err,error,*999)
      DO iterationIdx=startIteration,stopIteration,iterationIncrement
        IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Fixed control loop iteration: ",iterationIdx,err,error,*999)
        ENDIF
        IF(diagnostics1) THEN
          CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Fixed control loop iteration: ",iterationIdx,err,error,*999)
        ENDIF
        fixedLoop%iterationNumber=iterationIdx
        CALL Problem_PreLoop(controlLoop,err,error,*999)
        IF(numberOfSubLoops==0) THEN
          !If there are no sub loops then solve
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
          DO solverIdx=1,numberOfSolvers
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        CALL Problem_PostLoop(controlLoop,err,error,*999)
      ENDDO !iterationIdx
    CASE(CONTROL_TIME_LOOP_TYPE)
      !Get the time loop
      NULLIFY(timeLoop)
      CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
      
      !Precompute the number of iterations from total time span and time increment if it was not specified explicitely 
      IF(timeLoop%numberOfIterations==0) THEN
        timeLoop%numberOfIterations=CEILING((timeLoop%stopTime-timeLoop%startTime)/timeLoop%timeIncrement)
        !If number of iterations was specified but does not match timeIncrement, e.g. timeIncrement is still at the default
        !value, compute correct timeIncrement.
      ELSE IF(CEILING((timeLoop%stopTime-timeLoop%startTime)/timeLoop%timeIncrement) /= timeLoop%numberOfIterations) THEN
        timeLoop%timeIncrement = (timeLoop%stopTime-timeLoop%startTime)/timeLoop%numberOfIterations
      ENDIF
      
      !Initialise the current time and iteration number. Solvers have been initialised in Problem_ControlLoopSolversSetup for the
      !start time.
      CALL ControlLoop_CurrentTimeSetup(controlLoop,err,error,*999)
      timeLoop%iterationNumber=0
      
      DO WHILE(timeLoop%iterationNumber<timeLoop%numberOfIterations)
        !Increment loop counter and time
        timeLoop%iterationNumber=timeLoop%iterationNumber+1
        timeLoop%globalIterationNumber=timeLoop%globalIterationNumber+1
        timeLoop%currentTime=timeLoop%currentTime+timeLoop%timeIncrement
        IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
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
        CALL Problem_PreLoop(controlLoop,err,error,*999)
        IF(numberOfSubLoops==0) THEN
          !If there are no sub loops then solve.
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
          DO solverIdx=1,numberOfSolvers
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        !Perform any post loop actions.
        CALL Problem_PostLoop(controlLoop,err,error,*999)
      ENDDO !time loop
    CASE(CONTROL_WHILE_LOOP_TYPE)
      NULLIFY(whileLoop)
      CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,ERR,ERROR,*999)
      whileLoop%iterationNumber=0
      whileLoop%continueLoop=.TRUE.
      DO WHILE(whileLoop%continueLoop.AND.whileLoop%iterationNumber<whileLoop%maximumNumberOfIterations)
        whileLoop%iterationNumber=whileLoop%iterationNumber+1
        IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
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
        CALL Problem_PreLoop(controlLoop,err,error,*999)
        IF(numberOfSubLoops==0) THEN
          !If there are no sub loops then solve
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
          DO solverIdx=1,numberOfSolvers
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsExists(solver,solverEquations,err,error,*999)
            IF(ASSOCIATED(solverEquations)) CALL Problem_SolverLoadIncrementApply(solverEquations,1,1,err,error,*999)
            
            CALL Problem_SolverSolve(solver,err,error,*999)
            
          ENDDO !solverIdx
        ELSE
          !If there are sub loops the recursively solve those control loops
          DO loopIdx=1,numberOfSubLoops
            NULLIFY(controlLoop2)
            CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
            CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
          ENDDO !loopIdx
        ENDIF
        CALL Problem_PostLoop(controlLoop,err,error,*999)
      ENDDO !while loop
    CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
!!TODO: Should probably do something like the time loops for nested load increment loops and any upstream increment.
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
          IF(controlOutputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
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
          CALL Problem_PreLoop(controlLoop,err,error,*999)
          IF(numberOfSubLoops==0) THEN
            !If there are no sub loops then solve
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
            DO solverIdx=1,numberOfSolvers
              NULLIFY(solver)
              CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
              !Apply incremented boundary conditions here =>
              NULLIFY(solverEquations)
              CALL Solver_SolverEquationsExists(solver,solverEquations,err,error,*999)
              IF(ASSOCIATED(solverEquations)) &
                & CALL Problem_SolverLoadIncrementApply(solverEquations,loadIncrementLoop%iterationNumber, &
                & loadIncrementLoop%maximumNumberOfIterations,err,error,*999)
              
              CALL Problem_SolverSolve(solver,err,error,*999)
                
            ENDDO !solverIdx
          ELSE
            !If there are sub loops the recursively solve those control loops
            DO loopIdx=1,numberOfSubLoops
              NULLIFY(controlLoop2)
              CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
              CALL Problem_ControlLoopSolve(controlLoop2,err,error,*999)
            ENDDO !loopIdx
          ENDIF
          CALL Problem_PostLoop(controlLoop,err,error,*999)
        ENDDO !while loop
      ENDIF
    CASE DEFAULT
      localError="The control loop loop type of "//TRIM(NumberToVString(loopType,"*",err,error))//" is invalid."
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

  !>Sets up the solvers in a problem control loop
  RECURSIVE SUBROUTINE Problem_ControlLoopSolversSetup(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to setup the solvers for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: iterationIdx,loopIdx,numberOfSolvers,numberOfSubLoops,solverIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_ControlLoopSolversSetup",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)

    !Solve this control loop
    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Control loop: ",controlLoop%label,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Control loop level = ",controlLoop%controlLoopLevel,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sub loop index     = ",controlLoop%subLoopIndex,err,error,*999)
    ENDIF
    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    IF(numberOfSubLoops==0) THEN
      !If there are no sub loops then setup the solvers.
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
      DO solverIdx=1,numberOfSolvers
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
        
        CALL Problem_SolverSetup(solver,err,error,*999)
        
      ENDDO !solverIdx
    ELSE
      !If there are sub loops then recursively setup the solvers in those control loops
      DO loopIdx=1,numberOfSubLoops
        NULLIFY(controlLoop2)
        CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
        CALL Problem_ControlLoopSolversSetup(controlLoop2,err,error,*999)
      ENDDO !loopIdx
    ENDIF
       
    EXITS("Problem_ControlLoopSolversSetup")
    RETURN
999 ERRORSEXITS("Problem_ControlLoopSolversSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_ControlLoopSolversSetup

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
    problem%specificationLength=0
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
      CALL Bioelectric_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_FITTING_CLASS)
      CALL Fitting_ProblemSetup(problem,problemSetupInfo,err,error,*999)
    CASE(PROBLEM_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_ProblemSetup(problem,problemSetupInfo,err,error,*999)
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,solverMatrixIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(SolverType), POINTER :: cellMLSolver,linkingSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
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
      DO solverMatrixIdx=1,solverMatrices%numberOfMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,solverMatrix%solverVector,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
    !Check if the nonlinear solver is linked to a dynamic solver 
    linkingSolver=>solver%linkingSolver
    IF(ASSOCIATED(linkingSolver)) THEN
      IF(linkingSolver%solveType/=SOLVER_DYNAMIC_TYPE) &
        & CALL FlagError("Solver equations linking solver mapping is not dynamic.",err,error,*999)
      !Update the field values from the dynamic factor * current solver values AND add in mean predicted displacements/
      CALL Solver_VariablesDynamicNonlinearUpdate(solver,err,error,*999)
      !check for a linked CellML solver 
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      SELECT CASE(nonlinearSolver%nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        cellMLSolver=>newtonSolver%cellMLEvaluatorSolver
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        cellMLSolver=>quasiNewtonSolver%cellMLEvaluatorSolver
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//"."
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
      CALL Solver_VariablesFieldUpdate(solver,err,error,*999)
      !check for a linked CellML solver 
!!TODO: This should be generalised for nonlinear solvers in general and not just Newton solvers.
      newtonSolver=>solver%nonlinearSolver%newtonSolver
      IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Nonlinear solver Newton solver is not associated.",err,error,*999)
      cellMLSolver=>newtonSolver%cellMLEvaluatorSolver
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
      !  CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsLinearity,equationsSetIdx,interfaceConditionIdx,linkingSolveType,nonlinearSolveType,numberOfEquationsSets,numberOfInterfaceConditions,numberOfMatrices,outputType,solverMatrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(SolverType), POINTER :: cellMLSolver,linkingSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix    
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverResidualEvaluate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)

    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    IF(outputType>=SOLVER_MATRIX_OUTPUT) THEN
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver vector values:",err,error,*999)
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
      DO solverMatrixIdx=1,numberOfMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,solverMatrix%solverVector,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
    !Check if the nonlinear solver is linked to a dynamic solver
    NULLIFY(linkingSolver)
    CALL Solver_LinkingSolverExists(solver,linkingSolver,err,error,*999)
    IF(ASSOCIATED(linkingSolver)) THEN
      CALL Solver_SolverTypeGet(linkingSolver,linkingSolveType,err,error,*999)
      IF(linkingSolveType/=SOLVER_DYNAMIC_TYPE) &
        & CALL FlagError("Solver equations linking solver mapping is not dynamic.",err,error,*999)
      !Update the field values from the dynamic factor*current solver values AND add in predicted displacements
      CALL Solver_VariablesDynamicNonlinearUpdate(solver,err,error,*999)
      !Caculate the strain field for an CellML evaluator solver
      CALL Problem_PreResidualEvaluate(solver,err,error,*999)
      !Evaluate residual
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver residual evaluate: ",solver%label,err,error,*999)
      ENDIF
      !check for a linked CellML solver
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      NULLIFY(cellMLSolver)
      CALL SolverNonlinear_SolverTypeGet(nonlinearSolver,nonlinearSolveType,err,error,*999)
      SELECT CASE(nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        CALL SolverNonlinearNewton_LinkedCellMLSolverExists(newtonSolver,cellMLSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        CALL SolverNonlinearQuasiNewton_LinkedCellMLSolverExists(quasiNewtonSolver,cellMLSolver,err,error,*999)
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolveType,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Calculate the residual for each element (M, C, K and g)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeGet(equations,equationsLinearity,err,error,*999)
        SELECT CASE(equationsLinearity)
        CASE(EQUATIONS_LINEAR)
          !Assemble the equations for linear equations
          CALL EquationsSet_Assemble(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR)
          !Evaluate the residual for nonlinear equations
          CALL EquationsSet_ResidualEvaluate(equationsSet,err,error,*999)
        CASE DEFAULT
          localError="The equations linearity of "//TRIM(NumberToVString(equationsLinearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
      !Assemble the final solver residual.
      CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,err,error,*999)
    ELSE
      !Perform as normal nonlinear solver
      !Copy the current solution vector to the dependent field
      CALL Solver_VariablesFieldUpdate(solver,err,error,*999)
      !Caculate the strain field for an CellML evaluator solver
      CALL Problem_PreResidualEvaluate(solver,err,error,*999)
      !Evaluate residual
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solver residual evaluate: ",solver%label,err,error,*999)
      ENDIF
      !check for a linked CellML solver
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      NULLIFY(cellMLSolver)
      CALL SolverNonlinear_SolverTypeGet(nonlinearSolver,nonlinearSolveType,err,error,*999)
      SELECT CASE(nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(newtonSolver)
        CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
        CALL SolverNonlinearNewton_LinkedCellMLSolverExists(newtonSolver,cellMLSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        NULLIFY(quasiNewtonSolver)
        CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
        CALL SolverNonlinearQuasiNewton_LinkedCellMLSolverExists(quasiNewtonSolver,cellMLSolver,err,error,*999)
      CASE DEFAULT
        localError="Linked CellML solver is not implemented for nonlinear solver type " &
          & //TRIM(NumberToVString(nonlinearSolveType,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(ASSOCIATED(cellMLSolver)) CALL Solver_Solve(cellMLSolver,err,error,*999)
      !Make sure the equations sets are up to date
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeGet(equations,equationsLinearity,err,error,*999)
        SELECT CASE(equationsLinearity)
        CASE(EQUATIONS_LINEAR)
          !Assemble the equations for linear equations
          CALL EquationsSet_Assemble(equationsSet,err,error,*999)
        CASE(EQUATIONS_NONLINEAR)
          !Evaluate the residual for nonlinear equations
          CALL EquationsSet_ResidualEvaluate(equationsSet,err,error,*999)
        CASE DEFAULT
          localError="The equations linearity of "//TRIM(NumberToVString(equationsLinearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx
      !Update interface matrices
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition, &
          & err,error,*999)
        interfaceCondition=>solverMapping%interfaceConditions(interfaceConditionIdx)%ptr
        !Assemble the interface condition for the Jacobian LHS
        CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
      ENDDO !interfaceConditionIdx
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
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to pre-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
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
            CASE(EQUATIONS_SET_FITTING_CLASS)
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
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to post-evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
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
            CASE(EQUATIONS_SET_FITTING_CLASS)
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
    TYPE(ControlLoopType), POINTER :: rootControlLoop
    
    ENTERS("Problem_Solve",err,error,*999)

    CALL Problem_AssertIsFinished(problem,err,error,*999)

    NULLIFY(rootControlLoop)
    CALL Problem_ControlLoopRootGet(problem,rootControlLoop,err,error,*999)
    CALL Problem_RootControlLoopSetup(rootControlLoop,err,error,*999)
    CALL Problem_ControlLoopSolve(rootControlLoop,err,error,*999)
       
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(SolverMappingType), POINTER :: solverMapping
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
      CALL EquationsSet_LoadIncrementApply(equationsSet,solverEquations%boundaryConditions,iterationNumber, &
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
  SUBROUTINE Problem_PreLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Problem_PreLoop",err,error,*999)

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
      CALL Elasticity_PreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_BIOELECTRICS_CLASS)
      !do nothing
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_PreLoop(controlLoop,err,error,*999)
    CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
      !do nothing
    CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
      !do nothing
    CASE(PROBLEM_FITTING_CLASS)
      !do nothing
    CASE(PROBLEM_MODAL_CLASS)
      !do nothing
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_PreLoop(controlLoop,err,error,*999)
    CASE DEFAULT
      localError="Problem class "//TRIM(NumberToVString(problem%specification(1),"*",err,error))//" &
        & is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Problem_PreLoop")
    RETURN
999 ERRORSEXITS("Problem_PreLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_PreLoop

  !
  !================================================================================================================================
  !

  !>Executes after each loop of a control loop, ie after each time step for a time loop
  SUBROUTINE Problem_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Problem_PostLoop",err,error,*999)

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
     CALL Elasticity_PostLoop(controlLoop,err,error,*999)
   CASE(PROBLEM_BIOELECTRICS_CLASS)
     CALL Bioelectric_PostLoop(controlLoop,err,error,*999)
   CASE(PROBLEM_FLUID_MECHANICS_CLASS)
     CALL FluidMechanics_PostLoop(controlLoop,err,error,*999)
   CASE(PROBLEM_ELECTROMAGNETICS_CLASS)
     !Do nothing
   CASE(PROBLEM_CLASSICAL_FIELD_CLASS)
     CALL ClassicalField_PostLoop(controlLoop,err,error,*999)        
   CASE(PROBLEM_FITTING_CLASS)
     !Do nothing
   CASE(PROBLEM_MODAL_CLASS)
     !Do nothing
   CASE(PROBLEM_MULTI_PHYSICS_CLASS)
     CALL MultiPhysics_PostLoop(controlLoop,err,error,*999)
   CASE DEFAULT
     localError="The first problem specification of "// &
       & TRIM(NumberToVString(problem%specification(1),"*",err,error))//" is not valid."
     CALL FlagError(localError,err,error,*999)
   END SELECT
    
    EXITS("Problem_PostLoop")
    RETURN
999 ERRORSEXITS("Problem_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_PostLoop

  !
  !================================================================================================================================
  !

  !>Executes pre solver routines for a problem.
  SUBROUTINE Problem_SolverPreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
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
      CALL Bioelectric_PreSolve(solver,err,error,*999)
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
      CALL MultiPhysics_PreSolve(solver,err,error,*999)
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
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
      CALL Bioelectric_PostSolve(solver,err,error,*999)
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
      CALL MultiPhysics_PostSolve(solver,err,error,*999)
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets,numberOfInterfaceConditions
    REAL(DP) :: currentTime,timeIncrement
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    
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
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Set the solver time
    CALL Solver_DynamicTimesSet(solver,currentTime,timeIncrement,err,error,*999)
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsLinearity,equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets,numberOfInterfaceConditions
    REAL(DP) :: currentTime,timeIncrement
    LOGICAL :: restart,initialised
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverEquationsDynamicNonlinearSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    CALL SolverDynamic_RestartGet(dynamicSolver,restart,err,error,*999)
    CALL SolverDynamic_SolverInitialisedGet(dynamicSolver,initialised,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      IF(restart.OR..NOT.initialised) THEN
        !If we need to restart or we haven't initialised yet, make sure the equations sets are up to date
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeGet(equations,equationsLinearity,err,error,*999)
        SELECT CASE(equationsLinearity)
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
            & TRIM(NumberToVString(equationsLinearity,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Set the solver time
    CALL Solver_DynamicTimesSet(solver,currentTime,timeIncrement,err,error,*999)
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
    ENDDO !equationsSetIdx
    
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets,numberOfInterfaceConditions
    REAL(DP) :: currentTime,timeIncrement
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
     
    ENTERS("Problem_SolverEquationsQuasistaticLinearSolve",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    !Make sure the equations sets are up to date
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the current times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)    
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Solve for the current time
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets,numberOfInterfaceConditions
    REAL(DP) :: currentTime,timeIncrement
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Problem_SolverEquationsQuasistaticNonlinearSolve",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    !Get current control loop times
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    !Make sure the equations sets are up to date
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Set the equations set times
      CALL EquationsSet_TimesSet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Solve for the next time i.e., current time + time increment
    CALL Solver_Solve(solver,err,error,*999)
    !Update the field variables with residuals or backsubstitute 
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets,numberOfInterfaceConditions
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    
    ENTERS("Problem_SolverEquationsStaticLinearSolve",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif    
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Make sure the equations sets are up to date
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !CALL EQUATIONS_SET_FIXED_CONDITIONS_APPLY(equationsSet,err,error,*999)
      !Assemble the equations for linear problems
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    DO interfaceConditionIdx=1,numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Solve
    CALL Solver_Solve(solver,err,error,*999)
    !Back-substitute to find flux values for linear problems
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
    ENDDO !equationsSetIdx
    
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,interfaceConditionIdx,numberOfEquationsSets
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverEquationsStaticNonlinearSolve",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
#endif
    
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Apply boundary conditition
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      !Assemble the equations set
      CALL EquationsSet_Assemble(equationsSet,err,error,*999)
    ENDDO !equationsSetIdx
    !Make sure the interface matrices are up to date
    DO interfaceConditionIdx=1,solverMapping%numberOfInterfaceConditions
      NULLIFY(interfaceCondition)
      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
      CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
    ENDDO !interfaceConditionIdx
    !Solve
    CALL Solver_Solve(solver,err,error,*999)
    !Update the rhs field variables with residuals or backsubstitute
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    DO equationsSetIdx=1,numberOfEquationsSets
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
      CALL EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*999)
    ENDDO !equationsSetIdx
    
    EXITS("Problem_SolverEquationsStaticNonlinearSolve")
    RETURN
999 ERRORSEXITS("Problem_SolverEquationsStaticNonlinearSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverEquationsStaticNonlinearSolve

  !
  !================================================================================================================================
  !


  !>Sets up a solver for a problem.
  SUBROUTINE Problem_SolverSetup(solver,err,error,*)

   !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to setup
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsLinearity,equationsSetIdx,inputIterationNumber,interfaceConditionIdx,iterationNumber, &
      & numberOfEquationsSets,numberOfInterfaceConditions,outputIterationNumber,solverDegree,solverEquationsLinearity, &
      & solverEquationsTimeDependence,solverOrder,solveType
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    LOGICAL :: initSolver,nonlinear,setup,setupFinished,solverInitialised
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Problem_SolverSetup",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif
    
    CALL Solver_SolverSetupGet(solver,setupFinished,err,error,*999)
    IF(.NOT.setupFinished) THEN
      !Setup the solver
      CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
      SELECT CASE(solveType)
      CASE(SOLVER_LINEAR_TYPE)
        !Do nothing
      CASE(SOLVER_NONLINEAR_TYPE)
        !Do nothing
      CASE(SOLVER_DYNAMIC_TYPE)
        NULLIFY(dynamicSolver)
        CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
        CALL SolverDynamic_DegreeGet(dynamicSolver,solverDegree,err,error,*999)
        CALL SolverDynamic_OrderGet(dynamicSolver,solverOrder,err,error,*999)
        CALL SolverDynamic_SolverInitialisedGet(dynamicSolver,solverInitialised,err,error,*999)
        CALL SolverDynamic_DegreeGet(dynamicSolver,solverDegree,err,error,*999)
        CALL SolverDynamic_OrderGet(dynamicSolver,solverOrder,err,error,*999)
        !Get the solver equations linearity and time dependence.
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeGet(solverEquations,solverEquationsLinearity,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeGet(solverEquations,solverEquationsTimeDependence,err,error,*999)
        initSolver=(.NOT.solverInitialised.AND. &
          & ((solverOrder==SOLVER_DYNAMIC_FIRST_ORDER.AND.solverDegree>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
          & (solverOrder==SOLVER_DYNAMIC_SECOND_ORDER.AND.solverDegree>SOLVER_DYNAMIC_SECOND_DEGREE)))
        nonlinear=(solverEquationsLinearity==SOLVER_EQUATIONS_NONLINEAR)
        setup=(initSolver.OR.nonlinear)
        IF(setup) THEN
          !Need to setup solvers as we have a residual to evaluate at the start time or the dynamic solver needs initialisation
          !Find the start time
          NULLIFY(controlLoop)
          CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
          CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,startTime,stopTime,timeIncrement,iterationNumber, &
            & outputIterationNumber,inputIterationNumber,err,error,*999)
          CALL Solver_DynamicTimesSet(solver,startTime,timeIncrement,err,error,*999)
          !Get solver matrices
          NULLIFY(solverMatrices)
          CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
          NULLIFY(solverMatrix)
          CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
          NULLIFY(solverVector)
          CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
          !Nullify the solver vector so that alpha is zero.
          CALL DistributedVector_AllValuesSet(solverVector,0.0_DP,err,error,*999)
          !Solve the solver equations
          CALL Problem_SolverEquationsSolve(solverEquations,err,error,*999)
        ENDIF !setup  
      CASE(SOLVER_DAE_TYPE)
        !Do nothing
      CASE(SOLVER_EIGENPROBLEM_TYPE)
        !Do nothing
      CASE(SOLVER_OPTIMISER_TYPE)
        !Do nothing
      CASE(SOLVER_CELLML_EVALUATOR_TYPE)
        !Do nothing
      CASE(SOLVER_STATE_ITERATION_TYPE)
        !Do nothing
      CASE(SOLVER_GEOMETRIC_TRANSFORMATION_TYPE)
        !Do nothing
      CASE DEFAULT
        localError="The solver type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
          & " is invalid or not implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finish the setup
      CALL Solver_SolverSetupSet(solver,.TRUE.,err,error,*999)      
    ENDIF
     
    EXITS("Problem_SolverSetup")
    RETURN
999 ERRORSEXITS("Problem_SolverSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverSetup

  !
  !================================================================================================================================
  !


  !>Solves a solver for a problem.
  SUBROUTINE Problem_SolverSolve(solver,err,error,*)

   !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to solve
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
           
    CALL Problem_SolverPreSolve(solver,err,error,*999)
          
    IF(ASSOCIATED(solver%solverEquations)) THEN
      !A solver with solver equations.
      CALL Problem_SolverEquationsSolve(solver%solverEquations,err,error,*999)
    ELSE
      !Check for other equations.
      IF(ASSOCIATED(solver%cellMLEquations)) THEN
        !A solver with CellML equations.
        CALL Problem_CellMLEquationsSolve(solver%cellMLEquations,err,error,*999)
      ELSE IF(solver%solveType==SOLVER_GEOMETRIC_TRANSFORMATION_TYPE) THEN
        CALL Problem_SolverGeometricTransformationSolve(solver%geometricTransformationSolver,err,error,*999)
      ELSE
        !Do nothing now. 
        !CALL FlagError("Solver does not have any equations associated.",err,error,*999)
      ENDIF
    ENDIF

    CALL Problem_SolverPostSolve(solver,err,error,*999)
         
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
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer to the solver equations to get the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(SolverMappingType), POINTER :: solverMapping
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
    TYPE(SolverType), POINTER :: solver
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
    CALL FieldVariable_ParameterSetExists(fieldVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,boundaryConditionsParameterSet, &
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to monitor
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The number of iterations
    REAL(DP), INTENT(IN) :: residualNorm !<The residual norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: interfaceConditionIdx
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
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
        IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
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
                CALL InterfaceCondition_Assemble(interfaceCondition,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !interfaceConditionIdx
        ENDIF !Reproject
        !Output meshes at iterations
        IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
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
    IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      NULLIFY(nonlinearSolver)
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
      CALL SolverNonlinear_Monitor(nonlinearSolver,iterationNumber,residualNorm,err,error,*999)
    ELSE
      localError="Invalid solve type. The solve type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
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
    TYPE(SolverType), POINTER :: solver !<A pointer to solver to output the fields for
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<Iteration number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component,coupledElementNumber,coupledLineFaceNumber,coupledMeshFaceLineNumber,coupledMeshIdx,dataPointIdx, &
      & exitTag,globalDataPointNumber,interfaceConditionIdx,interfaceElementNumber,iUnit,loadStep,numberOfCoupledElements, &
      & numberOfCoupledMeshes,numberOfDataDimensions,numberOfInterfaceConditions,numberOfProjectedData,numberOfReducedXi, &
      & pSpecification(4),solveCall
    REAL(DP) :: position(3),reducedXi(3)
    CHARACTER(LEN=100) :: filenameOutput,groupname
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DataPointsType), POINTER :: interfaceDatapoints
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DecompositionType), POINTER :: dependentDecomposition,lagrangeDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: lagrangeDecompositionDataPoints
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData
    TYPE(DecompositionElementsType), POINTER :: dependentDecompositionElements,lagrangeDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: dependentDecompositionTopology,lagrangeDecompositionTopology
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceDependentType), POINTER :: interfaceDependent
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfacePointConnectivityType), POINTER :: pointConnectivity
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(FieldType), POINTER :: coupledMeshDependentField,lagrangeField
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldVariableType), POINTER :: coupledMeshDependentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: directory,localError

    ENTERS("Problem_SolverNewtonFieldsOutput",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      SELECT CASE(pSpecification(2))
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
        !   maxSolveCalls=100
        !   coupledMeshIdx=1
        !   loadStep=1
        !   firstIterationNumber=0
        !   DO solveCall=1,maxSolveCalls
        !     fileToCheck=directory// &
        !       & "mesh"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))// &
        !       & "_solveCall"//TRIM(NumberToVString(solveCall,"*",err,error))// &
        !       & "_load"//TRIM(NumberToVString(loadStep,"*",err,error))// &
        !       & "_iter"//TRIM(NumberToVString(firstIterationNumber,"*",err,error))//".part0.exnode"
        !     INQUIRE(FILE=CHAR(fileToCheck),EXIST=fileExists)
        !     IF(.NOT.fileExists) THEN
        !       EXIT
        !     ENDIF
        !   ENDDO
        
        !   loadStep=solver%solvers%controlLoop%loadIncrementLoop%iterationNumber
        
        !   IF((iterationNumber > 0).OR.(loadStep > 1))THEN
        !     solveCall = solveCall - 1
        !   ENDIF
        
        !   WRITE(*,'(1X,''SolveCall: '',I4)') solveCall
        !   WRITE(*,'(1X,''  LoadStep: '',I4)') loadStep
        !   WRITE(*,'(1X,''    Iteration: '',I4)') iterationNumber
        
        !   DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        !     region=>solverMapping%equationsSets(equationsSetIdx)%ptr%REGION
        !     IF(ASSOCIATED(region))THEN
        !       NULLIFY(fields)
        !       fields=>region%FIELDS
        !       fileName=directory//"mesh"//TRIM(NumberToVString(equationsSetIdx,"*",err,error))// &
        !         & "_solveCall"//TRIM(NumberToVString(solveCall,"*",err,error))// &
        !         & "_load"//TRIM(NumberToVString(loadStep,"*",err,error))// &
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
        localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
      & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
      !Do nothing???
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(pSpecification(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    SELECT CASE(pSpecification(1))
    CASE(PROBLEM_ELASTICITY_CLASS)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_LINEAR_ELASTICITY_TYPE,PROBLEM_FINITE_ELASTICITY_TYPE)
        ! Pass
      CASE(PROBLEM_LINEAR_ELASTICITY_CONTACT_TYPE,PROBLEM_FINITE_ELASTICITY_CONTACT_TYPE)
        
        IF(diagnostics1) THEN
          iUnit = 300
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
          DO interfaceConditionIdx=1,numberOfInterfaceConditions
            NULLIFY(interfaceCondition)
            CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
            NULLIFY(INTERFACE)
            CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
            NULLIFY(pointsConnectivity)
            CALL Interface_PointsConnectivityExists(INTERFACE,pointsConnectivity,err,error,*999)
            IF(ASSOCIATED(pointsConnectivity)) THEN
              NULLIFY(interfaceLagrange)
              CALL InterfaceCondition_InterfaceLagrangeGet(interfaceCondition,interfaceLagrange,err,error,*999)
              NULLIFY(lagrangeField)
              CALL InterfaceLagrange_LagrangeFieldGet(interfaceLagrange,lagrangeField,err,error,*999)
              NULLIFY(lagrangeDecomposition)
              CALL Field_DecompositionGet(lagrangeField,lagrangeDecomposition,err,error,*999)
              NULLIFY(lagrangeDecompositionTopology)
              CALL Decomposition_DecompositionTopologyGet(lagrangeDecomposition,lagrangeDecompositionTopology,err,error,*999)
              NULLIFY(lagrangeDecompositionElements)
              CALL DecompositionTopology_DecompositionElementsGet(lagrangeDecompositionTopology,lagrangeDecompositionElements, &
                & err,error,*999)
              NULLIFY(lagrangeDecompositionDataPoints)
              CALL DecompositionTopology_DecompositionDataPointsGet(lagrangeDecompositionTopology,lagrangeDecompositionDataPoints, &
                & err,error,*999)
              NULLIFY(interfaceDependent)
              CALL InterfaceCondition_InterfaceDependentGet(interfaceCondition,interfaceDependent,err,error,*999)
              NULLIFY(interfaceDataPoints)
              CALL InterfacePointsConnectivity_DataPointsGet(pointsConnectivity,interfaceDataPoints,err,error,*999)
              CALL DataPoints_NumberOfDimensionsGet(interfaceDataPoints,numberOfDataDimensions,err,error,*999)
              CALL Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*999)
              DO coupledMeshIdx=1,numberOfCoupledMeshes
               filenameOutput=directory//"PointsConnectivity"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))// &
                  & "_solveCall"//TRIM(NumberToVString(solveCall,"*",err,error))// &
                  & "_load"//TRIM(NumberToVString(loadStep,"*",err,error))// &
                  & "_iter"//TRIM(NumberToVString(iterationNumber,"*",err,error))//".exdata"
                OPEN(UNIT=iUnit,FILE=filenameOutput,STATUS="UNKNOWN",ACTION="WRITE",IOSTAT=err)
                groupname="PointsConnectivity"//TRIM(NumberToVString(coupledMeshIdx,"*",err,error))
                WRITE(iUnit,'( '' Group name: '',A)') groupname
                WRITE(iUnit,'(1X,''#Fields=4'')')
                WRITE(iUnit,'(1X,''1) coordinates, coordinate, rectangular cartesian, #Components=3'')')
                WRITE(iUnit,'(1X,''  x.  Value index= 1, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  y.  Value index= 2, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  z.  Value index= 3, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''2) error, field, rectangular cartesian, #Components=3'')')
                WRITE(iUnit,'(1X,''  x.  Value index= 4, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  y.  Value index= 5, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  z.  Value index= 6, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''3) projectedCoordinate, field, rectangular cartesian, #Components=3'')')
                WRITE(iUnit,'(1X,''  x.  Value index= 7, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  y.  Value index= 8, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''  z.  Value index= 9, #Derivatives=0'')')
                WRITE(iUnit,'(1X,''4) exitTag, field, rectangular cartesian, #Components=1'')')
                WRITE(iUnit,'(1X,''  tag.  Value index= 10, #Derivatives=0'')')
                NULLIFY(equationsSet)
                CALL InterfaceDependent_EquationsSetGet(interfaceDependent,coupledMeshIdx,equationsSet,err,error,*999)
                NULLIFY(coupledMeshDependentField)
                CALL EquationsSet_DependentFieldGet(equationsSet,coupledMeshDependentField,err,error,*999)
                NULLIFY(coupledMeshDependentVariable)
                CALL Field_VariableGet(coupledMeshDependentField,FIELD_U_VARIABLE_TYPE,coupledMeshDependentVariable,err,error,*999)
                NULLIFY(dependentDecomposition)
                CALL Field_DecompositionGet(coupledMeshDependentField,dependentDecomposition,err,error,*999)
                NULLIFY(dependentDecompositionTopology)
                CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
                NULLIFY(dependentDecompositionElements)
                CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
                  & err,error,*999)
                NULLIFY(interpolationParameters)
                CALL FieldVariable_InterpolationParameterInitialise(coupledMeshDependentVariable,interpolationParameters, &
                  & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                NULLIFY(interpolatedPoint)
                CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999, &
                  & FIELD_GEOMETRIC_COMPONENTS_TYPE)
                NULLIFY(dataProjection)
                CALL DataPoints_DataProjectionIndexGet(interfaceDataPoints,coupledMeshIdx+1,dataProjection,err,error,*999)
                CALL InterfacePointsConnectivity_MaximumCoupledElementsGet(pointsConnectivity,coupledMeshIdx, &
                  & numberOfCoupledElements,err,error,*999)
                DO interfaceElementNumber=1,numberOfCoupledElements
                  NULLIFY(decompositionElementData)
                  CALL DecompositionDataPoints_ElementDataPointsGet(lagrangeDecompositionDataPoints,interfaceElementNumber, &
                    & decompositionElementData,err,error,*999)
                  CALL DecompositionElementDataPoints_NumberOfDataPointsGet(decompositionElementData,numberOfProjectedData, &
                    & err,error,*999)
                  DO dataPointIdx=1,numberOfProjectedData
                    CALL DecompositionElementDataPoints_GlobalNumberGet(decompositionElementData,dataPointIdx, &
                      & globalDataPointNumber,err,error,*999)
                    CALL DataPoints_DataPositionGet(interfaceDataPoints,globalDataPointNumber,position(1:numberOfDataDimensions), &
                      & err,error,*999)
                    WRITE(iUnit,'(1X,''Node:'',I4)') globalDataPointNumber
                    DO component=1,numberOfDataDimensions
                      WRITE(iUnit,'(1X,3E25.15)') position(component)
                    ENDDO !component
                    NULLIFY(pointConnectivity)
                    CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,globalDataPointNumber,coupledMeshIdx, &
                      & pointConnectivity,err,error,*999)
                    CALL InterfacePointConnectivity_CoupledElementNumberGet(pointConnectivity,coupledElementNumber, &
                      & err,error,*999)
                    CALL InterfacePointConnectivity_CoupledLineFaceNumberGet(pointConnectivity,coupledLineFaceNumber, &
                      & err,error,*999)
                    CALL InterfacePointConnectivity_ReducedXiGet(pointConnectivity,numberOfReducedXi,reducedXi,err,error,*999)
                    CALL DecompositionElements_ElementFaceNumberGet(dependentDecompositionElements,coupledLineFaceNumber, &
                      & coupledElementNumber,coupledMeshFaceLineNumber,err,error,*999)
                    CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,coupledMeshFaceLineNumber, &
                      & interpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
                    CALL Field_InterpolateXi(NO_PART_DERIV,reducedXi(1:numberOfReducedXi),interpolatedPoint, &
                      & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
                    DO component=1,numberOfDataDimensions
                      WRITE(iUnit,'(1X,3E25.15)') interpolatedPoint%values(component,NO_PART_DERIV) - position(component)
                    ENDDO !component                    
                    DO component=1,numberOfDataDimensions
                      WRITE(iUnit,'(1X,3E25.15)') interpolatedPoint%values(component,NO_PART_DERIV)
                    ENDDO !component
                    CALL DataProjection_ResultExitTagGet(dataProjection,globalDataPointNumber,exitTag,err,error,*999)
                    WRITE(iUnit,'(1X,I2)') exitTag
                  ENDDO !dataPointIdx
                ENDDO !interfaceElementNumber
                CALL Field_InterpolatedPointFinalise(interpolatedPoint,err,error,*999)
                CALL FieldVariable_InterpolationParameterFinalise(interpolationParameters,err,error,*999)
                OPEN(UNIT=iUnit)
              ENDDO !coupledMeshIdx
            ENDIF
          ENDDO !interfaceConditionIdx
        ENDIF
        
      CASE DEFAULT
        localError="The problem type of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_BIOELECTRICS_CLASS,PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_ELECTROMAGNETICS_CLASS, &
      & PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_FITTING_CLASS,PROBLEM_MODAL_CLASS,PROBLEM_MULTI_PHYSICS_CLASS)
      !Do nothing???
    CASE DEFAULT
      localError="The problem class of "//TRIM(NumberToVString(pSpecification(1),"*",err,error))//" is invalid."
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to monitor
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
    CALL SolverOptimiser_Monitor(optimiserSolver,err,error,*999)
    
    EXITS("Problem_SolverOptimiserMonitor")
    RETURN
999 ERRORSEXITS("Problem_SolverOptimiserMonitor",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SolverOptimiserMonitor
  
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
    INTEGER(INTG) :: problemClass,specificationIdx
    TYPE(VARYING_STRING) :: localError

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
    !Set the specification length
    problem%specificationLength=0
    DO specificationIdx=1,SIZE(problem%specification,1)
      IF(problem%specification(specificationIdx)>0) THEN
        problem%specificationLength=specificationIdx
      ENDIF
    ENDDO !specificationIdx

    EXITS("Problem_SpecificationSet")
    RETURN
999 ERRORSEXITS("Problem_SpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Problem_SpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the work group for a problem. \see OpenCMISS::Iron::cmfe_Problem_WorkGroupSet
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
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMatricesRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
 
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc snes
  TYPE(PetscVecType), INTENT(INOUT) :: X !<The PETSc x Vec
  TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr,nonlinearSolveType,numberOfMatrices
  TYPE(DistributedVectorType), POINTER :: solverVector
  TYPE(NewtonSolverType), POINTER :: newtonSolver
  TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
  TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
  TYPE(SolverEquationsType), POINTER :: solverEquations
  TYPE(SolverMatricesType), POINTER :: solverMatrices
  TYPE(SolverMatrixType), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*998)

  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*999)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
  CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
  IF(numberOfMatrices/=1) THEN
    localError="The number of solver matrices of "//TRIM(NumberToVString(numberOfMatrices,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*998)
  ENDIF
  
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
  
  CALL DistributedVector_OverrideSetOn(solverVector,x,err,error,*999)
  
  CALL Problem_SolverJacobianEvaluate(ctx,err,error,*999)
  
  CALL DistributedVector_OverrideSetOff(solverVector,err,error,*999)
  
!!TODO: move this to Problem_SolverJacobianEvaluate or elsewhere?
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)
  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*999)
  CALL SolverNonlinear_SolverTypeGet(nonlinearSolver,nonlinearSolveType,err,error,*999)
  SELECT CASE(nonlinearSolveType)
  CASE(SOLVER_NONLINEAR_NEWTON)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    newtonSolver%totalNumberOfJacobianEvaluations=newtonSolver%totalNumberOfJacobianEvaluations+1
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    quasiNewtonSolver%totalNumberOfJacobianEvaluations=quasiNewtonSolver%totalNumberOfJacobianEvaluations+1
  CASE DEFAULT
    localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolveType,"*",err,error))//" is invalid."
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
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMatricesRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE

  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec
  TYPE(PetscMatType), INTENT(INOUT) :: A !<The PETSc A Mat
  TYPE(PetscMatType), INTENT(INOUT) :: B !<The PETSc B Mat
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr,nonlinearSolveType,numberOfMatrices,sparsityType
  TYPE(NewtonSolverType), POINTER :: newtonSolver
  TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
  TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
  TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
  TYPE(QuasiNewtonLinesearchSolverType), POINTER :: quasiNewtonLinesearchSolver
  TYPE(SolverEquationsType), POINTER :: solverEquations
  TYPE(SolverMatricesType), POINTER :: solverMatrices
  TYPE(SolverMatrixType), POINTER :: solverMatrix
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
  CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
  IF(numberOfMatrices/=1) THEN
    localError="The number of solver matrices of "//TRIM(NumberToVString(numberOfMatrices,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*998)
  ENDIF
  
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
  CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
  SELECT CASE(sparsityType)
  CASE(SOLVER_SPARSE_MATRICES)
    CALL SolverNonlinear_SolverTypeGet(nonlinearSolver,nonlinearSolveType,err,error,*999)
    SELECT CASE(nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      NULLIFY(linesearchSolver)
      CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
      jacobianMatFDColoring=>linesearchSolver%jacobianMatFDColoring
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      NULLIFY(quasiNewtonLinesearchSolver)
      CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,quasiNewtonLinesearchSolver,err,error,*999)
      jacobianMatFDColoring=>quasiNewtonLinesearchSolver%jacobianMatFDColoring
    CASE DEFAULT
      localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(.NOT.ASSOCIATED(jacobianMatFDColoring)) CALL FlagError("Linesearch solver FD colouring is not associated.",err,error,*998)
    CALL Petsc_SnesComputeJacobianDefaultColor(snes,x,A,B,jacobianMatFDColoring,err,error,*999)
  CASE(SOLVER_FULL_MATRICES)
    CALL Petsc_SnesComputeJacobianDefault(snes,x,A,B,ctx,err,error,*999)
  CASE DEFAULT
    localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
    CALL FlagError(localError,err,error,*999)
  END SELECT
  IF(ctx%outputType>=SOLVER_MATRIX_OUTPUT) THEN
    CALL DistributedMatrix_OverrideSetOn(solverMatrix%matrix,A,err,error,*999)
    CALL SolverMatrices_Output(GENERAL_OUTPUT_TYPE,SOLVER_MATRICES_JACOBIAN_ONLY,solverMatrices,err,error,*998)
    CALL DistributedMatrix_OverrideSetOff(solverMatrix%matrix,err,error,*999)
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
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscTaoType), INTENT(INOUT) :: tao !<The PETSc tao type
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec type
  REAL(DP), INTENT(OUT) :: f !<On exit, the evaluated objective
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(DistributedVectorType), POINTER :: solverVector
  TYPE(SolverEquationsType), POINTER :: solverEquations
  TYPE(SolverMatricesType), POINTER :: solverMatrices
  TYPE(SolverMatrixType), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*997)
  CALL Solver_AssertIsOptimiser(ctx,err,error,*999)
  
  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*999)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)  
  IF(solverMatrices%numberOfMatrices/=1) THEN
    localError="The number of solver matrices of "//TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*997)          
  ENDIF
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)

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
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE
  
  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc snes type
  TYPE(PetscVecType), INTENT(INOUT) :: x !<The PETSc x Vec type
  TYPE(PetscVecType), INTENT(INOUT) :: f !<The PETSc f Vec type
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  INTEGER(INTG) :: dummyErr
  TYPE(DistributedVectorType), POINTER :: residualVector,solverVector
  TYPE(NewtonSolverType), POINTER :: newtonSolver
  TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
  TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
  TYPE(SolverEquationsType), POINTER :: solverEquations
  TYPE(SolverMatricesType), POINTER :: solverMatrices
  TYPE(SolverMatrixType), POINTER :: solverMatrix
  TYPE(VARYING_STRING) :: dummyError,error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*997)
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)

  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*997)
  SELECT CASE(nonLinearSolver%nonlinearSolveType)
  CASE(SOLVER_NONLINEAR_NEWTON)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*997)
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*997)
  CASE DEFAULT
    localError="The nonlinear solver type of "//TRIM(NumberToVString(nonLinearSolver%nonlinearSolveType,"*",err,error))// &
      & " is invalid."
    CALL FlagError(localError,err,error,*997)
  END SELECT
  NULLIFY(solverEquations)
  CALL Solver_SolverEquationsGet(ctx,solverEquations,err,error,*997)
  NULLIFY(solverMatrices)
  CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*997)
  IF(solverMatrices%numberOfMatrices/=1) THEN
    localError="The number of solver matrices of "// &
      & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))// &
      & " is invalid. There should be 1 solver matrix."
    CALL FlagError(localError,err,error,*997)          
  ENDIF
  NULLIFY(solverMatrix)
  CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*997)
  NULLIFY(solverVector)
  CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*997)
  NULLIFY(residualVector)
  CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualVector,err,error,*997)
  
  CALL DistributedVector_OverrideSetOn(solverVector,X,err,error,*998)
  CALL DistributedVector_OverrideSetOn(residualVector,F,err,error,*999)                
                    
  CALL Problem_SolverResidualEvaluate(ctx,err,error,*999)
  
  CALL DistributedVector_OverrideSetOff(residualVector,err,error,*999)
  CALL DistributedVector_OverrideSetOff(solverVector,err,error,*998)
  
!!TODO: move this to Problem_SolverResidualEvaluate or elsewhere?
  SELECT CASE(nonLinearSolver%nonlinearSolveType)
  CASE(SOLVER_NONLINEAR_NEWTON)
    newtonSolver%totalNumberOfFunctionEvaluations=newtonSolver%totalNumberOfFunctionEvaluations+1
  CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
    quasiNewtonSolver%totalNumberOfFunctionEvaluations=quasiNewtonSolver%totalNumberOfFunctionEvaluations+1
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
  USE SolverRoutines
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
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(PetscVecType) :: x,f,y,w,g
  TYPE(NewtonSolverType), POINTER :: newtonSolver
  TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
  TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
  TYPE(PetscSnesLinesearchType) :: lineSearch
  REAL(DP) :: energy,normalisedEnergy
  TYPE(VARYING_STRING) :: error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*999)
  CALL Solver_AssertIsNonlinear(ctx,err,error,*999)

  NULLIFY(nonlinearSolver)
  CALL Solver_NonlinearSolverGet(ctx,nonlinearSolver,err,error,*999)
  
  SELECT CASE(nonlinearSolver%nonlinearSolveType)
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
          IF(ABS(normalisedEnergy)<newtonSolver%absoluteTolerance) THEN
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
          IF(ABS(normalisedEnergy)<quasiNewtonSolver%absoluteTolerance) THEN
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
  TYPE(CellMLType), POINTER :: cellML
  TYPE(SolverType), POINTER :: solver
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
  TYPE(SolverType), POINTER :: context !<The passed through context
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
  TYPE(SolverType), POINTER :: context !<The passed through context (solver)
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
