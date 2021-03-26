!> \file
!> \author Chris Bradley
!> \brief This module handles all control loop routines.
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

!> \defgroup OpenCMISS_ControlLoop OpenCMISS::Iron::ControlLoop
!> This module handles all control loop routines.
MODULE ControlLoopRoutines

  USE BaseRoutines
  USE Constants
  USE ControlLoopAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE ProblemAccessRoutines
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesAccessRoutines  
  USE Strings
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE ControlLoop_LabelGet
    MODULE PROCEDURE ControlLoop_LabelGetC
    MODULE PROCEDURE ControlLoop_LabelGetVS
  END INTERFACE ControlLoop_LabelGet
  
  INTERFACE ControlLoop_LabelSet
    MODULE PROCEDURE ControlLoop_LabelSetC
    MODULE PROCEDURE ControlLoop_LabelSetVS
  END INTERFACE ControlLoop_LabelSet  

  PUBLIC ControlLoop_AbsoluteToleranceSet
  
  PUBLIC ControlLoop_CreateFinish,ControlLoop_CreateStart

  PUBLIC ControlLoop_Destroy

  PUBLIC ControlLoop_FieldVariablesCalculate
  
  PUBLIC ControlLoop_IterationsSet

  PUBLIC ControlLoop_LabelGet,ControlLoop_LabelSet

  PUBLIC ControlLoop_LoadOutputSet

  PUBLIC ControlLoop_MaximumIterationsSet

  PUBLIC ControlLoop_NumberOfIterationsSet

  PUBLIC ControlLoop_NumberOfSubLoopsSet

  PUBLIC ControlLoop_OutputTypeSet

  PUBLIC ControlLoop_PreviousValuesUpdate

  PUBLIC ControlLoop_RelativeToleranceSet

  PUBLIC ControlLoop_SolversDestroy

  PUBLIC ControlLoop_SolverEquationsDestroy

  PUBLIC ControlLoop_TimesGet,ControlLoop_TimesSet

  PUBLIC ControlLoop_TypeSet
  
  PUBLIC ControlLoop_TimeInputSet

  PUBLIC ControlLoop_TimeOutputSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Finish the process of creating a control loop
  RECURSIVE SUBROUTINE ControlLoop_CreateFinish(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to finish.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2
    
    ENTERS("ControlLoop_CreateFinish",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    
    !Finish the sub-loops first
    IF(controlLoop%numberOfSubLoops>0) THEN
      DO loopIdx=1,controlLoop%numberOfSubLoops
        controlLoop2=>controlLoop%subLoops(loopIdx)%ptr
        CALL ControlLoop_CreateFinish(controlLoop2,err,error,*999)
      ENDDO !loopIdx
    ENDIF
    !Finish this control loop
    controlLoop%controlLoopFinished=.TRUE.
       
    EXITS("ControlLoop_CreateFinish")
    RETURN
999 ERRORSEXITS("ControlLoop_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the process of creating a control loop for a problem.
  SUBROUTINE ControlLoop_CreateStart(problem,controlLoop,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER, INTENT(INOUT) :: problem !<A pointer to the problem to initialise the control for.
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<On exit, a pointer to the control loop. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ControlLoop_CreateStart",err,error,*998)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    
    CALL ControlLoop_Initialise(problem%controlLoop,err,error,*999)        
    problem%controlLoop%problem=>problem
    problem%controlLoop%controlLoopLevel=1
    controlLoop=>problem%controlLoop
    
    EXITS("ControlLoop_CreateStart")
    RETURN
999 NULLIFY(controlLoop)
998 ERRORSEXITS("ControlLoop_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_CreateStart

  !
  !================================================================================================================================
  !

  !>Destroy a control loop \see OpenCMISS::Iron::cmfe_ControlLoop_Destroy
  SUBROUTINE ControlLoop_Destroy(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    CALL ControlLoop_Finalise(controlLoop,err,error,*999)
       
    EXITS("ControlLoop_Destroy")
    RETURN
999 ERRORSEXITS("ControlLoop_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises a ControlLoopFieldVariableType and deallocates all memory
  SUBROUTINE ControlLoop_FieldVariableFinalise(controlLoopFieldVariable,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariableType), INTENT(INOUT) :: controlLoopFieldVariable !<The control loop field variable to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_FieldVariableFinalise",err,error,*999)

    NULLIFY(controlLoopFieldVariable%fieldVariable)
    controlLoopFieldVariable%timeDependence=0
    controlLoopFieldVariable%linearity=0
      
    EXITS("ControlLoop_FieldVariableFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariableFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a ControlLoopFieldVariableType
  SUBROUTINE ControlLoop_FieldVariableInitialise(controlLoopFieldVariable,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariableType), INTENT(INOUT) :: controlLoopFieldVariable !<The control loop field variable to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("ControlLoop_FieldVariableInitialise",err,error,*999)

    NULLIFY(controlLoopFieldVariable%fieldVariable)
    controlLoopFieldVariable%timeDependence=0
    controlLoopFieldVariable%linearity=0
      
    EXITS("ControlLoop_FieldVariableInitialise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariableInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableInitialise

  !
  !================================================================================================================================
  !

  !>Adds a field variable to the list of control loop field variables
  SUBROUTINE ControlLoop_FieldVariableAdd(controlLoopFieldVariables,variableLinearity,variableTimeDependence,fieldVariable, &
    & err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariablesType), POINTER :: controlLoopFieldVariables !<The control loop field variables to add the field variable to
    INTEGER(INTG), INTENT(IN) :: variableLinearity !<The linearity of the field variable
    INTEGER(INTG), INTENT(IN) :: variableTimeDependence !<The time dependence of the field variable
    TYPE(FieldVariableType), POINTER :: fieldVariable !<A pointer to the field variable to add.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableIdx
    LOGICAL :: found
    TYPE(ControlLoopFieldVariableType), ALLOCATABLE :: newFieldVariables(:)
    TYPE(FieldVariableType), POINTER :: controlLoopVariable
   
    ENTERS("ControlLoop_FieldVariableAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoopFieldVariables)) CALL FlagError("Control loop field variables is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(fieldVariable)) CALL FlagError("Field variable is not associated.",err,error,*999)
    
    !See if we already have this field variable in the list for the control loop
    found=.FALSE.
    DO variableIdx=1,controlLoopFieldVariables%numberOfFieldVariables
      controlLoopVariable=>controlLoopFieldVariables%fieldVariables(variableIdx)%fieldVariable
      IF(ASSOCIATED(fieldVariable,controlLoopVariable)) THEN
        found=.TRUE.
        EXIT
      ENDIF
    ENDDO !variableIdx
    IF(found) THEN
      !We have found the variable. Check if the time dependence and nonlinearity needs to be updated.
      IF(variableTimeDependence>controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence) &
        & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence=variableTimeDependence
      IF(variableLinearity>controlLoopFieldVariables%fieldVariables(variableIdx)%linearity) &
        & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence=variableLinearity
    ELSE
      !We have not found the field variable in the control loop list so add it to the list.
      !Reallocate the list to take the new variable
      ALLOCATE(newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate new field variables.",err,error,*999)
      DO variableIdx=1,controlLoopFieldVariables%numberOfFieldVariables
        newFieldVariables(variableIdx)%fieldVariable=> &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%fieldVariable
        newFieldVariables(variableIdx)%timeDependence= &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%timeDependence
        newFieldVariables(variableIdx)%linearity= &
          & controlLoopFieldVariables%fieldVariables(variableIdx)%linearity
      ENDDO !variableIdx
      !Add in the field variable and it's information
      CALL ControlLoop_FieldVariableInitialise(newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1),err,error,*999)
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%fieldVariable=>fieldVariable
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%timeDependence=variableTimeDependence
      newFieldVariables(controlLoopFieldVariables%numberOfFieldVariables+1)%linearity=variableLinearity
      !Move alloc the new list
      CALL MOVE_ALLOC(newFieldVariables,controlLoopFieldVariables%fieldVariables)
      controlLoopFieldVariables%numberOfFieldVariables=controlLoopFieldVariables%numberOfFieldVariables+1
    ENDIF
     
    EXITS("ControlLoop_FieldVariableAdd")
    RETURN
999 IF(ALLOCATED(newFieldVariables)) DEALLOCATE(newFieldVariables)
    ERRORSEXITS("ControlLoop_FieldVariableAdd",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariableAdd

  !
  !================================================================================================================================
  !

  !>Calculate the dependent field variables involved in a control loop
  RECURSIVE SUBROUTINE ControlLoop_FieldVariablesCalculate(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to calculate the field variables for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx,numberOfSolverMatrices,numberOfVariables,solverIdx,solverMatrixIdx,variableIdx,variableLinearity, &
      & variableTimeDependence
    TYPE(ControlLoopType), POINTER :: controlLoop2
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(FieldType), POINTER :: field
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("ControlLoop_FieldVariablesCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)

    IF(controlLoop%numberOfSubLoops>0) THEN
      !We have sub loops so recursively calculate the field variables for the underlying control loops
      DO loopIdx=1,controlLoop%numberOfSubLoops
        controlLoop2=>controlLoop%subLoops(loopIdx)%ptr
        CALL ControlLoop_FieldVariablesCalculate(controlLoop2,err,error,*999)
      ENDDO !loopIdx
      !Add all the variables from the sub-loops to this loop
      !Initialise this control loop field variables
      CALL ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*999)
      DO loopIdx=1,controlLoop%numberOfSubLoops
        controlLoop2=>controlLoop%subLoops(loopIdx)%ptr
        IF(.NOT.ASSOCIATED(controlLoop2%fieldVariables)) THEN
          localError="Field variables is not associated for sub loop index "// &
            & TRIM(NumberToVString(loopIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO variableIdx=1,controlLoop2%fieldVariables%numberOfFieldVariables
          fieldVariable=>controlLoop2%fieldVariables%fieldVariables(variableIdx)%fieldVariable
          variableLinearity=controlLoop2%fieldVariables%fieldVariables(variableIdx)%linearity
          variableTimeDependence=controlLoop2%fieldVariables%fieldVariables(variableIdx)%timeDependence
          CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
            & fieldVariable,err,error,*999)
        ENDDO !variableIdx
      ENDDO !loopIdx
    ELSE
      !There are no sub loops so process the solvers for this loop to calculate the field variables
      !Initialise this control loop field variables
      CALL ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      !Loop over the solvers
      DO solverIdx=1,solvers%numberOfSolvers
        !Get the solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
        solverEquations=>solver%solverEquations
!!TODO: Need to think about solvers that do not have solver equations.
        IF(ASSOCIATED(solverEquations)) THEN
          !If we have solver equations then find the variables.
!!TODO: could flag the linearity and time dependence of variable in equations.
          SELECT CASE(solverEquations%linearity)
          CASE(SOLVER_EQUATIONS_LINEAR)
            variableLinearity=CONTROL_LOOP_FIELD_VARIABLE_LINEAR
          CASE(SOLVER_EQUATIONS_NONLINEAR)
            variableLinearity=CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR
          CASE DEFAULT
            localError="The solver equations linearity type of "//TRIM(NumberToVString(solverEquations%linearity,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          SELECT CASE(solverEquations%timeDependence)
          CASE(SOLVER_EQUATIONS_STATIC)
            variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_STATIC
          CASE(SOLVER_EQUATIONS_QUASISTATIC)
            variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC
          CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
            dynamicSolver=>solver%dynamicSolver
            IF(ASSOCIATED(dynamicSolver)) THEN
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE)  THEN
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE)  THEN
                  variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC
                ELSE
                  variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC
                ENDIF
              ELSE
                variableTimeDependence=CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC
              ENDIF
            ELSE
              CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The solver equations time dependence type of "// &
              & TRIM(NumberToVString(solverEquations%timeDependence,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Get the solver mapping
          NULLIFY(solverMatrices)
          CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverMatrices_SolverMappingGet(solverMatrices,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*999)
          !Loop over the solver matrices
          DO solverMatrixIdx=1,numberOfSolverMatrices
            NULLIFY(solverMatrixToEquationsMap)
            CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
              & err,error,*999)
            NULLIFY(solverMappingVariables)
            CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
            !Loop over the field variables associated with the solver mapping
            CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
            DO variableIdx=1,numberOfVariables
              NULLIFY(solverMappingVariable)
              CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
              NULLIFY(fieldVariable)
              CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
              CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
                & fieldVariable,err,error,*999)
            ENDDO !variableIdx
          ENDDO !solverMatrixIdx
          !Add in the RHS
          NULLIFY(solverMappingVariables)
          CALL SolverMapping_RHSVariablesListGet(solverMapping,solverMappingVariables,err,error,*999)
          CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(solverMappingVariable)
            CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
            NULLIFY(fieldVariable)
            CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
            CALL ControlLoop_FieldVariableAdd(controlLoop%fieldVariables,variableLinearity,variableTimeDependence, &
              & fieldVariable,err,error,*999)
          ENDDO !equationSetIdx
        ENDIF
      ENDDO !solverIdx
    ENDIF

    IF(diagnostics1) THEN
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Control loop field variables:",err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Loop level = ",controlLoop%controlLoopLevel,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Sub loop index = ",controlLoop%subLoopIndex,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Label = ",controlLoop%label,err,error,*999)
      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Loop type = ",controlLoop%loopType,err,error,*999)
      CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"  Field Variables:",err,error,*999)
      IF(ASSOCIATED(controlLoop%fieldVariables)) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of field variables = ",controlLoop%fieldVariables% &
          & numberOfFieldVariables,err,error,*999)
        IF(ALLOCATED(controlLoop%fieldVariables%fieldVariables)) THEN
          DO variableIdx=1,controlLoop%fieldVariables%numberOfFieldVariables          
            fieldVariable=>controlLoop%fieldVariables%fieldVariables(variableIdx)%fieldVariable
            IF(ASSOCIATED(fieldVariable)) THEN
              field=>fieldVariable%field
              IF(ASSOCIATED(field)) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Variable index : ",variableIdx,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable field user number = ",field%userNumber, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable type = ",fieldVariable%variableType, &
                  & err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable linearity = ",controlLoop%fieldVariables% &
                  & fieldVariables(variableIdx)%linearity,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Variable time dependence = ",controlLoop%fieldVariables% &
                  & fieldVariables(variableIdx)%timeDependence,err,error,*999)
              ELSE
                localError="Field variable field is not associated for variable index "// &
                  & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            ELSE
              localError="Field variable is not associated for variable index "// &
                & TRIM(NumberToVString(variableIdx,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)              
            ENDIF
          ENDDO !variableIdx
        ELSE
          CALL FlagError("Control loop field variables is not allocated.",err,error,*999)
        ENDIF
      ELSE
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Number of field variables = ",0,err,error,*999)
      ENDIF      
    ENDIF
    
    EXITS("ControlLoop_FieldVariablesCalculate")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariablesCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesCalculate

  !
  !================================================================================================================================
  !

  !>Finalises a ControlLoopFieldVariablesType and deallocates all memory
  SUBROUTINE ControlLoop_FieldVariablesFinalise(controlLoopFieldVariables,err,error,*)

    !Argument variables
    TYPE(ControlLoopFieldVariablesType), POINTER :: controlLoopFieldVariables !<The control loop field variables to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: fieldVariableIdx
   
    ENTERS("ControlLoop_FieldVariablesFinalise",err,error,*999)

    IF(ASSOCIATED(controlLoopFieldVariables)) THEN
      DO fieldVariableIdx=1,SIZE(controlLoopFieldVariables%fieldVariables,1)
        CALL ControlLoop_FieldVariableFinalise(controlLoopFieldVariables%fieldVariables(fieldVariableIdx),err,error,*998)
      ENDDO !fieldVariableIdx
998   DEALLOCATE(controlLoopFieldVariables%fieldVariables)
      controlLoopFieldVariables%numberOfFieldVariables=0
    ENDIF
      
    EXITS("ControlLoop_FieldVariablesFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_FieldVariablesFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a ControlLoopFieldVariableType for a controlLoop
  SUBROUTINE ControlLoop_FieldVariablesInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<The control loop to initialsie the field variables for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
   
    ENTERS("ControlLoop_FieldVariablesInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    
    IF(ASSOCIATED(controlLoop%fieldVariables)) CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,err,error,*999)
    ALLOCATE(controlLoop%fieldVariables,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate control loop field variables.",err,error,*999)
    controlLoop%fieldVariables%numberOfFieldVariables=0
      
    EXITS("ControlLoop_FieldVariablesInitialise")
    RETURN
999 CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_FieldVariablesInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FieldVariablesInitialise

  !
  !=================================================================================================================================
  !

  !>Finalise a control loop and deallocate all memory
  RECURSIVE SUBROUTINE ControlLoop_Finalise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2
 
    ENTERS("ControlLoop_Finalise",err,error,*999)

    IF(ASSOCIATED(controlLoop)) THEN
      !Finalise any sub control loops first
      IF(controlLoop%numberOfSubLoops>0) THEN
        DO loopIdx=1,controlLoop%numberOfSubLoops
          controlLoop2=>controlLoop%subLoops(loopIdx)%ptr
          CALL ControlLoop_Finalise(controlLoop2,err,error,*999)
        ENDDO !loopIdx
        DEALLOCATE(controlLoop%subLoops)
      ENDIF
      !Finalise any solvers
      IF(ASSOCIATED(controlLoop%solvers)) CALL Solvers_Destroy(controlLoop%solvers,err,error,*999)
      !Now finalise this control loop
      controlLoop%label=""
      CALL ControlLoop_SimpleFinalise(controlLoop%simpleLoop,err,error,*999)
      CALL ControlLoop_FixedFinalise(controlLoop%fixedLoop,err,error,*999)
      CALL ControlLoop_LoadIncrementFinalise(controlLoop%loadIncrementLoop,err,error,*999)
      CALL ControlLoop_TimeFinalise(controlLoop%timeLoop,err,error,*999)
      CALL ControlLoop_WhileFinalise(controlLoop%whileLoop,err,error,*999)
      CALL ControlLoop_FieldVariablesFinalise(controlLoop%fieldVariables,err,error,*999)
      DEALLOCATE(controlLoop)
    ENDIF
       
    EXITS("ControlLoop_Finalise")
    RETURN
999 ERRORSEXITS("ControlLoop_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Finalise

  !
  !=================================================================================================================================
  !

  !>Initialise a control loop.
  SUBROUTINE ControlLoop_Initialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("ControlLoop_Initialise",err,error,*998)

    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is already associated.",err,error,*998)
    
    ALLOCATE(controlLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate control loop.",err,error,*999)
    NULLIFY(controlLoop%problem)
    NULLIFY(controlLoop%parentLoop)
    controlLoop%controlLoopFinished=.FALSE.
    controlLoop%label=" "
    controlLoop%controlLoopLevel=0
    controlLoop%subLoopIndex=0
    controlLoop%outputType=CONTROL_LOOP_NO_OUTPUT
    NULLIFY(controlLoop%simpleLoop)
    NULLIFY(controlLoop%fixedLoop)
    NULLIFY(controlLoop%timeLoop)
    NULLIFY(controlLoop%whileLoop)
    NULLIFY(controlLoop%loadIncrementLoop)
    controlLoop%numberOfSubLoops=0
    NULLIFY(controlLoop%fieldVariables)
    NULLIFY(controlLoop%solvers)
    controlLoop%loopType=CONTROL_SIMPLE_TYPE
    CALL ControlLoop_SimpleInitialise(controlLoop,err,error,*999)
              
    EXITS("ControlLoop_Initialise")
    RETURN
999 CALL ControlLoop_Finalise(controlLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_Initialise

  !
  !================================================================================================================================
  !

  !>Finalises a fixed control loop and deallocates all memory.
  SUBROUTINE ControlLoop_FixedFinalise(fixedLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopFixedType), POINTER, INTENT(INOUT) :: fixedLoop !<A pointer to the fixed control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_FixedFinalise",err,error,*999)

    IF(ASSOCIATED(fixedLoop)) THEN
      DEALLOCATE(fixedLoop)
    ENDIF
       
    EXITS("ControlLoop_FixedFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_FixedFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FixedFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a fixed loop for a control loop.
  SUBROUTINE ControlLoop_FixedInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise the fixed loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("ControlLoop_FixedInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%fixedLoop)) &
      & CALL FlagError("The fixed loop is already associated for this control loop.",err,error,*998)
    
    ALLOCATE(controlLoop%fixedLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate fixed loop for the control loop.",err,error,*999)
    controlLoop%fixedLoop%controlLoop=>controlLoop
    controlLoop%fixedLoop%iterationNumber=0
    controlLoop%fixedLoop%startIteration=1
    controlLoop%fixedLoop%stopIteration=100
    controlLoop%fixedLoop%iterationIncrement=1
    
    EXITS("ControlLoop_FixedInitialise")
    RETURN
999 CALL ControlLoop_FixedFinalise(controlLoop%fixedLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_FixedInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_FixedInitialise

  !
  !================================================================================================================================
  !

  !>Sets the iteration parameters for a fixed control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_IterationsSet
  SUBROUTINE ControlLoop_IterationsSet(controlLoop,startIteration,stopIteration,iterationIncrement,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to fixed control loop to set the iterations for
    INTEGER(INTG), INTENT(IN) :: startIteration !<The start iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: stopIteration !<The stop iteration for the fixed control loop.
    INTEGER(INTG), INTENT(IN) :: iterationIncrement !<The iteration increment for the fixed control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopFixedType), POINTER :: fixedLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_IterationsSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsFixedLoop(controlLoop,err,error,*999)
    IF(iterationIncrement==0) THEN
      localError="The specified time increment of "//TRIM(NumberToVString(iterationIncrement,"*",err,error))// &
        & " is invalid. The iteration increment must not be zero."          
      CALL FlagError(localError,err,error,*999)
    ELSE
      IF(iterationIncrement>0) THEN
        IF(stopIteration<=startIteration) THEN
          localError="The specified stop iteration of "//TRIM(NumberToVString(stopIteration,"*",err,error))// &
            & " is incompatiable with a specified start increment of "// &
            & TRIM(NumberToVString(startIteration,"*",err,error))// &
            & ". For a positive iteration increment the stop iteration must be > than the start iteration."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        IF(startIteration<=stopIteration) THEN
          localError="The specified start iteration of "//TRIM(NumberToVString(startIteration,"*",err,error))// &
            & " is incompatiable with a specified stop iteration of "// &
            & TRIM(NumberToVString(stopIteration,"*",err,error))// &
            & ". For a negative iteration increment the stop iteration must be < than the start iteration."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF

    NULLIFY(fixedLoop)
    CALL ControlLoop_FixedLoopGet(controlLoop,fixedLoop,err,error,*999)
    
    fixedLoop%startIteration=startIteration
    fixedLoop%stopIteration=stopIteration
    fixedLoop%iterationIncrement=iterationIncrement
    
    EXITS("ControlLoop_IterationsSet")
    RETURN
999 ERRORSEXITS("ControlLoop_IterationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_IterationsSet
  
  !
  !================================================================================================================================
  !

  !>Returns the label of a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_LabelGet
  SUBROUTINE ControlLoop_LabelGetC(controlLoop,label,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to get the label for
    CHARACTER(LEN=*), INTENT(OUT) :: label !<On return, the control loop label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cLength,vsLength

    ENTERS("ControlLoop_LabelGetC",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    cLength=LEN(label)
    vsLength=LEN_TRIM(controlLoop%label)
    IF(cLength>vsLength) THEN
      label=CHAR(controlLoop%label,vsLength)
    ELSE
      label=CHAR(controlLoop%label,cLength)
    ENDIF
    
    EXITS("ControlLoop_LabelGetC")
    RETURN
999 ERRORSEXITS("ControlLoop_LabelGetC",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LabelGetC

   !
  !================================================================================================================================
  !

  !>Returns the label of a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_LabelGet
  SUBROUTINE ControlLoop_LabelGetVS(controlLoop,label,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to get the label for
    TYPE(VARYING_STRING), INTENT(OUT) :: label !<On return, the control loop label.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ControlLoop_LabelGetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    label=VAR_STR(CHAR(controlLoop%label))
     
    EXITS("ControlLoop_LabelGetVS")
    RETURN
999 ERRORSEXITS("ControlLoop_LabelGetVS",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LabelGetVS

  !
  !================================================================================================================================
  !

  !>Sets the label of a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_LabelSet
  SUBROUTINE ControlLoop_LabelSetC(controlLoop,label,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ControlLoop_LabelSetC",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    
    controlLoop%label=label
    
    EXITS("ControlLoop_LabelSetC")
    RETURN
999 ERRORSEXITS("ControlLoop_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_LabelSet
  SUBROUTINE ControlLoop_LabelSetVS(controlLoop,label,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("ControlLoop_LabelSetVS",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    
    controlLoop%label=label
    
    EXITS("ControlLoop_LabelSetVS")
    RETURN
999 ERRORSEXITS("ControlLoop_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Sets the maximum number of iterations for a while or load increment control loop. \see OpenCMISS_cmfe_ControlLoop_MaximumIterationsSet
  SUBROUTINE ControlLoop_MaximumIterationsSet(controlLoop,maximumIterations,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to the while or load incremented control loop to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: maximumIterations !<The maximum number of iterations for the while or load incremented control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_MaximumIterationsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    IF(maximumIterations<=0) THEN
      localError="The specified maximum number of iterations of "// &
        & TRIM(NumberToVString(maximumIterations,"*",err,error))// &
        & " is invalid. The maximum number of iterations must be > 0."          
      CALL FlagError(localError,err,error,*999)            
    ENDIF

    IF(controlLoop%controlLoopFinished) THEN
      IF(controlLoop%loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
        !allow to update the maximum number of iterations at a later time for the load increment loop type
        NULLIFY(loadIncrementLoop)
        CALL ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*999)
        loadIncrementLoop%maximumNumberOfIterations=maximumIterations
      ELSE
        CALL FlagError("Control loop has already been finished.",err,error,*999)
      ENDIF
    ELSE
      IF(controlLoop%loopType==CONTROL_WHILE_LOOP_TYPE) THEN
        NULLIFY(whileLoop)
        CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,err,error,*999)
        whileLoop%maximumNumberOfIterations=maximumIterations
      ELSE IF(controlLoop%loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
        NULLIFY(loadIncrementLoop)
        CALL ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*999)
        loadIncrementLoop%maximumNumberOfIterations=maximumIterations
      ELSE
        CALL FlagError("The specified control loop is not a while or load increment control loop.",err,error,*999)
      ENDIF
    ENDIF
       
    EXITS("ControlLoop_MaximumIterationsSet")
    RETURN
999 ERRORSEXITS("ControlLoop_MaximumIterationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_MaximumIterationsSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the output for a load incremented control loop identified by an object. \see OpenCMISS_cmfe_ControlLoop_LoadOutputSet
  SUBROUTINE ControlLoop_LoadOutputSet(controlLoop,outputFrequency,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to the load incremented control loop to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: outputFrequency !<The output frequency modulo to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
 
    ENTERS("ControlLoop_LoadOutputSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsLoadIncrementLoop(controlLoop,err,error,*999)
    
    NULLIFY(loadIncrementLoop)
    CALL ControlLoop_LoadIncrementLoopGet(controlLoop,loadIncrementLoop,err,error,*999)
    loadIncrementLoop%outputNumber=outputFrequency
        
    EXITS("ControlLoop_LoadOutputSet")
    RETURN
999 ERRORSEXITS("ControlLoop_LoadOutputSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LoadOutputSet

  !
  !================================================================================================================================
  !

  !>Sets the absolute tolerance (convergence condition tolerance) for a while control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_AbsoluteToleranceSet
  SUBROUTINE ControlLoop_AbsoluteToleranceSet(controlLoop,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to while control loop to set the maximum iterations for
    REAL(DP), INTENT(IN) :: absoluteTolerance !<The absolute tolerance value for a while control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_AbsoluteToleranceSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsWhileLoop(controlLoop,err,error,*999)
    IF(absoluteTolerance<=0) THEN
      localError="The specified absolute tolerance of "// &
        & TRIM(NumberToVString(absoluteTolerance,"*",err,error))// &
        & " is invalid for a while loop. The tolerance must be greater than zero."          
      CALL FlagError(localError,err,error,*999)            
    ENDIF

    NULLIFY(whileLoop)
    CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,err,error,*999)
    whileLoop%absoluteTolerance=absoluteTolerance
       
    EXITS("ControlLoop_AbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("ControlLoop_AbsoluteToleranceSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_AbsoluteToleranceSet

  !
  !================================================================================================================================
  !

  !>Sets the relative tolerance (convergence condition tolerance) for a while control loop. \see OpenCMISS::Iron::cmfe_ControlLoopRelativeToleranceSet
  SUBROUTINE ControlLoop_RelativeToleranceSet(controlLoop,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to while control loop to set the maximum iterations for
    REAL(DP), INTENT(IN) :: relativeTolerance !<The relative tolerance value for a while control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_RelativeToleranceSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsWhileLoop(controlLoop,err,error,*999)
    IF(relativeTolerance<=0) THEN
      localError="The specified relative tolerance of "// &
        & TRIM(NumberToVString(relativeTolerance,"*",err,error))// &
        & " is invalid for a while loop. The tolerance must be greater than zero."          
      CALL FlagError(localError,err,error,*999)            
    ENDIF

    NULLIFY(whileLoop)
    CALL ControlLoop_WhileLoopGet(controlLoop,whileLoop,err,error,*999)
    whileLoop%relativeTolerance=relativeTolerance
       
    EXITS("ControlLoop_RelativeToleranceSet")
    RETURN
999 ERRORSEXITS("ControlLoop_RelativeToleranceSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_RelativeToleranceSet

  !
  !================================================================================================================================
  !

  !>Sets the number of iterations for a time type control loop. If set to 0 (default), it will be computed by start and stop time and time increment. \see OpenCMISS_ControlLoop_NumberOfIterationsSet
  SUBROUTINE ControlLoop_NumberOfIterationsSet(controlLoop,numberOfIterations,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to time control loop to set the number of iterations for
    INTEGER(INTG), INTENT(IN) :: numberOfIterations !<The number of iterations for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("ControlLoop_NumberOfIterationsSet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    CALL ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*999)
    IF(numberOfIterations<0) THEN
      localError="The specified number of iterations of "//TRIM(NumberToVString(numberOfIterations,"*",err,error))// &
        & " is invalid. The number of iterations must be >= 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(timeLoop)
    CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
    timeLoop%numberOfIterations=numberOfIterations          
    !Update time increment if number of iterations differs from time stepping settings
    IF (CEILING((timeLoop%stopTime-timeLoop%startTime)/timeLoop%timeIncrement) /= timeLoop%numberOfIterations) &
      & timeLoop%timeIncrement = (timeLoop%stopTime-timeLoop%startTime)/timeLoop%numberOfIterations
       
    EXITS("ControlLoop_NumberOfIterationsSet")
    RETURN
999 ERRORSEXITS("ControlLoop_NumberOfIterationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_NumberOfIterationsSet
  
  !
  !================================================================================================================================
  !
  
  !>Sets/changes the number of sub loops in a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_NumberOfSubLoopsSet
  SUBROUTINE ControlLoop_NumberOfSubLoopsSet(controlLoop,numberOfSubLoops,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to control loop to set the number of sub loops for
    INTEGER(INTG), INTENT(IN) :: numberOfSubLoops !<The number of sub loops to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx
    TYPE(ControlLoopPtrType), ALLOCATABLE :: oldSubLoops(:)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ControlLoop_NumberOfSubLoopsSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    IF(numberOfSubLoops<0) THEN
      localError="The specified number of sub loops of "//TRIM(NumberToVString(numberOfSubLoops,"*",err,error))// &
        & " is invalid. The number of sub loops must be >= 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(numberOfSubLoops/=controlLoop%numberOfSubLoops) THEN
      IF(controlLoop%numberOfSubLoops>0) THEN
        ALLOCATE(oldSubLoops(controlLoop%numberOfSubLoops),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate old sub loops.",err,error,*999)
        DO loopIdx=1,controlLoop%numberOfSubLoops
          oldSubLoops(loopIdx)%ptr=>controlLoop%subLoops(loopIdx)%ptr
        ENDDO !loopIdx
        DEALLOCATE(controlLoop%subLoops)
      ENDIF
      ALLOCATE(controlLoop%subLoops(numberOfSubLoops),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate control loop sub loops.",err,error,*999)
      IF(numberOfSubLoops>controlLoop%numberOfSubLoops) THEN
        DO loopIdx=1,controlLoop%numberOfSubLoops
          controlLoop%subLoops(loopIdx)%ptr=>oldSubLoops(loopIdx)%ptr
        ENDDO !loopIdx
        DO loopIdx=controlLoop%numberOfSubLoops+1,numberOfSubLoops
          NULLIFY(controlLoop%subLoops(loopIdx)%ptr)
          CALL ControlLoop_Initialise(controlLoop%subLoops(loopIdx)%ptr,err,error,*999)
          controlLoop%subLoops(loopIdx)%ptr%problem=>controlLoop%problem
          controlLoop%subLoops(loopIdx)%ptr%parentLoop=>controlLoop
          controlLoop%subLoops(loopIdx)%ptr%controlLoopLevel=controlLoop%controlLoopLevel+1
          controlLoop%subLoops(loopIdx)%ptr%subLoopIndex=loopIdx
        ENDDO !loopIdx
      ELSE
        DO loopIdx=1,numberOfSubLoops
          controlLoop%subLoops(loopIdx)%ptr=>oldSubLoops(loopIdx)%ptr
        ENDDO !loopIdx
        DO loopIdx=numberOfSubLoops+1,controlLoop%numberOfSubLoops
          CALL ControlLoop_Finalise(oldSubLoops(loopIdx)%ptr,err,error,*999)
        ENDDO !loopIdx
      ENDIF
      IF(ALLOCATED(oldSubLoops)) DEALLOCATE(oldSubLoops)
      controlLoop%numberOfSubLoops=numberOfSubLoops
    ENDIF
    
    EXITS("ControlLoop_NumberOfSubLoopsSet")
    RETURN
999 IF(ALLOCATED(oldSubLoops)) DEALLOCATE(oldSubLoops)
    ERRORSEXITS("ControlLoop_NumberOfSubLoopsSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_NumberOfSubLoopsSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for a control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_OutputTypeSet
  SUBROUTINE ControlLoop_OutputTypeSet(controlLoop,outputType,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer the control loop to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The type of control loop output to be set \see ControlLoopRoutines_OutputTypes,ControlLoopRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ControlLoop_OutputTypeSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)

    SELECT CASE(outputType)
    CASE(CONTROL_LOOP_NO_OUTPUT)
      controlLoop%outputType=CONTROL_LOOP_NO_OUTPUT
    CASE(CONTROL_LOOP_PROGRESS_OUTPUT)
      controlLoop%outputType=CONTROL_LOOP_PROGRESS_OUTPUT
    CASE(CONTROL_LOOP_TIMING_OUTPUT)
      controlLoop%outputType=CONTROL_LOOP_TIMING_OUTPUT
    CASE DEFAULT
      localError="The specified control loop output type of "// &
        & TRIM(NumberToVString(outputType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("ControlLoop_OutputTypeSet")
    RETURN
999 ERRORSEXITS("ControlLoop_OutputTypeSet",err,error)    
    RETURN 1
   
  END SUBROUTINE ControlLoop_OutputTypeSet
        
  !
  !================================================================================================================================
  !

  !>Updates the previous values for dependent variables under a time control loop
  SUBROUTINE ControlLoop_PreviousValuesUpdate(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer the time control loop to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearity,timeDependence,variableIdx
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("ControlLoop_PreviousValuesUpdate",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)    
    IF(.NOT.ASSOCIATED(controlLoop%fieldVariables)) &
      & CALL FlagError("Control loop field variables is not associated.",err,error,*999)
    
    DO variableIdx=1,controlLoop%fieldVariables%numberOfFieldVariables
      fieldVariable=>controlLoop%fieldVariables%fieldVariables(variableIdx)%fieldVariable
      linearity=controlLoop%fieldVariables%fieldVariables(variableIdx)%linearity
      timeDependence=controlLoop%fieldVariables%fieldVariables(variableIdx)%timeDependence
      IF(.NOT.ASSOCIATED(fieldVariable)) THEN
        localError="The field variable for variable index "//TRIM(NumberToVString(variableIdx,"*",err,error))// &
          & " is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      SELECT CASE(timeDependence)
      CASE(CONTROL_LOOP_FIELD_VARIABLE_STATIC)
        !Do nothing additional
      CASE(CONTROL_LOOP_FIELD_VARIABLE_QUASISTATIC)
        !Copy field variable values to the previous values
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & 1.0_DP,err,error,*999)
        SELECT CASE(linearity)
        CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
          !Do nothing additional
        CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
          !Copy residual to the previous residual
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
            & 1.0_DP,err,error,*999)
        CASE DEFAULT
          localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_LOOP_FIELD_VARIABLE_FIRST_DEGREE_DYNAMIC)
        !Copy field variable values to the previous values
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS2_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS3_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS2_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        !Copy velocity values to the previous velocity
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VELOCITY_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
        SELECT CASE(linearity)
        CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
          !Do nothing additional
        CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
          !Copy residuals to the previous residuals
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS2_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS3_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS2_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
        CASE DEFAULT
          localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_LOOP_FIELD_VARIABLE_SECOND_DEGREE_DYNAMIC)
        !Copy field variable values to the previous values
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS2_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS3_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS2_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        !Copy velocity values to the previous velocity
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_VELOCITY_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
        !Copy acceleration values to the previous velocity
        CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
        SELECT CASE(linearity)
        CASE(CONTROL_LOOP_FIELD_VARIABLE_LINEAR)
          !Do nothing additional
        CASE(CONTROL_LOOP_FIELD_VARIABLE_NONLINEAR)
          !Copy residuals to the previous residuals
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS2_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS3_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS2_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
          CALL FieldVariable_ParameterSetsCopyIfExists(fieldVariable,FIELD_RESIDUAL_SET_TYPE, &
            & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
        CASE DEFAULT
          localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The control loop variable linearity of "//TRIM(NumberToVString(linearity,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDDO !variableIdx
    
    EXITS("ControlLoop_PreviousValuesUpdate")
    RETURN
999 ERRORSEXITS("ControlLoop_PreviousValuesUpdate",err,error)
    RETURN 1

  END SUBROUTINE ControlLoop_PreviousValuesUpdate

  !
  !================================================================================================================================
  !

  !>Finalises a simple control loop and deallocates all memory.
  SUBROUTINE ControlLoop_SimpleFinalise(simpleLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopSimpleType), POINTER, INTENT(INOUT) :: simpleLoop !<A pointer to the simple control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_SimpleFinalise",err,error,*999)

    IF(ASSOCIATED(simpleLoop)) THEN
      DEALLOCATE(simpleLoop)
    ENDIF
       
    EXITS("ControlLoop_SimpleFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_SimpleFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SimpleFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a simple loop for a control loop.
  SUBROUTINE ControlLoop_SimpleInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise the simple loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("ControlLoop_SimpleInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%simpleLoop)) &
      & CALL FlagError("The simple loop is already associated for this control loop.",err,error,*998)
    
    ALLOCATE(controlLoop%simpleLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate simple loop for the control loop.",err,error,*999)
    controlLoop%simpleLoop%controlLoop=>controlLoop
       
    EXITS("ControlLoop_SimpleInitialise")
    RETURN
999 CALL ControlLoop_SimpleFinalise(controlLoop%simpleLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_SimpleInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SimpleInitialise

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solvers for a control loop and all sub control loops.
  RECURSIVE SUBROUTINE ControlLoop_SolversDestroy(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to destroy the solvers for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2

    ENTERS("ControlLoop_SolversDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    
    !Destroy the solvers in any sub control loops first
    IF(controlLoop%numberOfSubLoops>0) THEN
      DO loopIdx=1,controlLoop%numberOfSubLoops
        NULLIFY(controlLoop2)
        CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
        CALL ControlLoop_SolversDestroy(controlLoop2,err,error,*999)
      ENDDO !loopIdx
    ENDIF
    !Destroy the solvers in this control loop
    IF(ASSOCIATED(controlLoop%solvers)) CALL Solvers_Destroy(controlLoop%solvers,err,error,*999)
       
    EXITS("ControlLoop_SolversDestroy")
    RETURN
999 ERRORSEXITS("ControlLoop_SolversDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SolversDestroy

  !
  !================================================================================================================================
  !

  !>Recursively destroys the solver equations for a control loop and all sub control loops. \todo Create solvers_solver_equations_destory and call?
  RECURSIVE SUBROUTINE ControlLoop_SolverEquationsDestroy(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to destroy the solver for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: loopIdx,solverIdx
    TYPE(ControlLoopType), POINTER :: controlLoop2
    TYPE(SolverType), POINTER :: solver
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquations
 
    ENTERS("ControlLoop_SolverEquationsDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    !Destroy the solver equations in any sub control loops first
    IF(controlLoop%numberOfSubLoops>0) THEN
      DO loopIdx=1,controlLoop%numberOfSubLoops
        NULLIFY(controlLoop2)
        CALL ControlLoop_SubLoopGet(controlLoop,loopIdx,controlLoop2,err,error,*999)
        CALL ControlLoop_SolverEquationsDestroy(controlLoop2,err,error,*999)
      ENDDO !loopIdx
    ENDIF
    !Destroy the solver equations in this control loop
    solvers=>controlLoop%solvers
    IF(ASSOCIATED(solvers)) THEN
      DO solverIdx=1,solvers%numberOfSolvers
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
        solverEquations=>solver%solverEquations
        IF(ASSOCIATED(solverEquations)) CALL SolverEquations_Destroy(solverEquations,err,error,*999)
      ENDDO !solverIdx
    ENDIF
       
    EXITS("ControlLoop_SolverEquationsDestroy")
    RETURN
999 ERRORSEXITS("ControlLoop_SolverEquationsDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_SolverEquationsDestroy

  !
  !================================================================================================================================
  !

  !>Finalises a time control loop and deallocates all memory.
  SUBROUTINE ControlLoop_TimeFinalise(timeLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopTimeType), POINTER, INTENT(INOUT) :: timeLoop !<A pointer to the time control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_TimeFinalise",err,error,*999)

    IF(ASSOCIATED(timeLoop)) THEN
      DEALLOCATE(timeLoop)
    ENDIF
       
    EXITS("ControlLoop_TimeFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_TimeFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimeFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a time loop for a control loop.
  SUBROUTINE ControlLoop_TimeInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise the time loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("ControlLoop_TimeInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%timeLoop)) &
      & CALL FlagError("The time loop is already associated for this control loop.",err,error,*998)
     
    ALLOCATE(controlLoop%timeLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate time loop for the control loop.",err,error,*999)
    controlLoop%timeLoop%controlLoop=>controlLoop
    controlLoop%timeLoop%iterationNumber=0
    controlLoop%timeLoop%numberOfIterations=0
    controlLoop%timeLoop%globalIterationNumber=0
    controlLoop%timeLoop%currentTime=0.0_DP
    controlLoop%timeLoop%startTime=0.0_DP
    controlLoop%timeLoop%stopTime=1.0_DP
    controlLoop%timeLoop%timeIncrement=0.01_DP
    controlLoop%timeLoop%outputNumber=0
    controlLoop%timeLoop%inputNumber=0
       
    EXITS("ControlLoop_TimeInitialise")
    RETURN
999 CALL ControlLoop_TimeFinalise(controlLoop%timeLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_TimeInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimeInitialise

  !
  !================================================================================================================================
  !
  
  !>Gets the current time parameters for a time control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_CurrentTimesGet
  SUBROUTINE ControlLoop_TimesGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
    & currentLoopIteration,outputIterationNumber,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop
    REAL(DP), INTENT(OUT) :: startTime
    REAL(DP), INTENT(OUT) :: stopTime
    REAL(DP), INTENT(OUT) :: currentTime
    REAL(DP), INTENT(OUT) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: currentLoopIteration
    INTEGER(INTG), INTENT(OUT) :: outputIterationNumber
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables    
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(ControlLoopType), POINTER :: parentLoop
    INTEGER(INTG), POINTER :: controlLoopLevel
    INTEGER(INTG) :: i

    ENTERS("ControlLoop_TimesGet",err,error,*999)

    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)
    
    controlLoopLevel=>controlLoop%controlLoopLevel
    parentLoop=>controlLoop
    DO i=controlLoopLevel,1,-1
      IF(controlLoopLevel==0) CALL FlagError("The specified control loop is not a time control loop.",err,error,*999)      
      IF(parentLoop%loopType==CONTROL_TIME_LOOP_TYPE) THEN
        NULLIFY(timeLoop)
        CALL ControlLoop_TimeLoopGet(parentLoop,timeLoop,err,error,*999)
        startTime=timeLoop%startTime
        stopTime=timeLoop%stopTime
        currentTime=timeLoop%currentTime
        timeIncrement=timeLoop%timeIncrement
        currentLoopIteration=timeLoop%iterationNumber
        outputIterationNumber=timeLoop%outputNumber
        EXIT
      ELSE
        parentLoop=>parentLoop%parentLoop
      ENDIF
    ENDDO !i
      
    EXITS("ControlLoop_TimesGet")
    RETURN
999 ERRORSEXITS("ControlLoop_TimesGet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimesGet

  !
  !================================================================================================================================
  !

  !>Sets the time parameters for a time control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_TimesSet
  SUBROUTINE ControlLoop_TimesSet(controlLoop,startTime,stopTime,timeIncrement,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to set the times for
    REAL(DP), INTENT(IN) :: startTime !<The start time for the time control loop.
    REAL(DP), INTENT(IN) :: stopTime !<The stop time for the time control loop.
    REAL(DP), INTENT(IN) :: timeIncrement !<The time increment for the time control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ControlLoop_TimesSet",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    CALL ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*999)
    NULLIFY(timeLoop)
    CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
    IF(ABS(timeIncrement)<=ZERO_TOLERANCE) THEN
      localError="The specified time increment of "//TRIM(NumberToVString(timeIncrement,"*",err,error))// &
        & " is invalid. The time increment must not be zero."          
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(timeIncrement>0.0_DP) THEN
      IF(stopTime<=startTime) THEN
        localError="The specified stop time of "//TRIM(NumberToVString(stopTime,"*",err,error))// &
          & " is incompatiable with a specified start time of "//TRIM(NumberToVString(startTime,"*",err,error))// &
          & ". For a positive time increment the stop time must be > than the start time."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      IF(startTime<=stopTime) THEN
        localError="The specified start time of "//TRIM(NumberToVString(startTime,"*",err,error))// &
          & " is incompatiable with a specified stop time of "//TRIM(NumberToVString(stopTime,"*",err,error))// &
          & ". For a negative time increment the stop time must be < than the start time."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
   
    timeLoop%startTime=startTime
    timeLoop%stopTime=stopTime
    timeLoop%timeIncrement=timeIncrement
    timeLoop%numberOfIterations=0    ! reset number of iterations
       
    EXITS("ControlLoop_TimesSet")
    RETURN
999 ERRORSEXITS("ControlLoop_TimesSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimesSet

  !
  !================================================================================================================================
  !

  !>Sets the output parameters for a time control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_TimeOutputSet
  SUBROUTINE ControlLoop_TimeOutputSet(controlLoop,outputFrequency,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to set the times for
    INTEGER(INTG) :: outputFrequency !<The output frequency modulo to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    
    ENTERS("ControlLoop_TimeOutputSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*999)
    NULLIFY(timeLoop)
    CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)
    IF(outputFrequency<0) THEN
      CALL FlagError("Invalid output frequency. The frequency should be greater than or equal to zero, but is "// &
        & TRIM(NumberToVString(outputFrequency,"*",err,error))//".",err,error,*999)
    END IF

    timeLoop%outputNumber=outputFrequency
      
    EXITS("ControlLoop_TimeOutputSet")
    RETURN
999 ERRORSEXITS("ControlLoop_TimeOutputSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimeOutputSet

  !
  !================================================================================================================================
  !

  !>Sets the input parameters for a time control loop. \see OpenCMISS::Iron::cmfe_ControlLoop_TimeInputSet
  SUBROUTINE ControlLoop_TimeInputSet(controlLoop,inputOption,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(IN) :: controlLoop !<A pointer to control loop to set the times for
    INTEGER(INTG) :: inputOption !<The input option modulo to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    TYPE(ControlLoopTimeType), POINTER :: timeLoop
    
    ENTERS("ControlLoop_TimeInputSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)
    CALL ControlLoop_AssertIsTimeLoop(controlLoop,err,error,*999)
    NULLIFY(timeLoop)
    CALL ControlLoop_TimeLoopGet(controlLoop,timeLoop,err,error,*999)

    timeLoop%inputNumber=inputOption
       
    EXITS("ControlLoop_TimeInputSet")
    RETURN
999 ERRORSEXITS("ControlLoop_TimeInputSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TimeInputSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the control loop type.
  SUBROUTINE ControlLoop_TypeSet(controlLoop,loopType,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to control loop to set the type of
    INTEGER(INTG), INTENT(IN) :: loopType !<The type of loop type to set \see ProblemRoutines_ControlLoopTypes,ProblemRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("ControlLoop_TypeSet",err,error,*999)

    CALL ControlLoop_AssertNotFinished(controlLoop,err,error,*999)

    IF(loopType/=controlLoop%loopType) THEN
      !Initialise the new loop type
      SELECT CASE(loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        CALL ControlLoop_SimpleInitialise(controlLoop,err,error,*999)
      CASE(CONTROL_FIXED_LOOP_TYPE)
        CALL ControlLoop_FixedInitialise(controlLoop,err,error,*999)
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL ControlLoop_TimeInitialise(controlLoop,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL ControlLoop_WhileInitialise(controlLoop,err,error,*999)
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        CALL ControlLoop_LoadIncrementInitialise(controlLoop,err,error,*999)
      CASE DEFAULT
        localError="The loop type of "//TRIM(NumberToVString(loopType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finialise the old loop type
      SELECT CASE(controlLoop%loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        CALL ControlLoop_SimpleFinalise(controlLoop%simpleLoop,err,error,*999)
      CASE(CONTROL_FIXED_LOOP_TYPE)
        CALL ControlLoop_FixedFinalise(controlLoop%fixedLoop,err,error,*999)
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL ControlLoop_TimeFinalise(controlLoop%timeLoop,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL ControlLoop_WhileFinalise(controlLoop%whileLoop,err,error,*999)
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        CALL ControlLoop_LoadIncrementFinalise(controlLoop%loadIncrementLoop,err,error,*999)
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(controlLoop%loopType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      controlLoop%loopType=loopType
    ENDIF
       
    EXITS("ControlLoop_TypeSet")
    RETURN
999 ERRORSEXITS("ControlLoop_TypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_TypeSet

  !
  !================================================================================================================================
  !

  !>Finalises a while control loop and deallocates all memory.
  SUBROUTINE ControlLoop_WhileFinalise(whileLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopWhileType), POINTER, INTENT(INOUT) :: whileLoop !<A pointer to the while control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_WhileFinalise",err,error,*999)

    IF(ASSOCIATED(whileLoop)) THEN
      DEALLOCATE(whileLoop)
    ENDIF
       
    EXITS("ControlLoop_WhileFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_WhileFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_WhileFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a while loop for a control loop.
  SUBROUTINE ControlLoop_WhileInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise the while loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("ControlLoop_WhileInitialise",err,error,*998)

    IF(ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%whileLoop)) &
      & CALL FlagError("The while loop is already associated for this control loop.",err,error,*998)
    
    ALLOCATE(controlLoop%whileLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate while loop for the control loop.",err,error,*999)
    controlLoop%whileLoop%controlLoop=>controlLoop
    controlLoop%whileLoop%iterationNumber=0
    controlLoop%whileLoop%maximumNumberOfIterations=100
    controlLoop%whileLoop%absoluteTolerance=1.0E-5_DP
    controlLoop%whileLoop%continueLoop=.TRUE.
      
    EXITS("ControlLoop_WhileInitialise")
    RETURN
999 CALL ControlLoop_WhileFinalise(controlLoop%whileLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_WhileInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_WhileInitialise

  !
  !================================================================================================================================
  !
  
  !>Finalises a load increment loop and deallocates all memory.
  SUBROUTINE ControlLoop_LoadIncrementFinalise(loadIncrementLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopLoadIncrementType), POINTER, INTENT(INOUT) :: loadIncrementLoop !<A pointer to the load increment control loop to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("ControlLoop_LoadIncrementFinalise",err,error,*999)

    IF(ASSOCIATED(loadIncrementLoop)) THEN
      DEALLOCATE(loadIncrementLoop)
    ENDIF
       
    EXITS("ControlLoop_LoadIncrementFinalise")
    RETURN
999 ERRORSEXITS("ControlLoop_LoadIncrementFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LoadIncrementFinalise

  !
  !================================================================================================================================
  !

  !>Initialises a load increment loop for a control loop.
  SUBROUTINE ControlLoop_LoadIncrementInitialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER, INTENT(INOUT) :: controlLoop !<A pointer to the control loop to initialise the load increment loop for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("ControlLoop_LoadIncrementInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%loadIncrementLoop)) &
      & CALL FlagError("The load increment loop is already associated for this control loop.",err,error,*998)
      
    ALLOCATE(controlLoop%loadIncrementLoop,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate load increment loop for the control loop.",err,error,*999)
    controlLoop%loadIncrementLoop%controlLoop=>controlLoop
    controlLoop%loadIncrementLoop%iterationNumber=0
    controlLoop%loadIncrementLoop%maximumNumberOfIterations=1 ! default is full load in one step
    controlLoop%loadIncrementLoop%outputNumber=0
    controlLoop%loadIncrementLoop%inputNumber=0
      
    EXITS("ControlLoop_LoadIncrementInitialise")
    RETURN
999 CALL ControlLoop_LoadIncrementFinalise(controlLoop%loadIncrementLoop,dummyErr,dummyError,*998)
998 ERRORSEXITS("ControlLoop_LoadIncrementInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE ControlLoop_LoadIncrementInitialise

  !
  !================================================================================================================================
  !

END MODULE ControlLoopRoutines

