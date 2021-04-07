!> \file 
!> \author Chris Bradley
!> \brief This module handles all solver routines.
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

!> This module handles all solver routines.
MODULE SolverRoutines

  USE BaseRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
#ifdef WITH_CELLML
  USE CELLML_MODEL_DEFINITION
#endif
  USE CmissCellML
  USE CellMLAccessRoutines
  USE CmissPetsc
  USE CmissPetscTypes
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE Kinds
  USE InputOutput
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingAccessRoutines
  USE interfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE Maths
  USE ProblemAccessRoutines
  USE ProfilingRoutines
  USE SolverAccessRoutines
  USE SolverMappingRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesRoutines
  USE SolverMatricesAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

#include "petscversion.h"
 
  !Module parameters
  
  !Module types

  !Module variables

  !Interfaces

  INTERFACE

    SUBROUTINE SolverDAEExternal_Integrate(numberOfDofs,startTime,endTime,initialStep, &
      & onlyOneModelIndex,modelsData,numberOfState,stateData,numberOfParameters, &
      & parametersData,numberOfIntermediate,intermediateData,err) BIND(C, NAME="SolverDAEExternalIntegrate")
      
      USE ISO_C_BINDING

      INTEGER(C_INT), VALUE, INTENT(IN) :: numberOfDofs
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: startTime
      REAL(C_DOUBLE), VALUE, INTENT(IN) :: endTime
      REAL(C_DOUBLE), INTENT(INOUT) :: initialStep
      INTEGER(C_INT), VALUE, INTENT(IN) :: onlyOneModelIndex
      INTEGER(C_INT), INTENT(IN) :: modelsData(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: numberOfState
      REAL(C_DOUBLE), INTENT(INOUT) :: stateData(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: numberOfParameters
      REAL(C_DOUBLE), INTENT(IN) :: parametersData(*)
      INTEGER(C_INT), VALUE, INTENT(IN) :: numberOfIntermediate
      REAL(C_DOUBLE), INTENT(OUT) :: intermediateData(*)
      INTEGER(C_INT), INTENT(OUT) :: err
      
    END SUBROUTINE SolverDAEExternal_Integrate
    
  END INTERFACE

  INTERFACE Solver_DynamicThetaSet
    MODULE PROCEDURE Solver_DynamicThetaSetDP0
    MODULE PROCEDURE Solver_DynamicThetaSetDP1
  END INTERFACE Solver_DynamicThetaSet

  INTERFACE Solver_LabelSet
    MODULE PROCEDURE Solver_LabelSetC
    MODULE PROCEDURE Solver_LabelSetVS
  END INTERFACE Solver_LabelSet

  PUBLIC CellMLEquations_CellMLAdd

  PUBLIC CellMLEquations_CreateFinish,CellMLEquations_CreateStart

  PUBLIC CellMLEquations_Destroy
  
  PUBLIC CellMLEquations_LinearityTypeSet

  PUBLIC CellMLEquations_TimeDependenceTypeSet

  PUBLIC CellMLEquations_TimeSet
  
  PUBLIC Solver_DAESolverTypeSet

  PUBLIC Solver_DAETimesSet,Solver_DAETimeStepSet
  
  PUBLIC Solver_DAEEulerSolverTypeSet

  PUBLIC Solver_DAECellMLRHSEvaluate
  
  PUBLIC Solver_Destroy
  
  PUBLIC Solver_DynamicDegreeSet

  PUBLIC Solver_DynamicLinearityTypeSet
  
  PUBLIC Solver_DynamicOrderSet

  PUBLIC Solver_DynamicSchemeSet

  PUBLIC Solver_DynamicRestartSet

  PUBLIC Solver_DynamicThetaSet

  PUBLIC Solver_DynamicALESet

  PUBLIC Solver_DynamicTimesSet

  PUBLIC SolverEquations_BoundaryConditionsCreateFinish,SolverEquations_BoundaryConditionsCreateStart

  PUBLIC SolverEquations_CreateFinish,SolverEquations_CreateStart

  PUBLIC SolverEquations_Destroy
  
  PUBLIC SolverEquations_EquationsSetAdd

  PUBLIC SolverEquations_InterfaceConditionAdd

  PUBLIC SolverEquations_LinearityTypeSet

  PUBLIC SolverEquations_SparsityTypeSet

  PUBLIC SolverEquations_SymmetryTypeSet

  PUBLIC SolverEquations_TimeDependenceTypeSet

  PUBLIC Solver_LabelSet
  
  PUBLIC Solver_LibraryTypeSet

  PUBLIC Solver_LinearTypeSet
  
  PUBLIC Solver_LinearDirectTypeSet

  PUBLIC Solver_MumpsSetIcntl,Solver_MumpsSetCntl

  PUBLIC Solver_LinearIterativeAbsoluteToleranceSet

  PUBLIC Solver_LinearIterativeDivergenceToleranceSet
  
  PUBLIC Solver_LinearIterativeGMRESRestartSet

  PUBLIC Solver_LinearIterativeMaximumIterationsSet

  PUBLIC Solver_LinearIterativePreconditionerTypeSet

  PUBLIC Solver_LinearIterativeRelativeToleranceSet
  
  PUBLIC Solver_LinearIterativeSolutionInitialiseTypeSet

  PUBLIC Solver_LinearIterativeTypeSet
  
  PUBLIC Solver_GeometricTransformationArbitraryPathSet,Solver_GeometricTransformationClear
  
  PUBLIC Solver_GeometricTransformationNumberOfLoadIncrementsSet
  
  PUBLIC Solver_GeometricTransformationScalingsSet
  
  PUBLIC Solver_GeometricTransformationFieldSet
  
  PUBLIC Solver_GeometricTransformationMatrixSet
  
  PUBLIC Solver_GeometricTransformationRotationSet,Solver_GeometricTransformationTranslationSet

  PUBLIC Solver_DynamicAssemble,Solver_StaticAssemble

  PUBLIC Solver_QuasiNewtonAbsoluteToleranceSet

  PUBLIC Solver_QuasiNewtonLineSearchMonitorOutputSet

  PUBLIC Solver_QuasiNewtonLinesearchTypeSet

  PUBLIC Solver_QuasiNewtonJacobianCalculationTypeSet
  
  PUBLIC Solver_QuasiNewtonConvergenceTestTypeSet

  PUBLIC Solver_QuasiNewtonLinesearchMaxStepSet

  PUBLIC Solver_QuasiNewtonLinesearchStepToleranceSet

  PUBLIC Solver_QuasiNewtonMaxNumberOfIterationsSet

  PUBLIC Solver_QuasiNewtonMaximumFunctionEvaluationsSet

  PUBLIC Solver_QuasiNewtonSolutionInitialiseTypeSet

  PUBLIC Solver_QuasiNewtonSolutionToleranceSet
  
  PUBLIC Solver_QuasiNewtonRelativeToleranceSet

  PUBLIC Solver_QuasiNewtonTrustregionDelta0Set

  PUBLIC Solver_QuasiNewtonTrustRegionToleranceSet

  PUBLIC Solver_QuasiNewtonTypeSet

  PUBLIC Solver_QuasiNewtonRestartNumberSet

  PUBLIC Solver_QuasiNewtonRestartTypeSet

  PUBLIC Solver_QuasiNewtonScaleTypeSet

  PUBLIC Solver_QuasiNewtonSolveTypeSet

  PUBLIC Solver_NewtonAbsoluteToleranceSet

  PUBLIC Solver_NewtonLinesearchMonitorOutputSet

  PUBLIC Solver_NewtonLinesearchAlphaSet

  PUBLIC Solver_NewtonLinesearchTypeSet

  PUBLIC Solver_NewtonJacobianCalculationTypeSet
  
  PUBLIC Solver_NewtonConvergenceTestTypeSet

  PUBLIC Solver_NewtonLinesearchMaxStepSet

  PUBLIC Solver_NewtonLinesearchStepToleranceSet

  PUBLIC Solver_NewtonMaxNumberOfIterationsSet

  PUBLIC Solver_NewtonMaximumFunctionEvaluationsSet

  PUBLIC Solver_NewtonSolutionInitialiseTypeSet

  PUBLIC Solver_NewtonSolutionToleranceSet
  
  PUBLIC Solver_NewtonRelativeToleranceSet

  PUBLIC Solver_NewtonTrustregionDelta0Set

  PUBLIC Solver_NewtonTrustregionToleranceSet

  PUBLIC Solver_NewtonTypeSet

  PUBLIC Solver_NonlinearDivergenceExit

  PUBLIC SolverNonlinear_Monitor

  PUBLIC Solver_NonlinearTypeSet

  PUBLIC SolverOptimiser_Monitor
  
  PUBLIC Solver_OutputTypeSet
  
  PUBLIC Solver_Solve
  
  PUBLIC Solver_SolverEquationsGet

  PUBLIC SolverDAE_TimeSteppingMonitor
  
  PUBLIC Solver_TypeSet

  PUBLIC Solver_VariablesDynamicFieldUpdate

  PUBLIC Solver_VariablesDynamicFieldPreviousValuesUpdate

  PUBLIC Solver_VariablesDynamicNonlinearUpdate

  PUBLIC Solver_VariablesFieldUpdate

 PUBLIC Solvers_CreateFinish,Solvers_CreateStart

  PUBLIC Solvers_Destroy

  PUBLIC Solvers_NumberOfSolversSet

  PUBLIC Solver_CellMLEvaluatorInitialise,Solver_CellMLEvaluatorFinalise

  PUBLIC Solver_NewtonCellMLEvaluatorCreate

  PUBLIC Solver_LinkedSolverAdd,Solver_LinkedSolverRemove
  
  PUBLIC Solver_SolutionUpdate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Adds a CellML environment to a solvers CellML equations. \see OpenCMISS::Iron::cmfe_CellMLEquations_CellMLAdd
  SUBROUTINE CellMLEquations_CellMLAdd(cellMLEquations,cellML,cellMLIndex,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to add the CellML environment to.
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to add
    INTEGER(INTG), INTENT(OUT) :: cellMLIndex !<On return, the index of the added CellML environment.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx
    TYPE(CellMLPtrType), ALLOCATABLE :: newCellMLEnvironments(:)
    TYPE(SolverType), POINTER :: solver
    
    ENTERS("CellMLEquations_CellMLAdd",err,error,*999)

    CALL CellMLEquations_AssertNotFinished(cellMLEquations,err,error,*999)
    CALL CellML_AssertIsFinished(cellML,err,error,*999)
    NULLIFY(solver)
    CALL CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*999)
    ALLOCATE(newCellMLEnvironments(cellMLEquations%numberOfCellMLEnvironments+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new CellML environments.",err,error,*999)
    DO cellmlIdx=1,cellMLEquations%numberOfCellMLEnvironments
      newCellMLEnvironments(cellmlIdx)%ptr=>cellMLEquations%cellMLEnvironments(cellmlIdx)%ptr
    ENDDO !cellmlIdx
    newCellMLEnvironments(cellMLEquations%numberOfCellMLEnvironments+1)%ptr=>cellML
    CALL MOVE_ALLOC(newCellMLEnvironments,cellMLEquations%cellMLEnvironments)
    cellMLEquations%numberOfCellMLEnvironments=cellMLEquations%numberOfCellMLEnvironments+1
    cellMLIndex=cellMLEquations%numberOfCellMLEnvironments
    
    EXITS("CellMLEquations_CellMLAdd")
    RETURN
999 IF(ALLOCATED(newCellMLEnvironments)) DEALLOCATE(newCellMLEnvironments)
    ERRORSEXITS("CellMLEquations_CellMLAdd",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_CellMLAdd
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating CellML equations
  SUBROUTINE CellMLEquations_CreateFinish(cellMLEquations,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver

    ENTERS("CellMLEquations_CreateFinish",err,error,*999)

    CALL CellMLEquations_AssertNotFinished(cellMLEquations,err,error,*999)
    NULLIFY(solver)
    CALL CellMLEquations_SolverGet(cellMLEquations,solver,err,error,*999)
    cellMLEquations%cellMLEquationsFinished=.TRUE.
        
    EXITS("CellMLEquations_CreateFinish")
    RETURN
999 ERRORSEXITS("CellMLEquations_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating CellML equations
  SUBROUTINE CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to start the creation of CellML equations on
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<On return, A pointer the CellML equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLEquations_CreateStart",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is already associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
#endif    
        
    CALL CellMLEquations_Initialise(solver,err,error,*999)
    cellMLEquations=>solver%cellMLEquations
        
    EXITS("CellMLEquations_CreateStart")
    RETURN
999 NULLIFY(cellMLEquations)
998 ERRORSEXITS("CellMLEquations_CreateStart",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_CreateStart
        
  !
  !================================================================================================================================
  !

  !>Destroys the CellML equations
  SUBROUTINE CellMLEquations_Destroy(cellMLEquations,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLEquations_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML equations is not associated.",err,error,*999)
    
    CALL CellMLEquations_Finalise(cellMLEquations,err,error,*999)
        
    EXITS("CellMLEquations_Destroy")
    RETURN
999 ERRORSEXITS("CellMLEquations_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_Destroy
        
  !
  !================================================================================================================================
  !

  !>Finalises the CellML equations and deallocates all memory.
  SUBROUTINE CellMLEquations_Finalise(cellMLEquations,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("CellMLEquations_Finalise",err,error,*999)

    IF(ASSOCIATED(cellMLEquations)) THEN
      IF(ALLOCATED(cellMLEquations%cellMLEnvironments)) DEALLOCATE(cellMLEquations%cellMLEnvironments)
      DEALLOCATE(cellMLEquations)
    ENDIF
        
    EXITS("CellMLEquations_Finalise")
    RETURN
999 ERRORSEXITS("CellMLEquations_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE CellMLEquations_Finalise
        
  !
  !================================================================================================================================
  !

  !>Initialises the CellML equations for a solver.
  SUBROUTINE CellMLEquations_Initialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the CellML equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("CellMLEquations_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%cellMLEquations)) CALL FlagError("CellML equations is already associated for this solver.",err,error,*998)
    
    ALLOCATE(solver%cellMLEquations,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate CellML equations.",err,error,*999)
    solver%cellMLEquations%solver=>solver
    solver%cellMLEquations%cellMLEquationsFinished=.FALSE.
    solver%cellMLEquations%linearity=CELLML_EQUATIONS_LINEAR
    solver%cellMLEquations%timeDependence=CELLML_EQUATIONS_STATIC
    solver%cellMLEquations%currentTime=0.0_DP
    solver%cellMLEquations%numberOfCellMLEnvironments=0
       
    EXITS("CellMLEquations_Initialise")
    RETURN
999 CALL CellMLEquations_Finalise(solver%cellMLEquations,dummyErr,dummyError,*998)
998 ERRORSEXITS("CellMLEquations_Initialise",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_Initialise
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for CellML equations \see OpenCMISS::Iron::cmfe_CellMLEquations_LinearityTypeSet
  SUBROUTINE CellMLEquations_LinearityTypeSet(cellMLEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to set the linearity type for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The type of linearity to be set \see SolverRoutines_CellMLEquationLinearityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellMLEquations_LinearityTypeSet",err,error,*999)

    CALL CellMLEquations_AssertNotFinished(cellMLEquations,err,error,*999)
   
    SELECT CASE(linearityType)
    CASE(CELLML_EQUATIONS_LINEAR)
      cellMLEquations%linearity=CELLML_EQUATIONS_LINEAR
    CASE(CELLML_EQUATIONS_NONLINEAR)
      cellMLEquations%linearity=CELLML_EQUATIONS_NONLINEAR
    CASE DEFAULT
      localError="The specified CellML equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CellMLEquations_LinearityTypeSet")
    RETURN
999 ERRORSEXITS("CellMLEquations_LinearityTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_LinearityTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for CellML equations \see OpenCMISS::Iron::cmfe_CellMLEquations_TimeDependenceTypeSet
  SUBROUTINE CellMLEquations_TimeDependenceTypeSet(cellMLEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to set the time dependence type for
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The type of time dependence to be set \see SolverRoutines_CellMLEquationTimeDependenceTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("CellMLEquations_TimeDependenceTypeSet",err,error,*999)

    CALL CellMLEquations_AssertNotFinished(cellMLEquations,err,error,*999)
   
    SELECT CASE(timeDependenceType)
    CASE(CELLML_EQUATIONS_STATIC)
      cellMLEquations%timeDependence=CELLML_EQUATIONS_STATIC
    CASE(CELLML_EQUATIONS_QUASISTATIC)
      cellMLEquations%timeDependence=CELLML_EQUATIONS_QUASISTATIC
    CASE(CELLML_EQUATIONS_DYNAMIC)
      cellMLEquations%timeDependence=CELLML_EQUATIONS_DYNAMIC
    CASE DEFAULT
      localError="The specified CellML equations time dependence type of "// &
        & TRIM(NumberToVString(timeDependenceType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("CellMLEquations_TimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("CellMLEquations_TimeDependenceTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_TimeDependenceTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the current time for CellML equations.
  SUBROUTINE CellMLEquations_TimeSet(cellMLEquations,time,err,error,*)

    !Argument variables
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations !<A pointer the CellML equations to set the time for.
    REAL(DP), INTENT(IN) :: time !<The time for the CellML equations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("CellMLEquations_TimeSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEquations)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif    

    cellMLEquations%currentTime=time
         
    EXITS("CellMLEquations_TimeSet")
    RETURN
999 ERRORSEXITS("CellMLEquations_TimeSet",err,error)
    RETURN 1
   
  END SUBROUTINE CellMLEquations_TimeSet

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a CellML evaluator solver 
  SUBROUTINE Solver_CellMLEvaluatorCreateFinish(cellMLEvaluatorSolver,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer to the CellML evaluator solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_CellMLEvaluatorCreateFinish",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluastor solver is not associated.",err,error,*999)
#endif    
    
    CALL FlagError("Not implemented.",err,error,*999)
       
    EXITS("Solver_CellMLEvaluatorCreateFinish")
    RETURN
999 ERRORSEXITS("Solver_CellMLEvaluatorCreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_CellMLEvaluatorCreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a CellML evaluator solver.
  SUBROUTINE Solver_CellMLEvaluatorFinalise(cellMLEvaluatorSolver,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_CellMLEvaluatorFinalise",err,error,*999)

    IF(ASSOCIATED(cellMLEvaluatorSolver)) THEN        
      DEALLOCATE(cellMLEvaluatorSolver)
    ENDIF
         
    EXITS("Solver_CellMLEvaluatorFinalise")
    RETURN
999 ERRORSEXITS("Solver_CellMLEvaluatorFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_CellMLEvaluatorFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a CellML evaluator solver for a solver.
  SUBROUTINE Solver_CellMLEvaluatorInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the CellML evaluator solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_CellMLEvaluatorInitialise",err,error,*998)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
#endif
    
    IF(ASSOCIATED(solver%cellMLEvaluatorSolver)) &
      & CALL FlagError("CellML evaluator solver is already associated for this solver.",err,error,*998)
    
    ALLOCATE(solver%cellMLEvaluatorSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver CellML evaluator solver.",err,error,*999)
    solver%cellMLEvaluatorSolver%solver=>solver
    solver%cellMLEvaluatorSolver%solverLibrary=SOLVER_CMISS_LIBRARY
       
    EXITS("Solver_CellMLEvaluatorInitialise")
    RETURN
999 CALL Solver_CellMLEvaluatorFinalise(solver%cellMLEvaluatorSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_CellMLEvaluatorInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_CellMLEvaluatorInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a CellML evaluator solver.
  SUBROUTINE SolverCellMLEvaluator_LibraryTypeSet(cellMLEvaluatorSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the CellML evaluator solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverCellMLEvaluator_LibraryTypeSet",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif    
    
    SELECT CASE(solverLibraryType)
    CASE(SOLVER_CMISS_LIBRARY)
      cellMLEvaluatorSolver%solverLibrary=SOLVER_CMISS_LIBRARY
    CASE DEFAULT
      localError="The specified solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
        & " is invalid for a CellML evaluator solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("SolverCellMLEvaluator_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverCellMLEvaluator_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverCellMLEvaluator_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Solve a CellML evaluator solver
  SUBROUTINE SolverCellMLEvaluator_Solve(cellMLEvaluatorSolver,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx,maximumNumberOfIntermediates,maximumNumberOfParameters,maximumNumberOfStates, &
      & numberOfCellMLEnvironments,onlyOneModelIndex,solverLibraryType,totalNumberOfDOFs
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP) :: time
    REAL(DP), POINTER :: intermediateData(:),parametersData(:),stateData(:)
    TYPE(CellMLType), POINTER :: cellML
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(FieldType), POINTER :: modelsField,stateField,parametersField,intermediateField
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverCellMLEvaluator_Solve",err,error,*999)

#ifdef WITH_PRECHECKS    
    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
#endif
    
    NULLIFY(solver)
    CALL SolverCellMLEvaluator_SolverGet(cellMLEvaluatorSolver,solver,err,error,*999)
    NULLIFY(cellMLEquations)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    CALL CellMLEquations_TimeGet(cellMLEquations,time,err,error,*999)
    CALL CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*999)
    DO cellmlIdx=1,numberOfCellMLEnvironments
      NULLIFY(cellML)
      CALL CellMLEquations_CellMLGet(cellMLEquations,cellMLIdx,cellML,err,error,*999)
      NULLIFY(cellMLModelsField)
      CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
      NULLIFY(modelsField)
      CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
      CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,totalNumberOfDOFs,err,error,*999)
      NULLIFY(modelsData)
      CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)

!!TODO: Maybe move this getting of fields earlier up the DAE solver chain? For now keep here.
                      
      !Make sure CellML fields have been updated to the current value of any mapped fields
      CALL CellML_FieldToCellMLUpdate(cellML,err,error,*999)
                      
      !Get the state information if this environment has any.
      NULLIFY(cellMLStateField)
      NULLIFY(stateField)
      NULLIFY(stateData)
      CALL CellML_CellMLStateFieldExists(cellML,cellMLStateField,err,error,*999)
      IF(ASSOCIATED(cellMLStateField)) THEN
        CALL CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*999)
        CALL Field_ParameterSetDataGet(stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,stateData,err,error,*999)
      ENDIF
                      
      !Get the parameters information if this environment has any.
      NULLIFY(cellMLParametersField)
      NULLIFY(parametersField)
      NULLIFY(parametersData)
      CALL CellML_CellMLParametersFieldExists(cellML,cellMLParametersField,err,error,*999)
      IF(ASSOCIATED(cellMLParametersField)) THEN
        CALL CellMLParametersField_ParametersFieldGet(cellMLParametersField,parametersField,err,error,*999)
        CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)
      ENDIF
                      
      !Get the intermediate information if this environment has any.
      NULLIFY(cellMLIntermediateField)
      NULLIFY(intermediateField)
      NULLIFY(intermediateData)
      CALL CellML_CellMLIntermediateFieldExists(cellML,cellMLIntermediateField,err,error,*999)
      IF(ASSOCIATED(cellMLIntermediateField)) THEN
        CALL CellMLIntermediateField_IntermediateFieldGet(cellMLIntermediateField,intermediateField,err,error,*999)
        CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData, &
          & err,error,*999)                            
      ENDIF

      !Solve these CellML equations
      CALL SolverCellMLEvaluator_LibraryTypeGet(cellMLEvaluatorSolver,solverLibraryType,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL CellMLModelsField_OnlyOneModelIndexGet(cellMLModelsField,onlyOneModelIndex,err,error,*999)
        CALL CellML_MaximumNumberOfIntermediateGet(cellML,maximumNumberOfIntermediates,err,error,*999)
        CALL CellML_MaximumNumberOfParametersGet(cellML,maximumNumberOfParameters,err,error,*999)
        CALL CellML_MaximumNumberOfStateGet(cellML,maximumNumberOfStates,err,error,*999)
        CALL SolverCellMLEvaluator_Evaluate(cellMLEvaluatorSolver,time,cellML,totalNumberOfDofs,onlyOneModelIndex, &
          & modelsData,maximumNumberOfStates,stateData,maximumNumberOfParameters,parametersData, &
          & maximumNumberOfIntermediates,intermediateData,err,error,*999)
      CASE DEFAULT
        localError="The CellML evaluator solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid or not implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT
                   
      !Restore field data
      CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      IF(ASSOCIATED(stateField)) CALL Field_ParameterSetDataRestore(stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & stateData,err,error,*999)
      IF(ASSOCIATED(parametersField)) CALL Field_ParameterSetDataRestore(parametersField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)
      IF(ASSOCIATED(intermediateField)) CALL Field_ParameterSetDataRestore(intermediateField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,intermediateData,err,error,*999)
                   
      !Make sure fields have been updated to the current value of any mapped CellML fields
      CALL CellML_CellMLToFieldUpdate(cellML,err,error,*999)
                    
    ENDDO !cellmlIdx
        
    EXITS("SolverCellMLEvaluator_Solve")
    RETURN
999 ERRORSEXITS("SolverCellMLEvaluator_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverCellMLEvaluator_Solve

  !
  !================================================================================================================================
  !
  
  !>Evaluate the CellML equations. 
  SUBROUTINE SolverCellMLEvaluator_Evaluate(cellMLEvaluatorSolver,time,cellML,N,onlyOneModelIndex,modelsData, &
    & maximumNumberOfStates,stateData,maximumNumberOfParameters,parametersData,maximumNumberOfIntermediates, &
    & intermediateData,err,error,*)

    !Argument variables
    TYPE(CellMLEvaluatorSolverType), POINTER :: cellMLEvaluatorSolver !<A pointer the CellML evaluator equation solver to evaluate
    REAL(DP), INTENT(IN) :: time !<The time for the CellML evaluate.
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to integrate the equations for.
    INTEGER(INTG), INTENT(IN) :: N !<The number of degrees-of-freedom
    INTEGER(INTG), INTENT(IN) :: onlyOneModelIndex !<If only one model is used in the models data the index of that model. 0 otherwise.
    INTEGER(INTG), POINTER :: modelsData(:) !<modelsData(dofIdx). The models data for the dofIdx'th dof.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfStates !<The maximum number of state variables per dof
    REAL(DP), POINTER :: stateData(:) !<stateData(stateIdx,dofIdx). The state data for the stateIdx'th state variable of the dofIdx'th dof. stateIdx varies from 1..numberOfStates.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfParameters !<The maximum number of parameter variables per dof.
    REAL(DP), POINTER :: parametersData(:) !<parametersData(parameterIdx,dofIdx). The parameters data for the parameterIdx'th parameter variable of the dofIdx'th dof. parameterIdx varies from 1..numberOfParameters.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIntermediates !<The maximum number of intermediate variables per dof.
    REAL(DP), POINTER :: intermediateData(:) !<intermediateData(intermediateIdx,dofIdx). The intermediate values data for the intermediateIdx'th intermediate variable of the dofIdx'th dof. intermediateIdx varies from 1..numberOfIntermediates
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,dofOrderType,intermediateEndDOF,intermediateIdx,intermediateStartDOF,modelIdx, &
      & numberOfIntermediates,numberOfParameters,numberOfStates,parameterEndDOF,parameterIdx,parameterStartDOF, &
      & stateEndDOF,stateIdx,stateStartDOF
    REAL(DP) :: intermediates(MAX(1,maximumNumberOfIntermediates)),parameters(MAX(1,maximumNumberOfParameters)), &
      & rates(MAX(1,maximumNumberOfStates)),states(MAX(1,maximumNumberOfStates))
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(FieldType), POINTER :: modelsField
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverCellMLEvaluator_Evaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(cellMLEvaluatorSolver)) CALL FlagError("CellML evaluator solver is not associated.",err,error,*999)
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    NULLIFY(modelsField)
    CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
    CALL Field_DOFOrderTypeGet(modelsField,FIELD_U_VARIABLE_TYPE,dofOrderType,err,error,*999)
    IF(dofOrderType==FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
      !Dof components are separated. Will need to copy data to temporary arrays.
      IF(onlyOneModelIndex/=0) THEN
        !We have CellML models on this rank
        IF(onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
          !Mulitple models
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN
              NULLIFY(cellMLModel)
              CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
              CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
              CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
              CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
              
              !Copy CellML data to temporary arrays
              DO stateIdx=1,numberOfStates
                states(stateIdx)=stateData((dofIdx-1)*N+stateIdx)
              ENDDO !stateIdx
              DO parameterIdx=1,numberOfParameters
                parameters(parameterIdx)=parametersData((dofIdx-1)*N+parameterIdx)
              ENDDO !parameterIdx

              ASSERT_WITH_CELLML()
              
#ifdef WITH_CELLML
              
              CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates,parameters)
              
#endif

              !Copy temporary data back to CellML arrays
              DO intermediateIdx=1,numberOfIntermediates
                intermediateData((dofIdx-1)*N+intermediateIdx)=intermediates(intermediateIdx)
              ENDDO !intermediateIdx
              DO stateIdx=1,numberOfStates
                stateData((dofIdx-1)*N+stateIdx)=states(stateIdx)
              ENDDO !stateIdx
              
            ENDIF !modelIdx                  
          ENDDO !dofIdx
        ELSE
          !One one model is used.
          NULLIFY(cellMLModel)
          CALL CellML_CellMLModelGet(cellML,onlyOneModelIndex,cellMLModel,err,error,*999)
          CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
          CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
          CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN
                      
              !Copy CellML data to temporary arrays
              DO stateIdx=1,numberOfStates
                states(stateIdx)=stateData((dofIdx-1)*N+stateIdx)
              ENDDO !stateIdx
              DO parameterIdx=1,numberOfParameters
                parameters(parameterIdx)=parametersData((dofIdx-1)*N+parameterIdx)
              ENDDO !parameterIdx

              ASSERT_WITH_CELLML()
              
#ifdef WITH_CELLML
              
              CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates,parameters)
              
#endif
                      
              !Copy temporary data back to CellML arrays
              DO intermediateIdx=1,numberOfIntermediates
                intermediateData((dofIdx-1)*N+intermediateIdx)=intermediates(intermediateIdx)
              ENDDO !intermediateIdx
              DO stateIdx=1,numberOfStates
                stateData((dofIdx-1)*N+stateIdx)=states(stateIdx)
              ENDDO !stateIdx
              
            ENDIF !modelIdx
          ENDDO !dofIdx
        ENDIF
      ENDIF
    ELSE
      !Dof components are continguous. Can pass data directly.
      IF(onlyOneModelIndex/=0) THEN
        !We have CellML models on this rank
        IF(onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
          !Mulitple models
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN
              NULLIFY(cellMLModel)
              CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
              CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
              CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
              CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
              !Call RHS. Note some models might not have state, rates, intermediate or parameter data so call accordingly
              !to avoid indexing in to null pointers
              IF(numberOfStates>0) THEN
                IF(numberOfIntermediates>0) THEN
                  IF(numberOfParameters>0) THEN
                    !We have state, intermediate and parameters in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1

                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML
                    
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, &
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF), &
                      & parametersData(parameterStartDOF:parameterEndDOF))

#endif                            
                    
                  ELSE
                    !We do not have parameters in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML                    
                    
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)
                    
#endif
                    
                  ENDIF
                ELSE
                  IF(numberOfParameters>0) THEN
                    !We do not have intermediates in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML                    
                                                 
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parametersData(parameterStartDOF:parameterEndDOF))

#endif                    
                            
                  ELSE
                    !We do not have intermediates or parameters in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML                    
                                                 
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parameters)

#endif                    
                    
                  ENDIF
                ENDIF
              ELSE
                IF(numberOfIntermediates>0) THEN
                  IF(numberOfParameters>0) THEN
                    !We do not have any states in the model
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML                                                                     
                            
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
                      & intermediateData(intermediateStartDOF:intermediateEndDOF),parametersData( &
                      & parameterStartDOF:parameterEndDOF))

#endif
                    
                  ELSE
                    !We do not have any states or parameters in the model
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML                                                                     
                            
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
                      & intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)

#endif                    
                    
                  ENDIF
                ELSE
                  CALL FlagError("Invalid CellML model - there are no states or intermediates.",err,error,*999)
                ENDIF
              ENDIF
            ENDIF  !modelIdx                
          ENDDO !dofIdx
        ELSE
          !One model is used.
          NULLIFY(cellMLModel)
          CALL CellML_CellMLModelGet(cellML,onlyOneModelIndex,cellMLModel,err,error,*999)
          CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
          CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
          CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
          !Call RHS. Note some models might not have state, rates, intermediate or parameter data so call accordingly
          !to avoid referencing null pointers
          IF(numberOfStates>0) THEN
            IF(numberOfIntermediates>0) THEN
              IF(numberOfParameters>0) THEN
                !We have states, intermediate and parameters for the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1                            
                    
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML
                    
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF),parametersData( &
                      & parameterStartDOF:parameterEndDOF))

#endif
                    
                  ENDIF !modelIdx
                ENDDO !dofIdx
              ELSE
                !We do not have parameters in the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                   
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML
                                        
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)

#endif                    
                    
                  ENDIF !modelIdx
                ENDDO !dofIdx                        
              ENDIF
            ELSE
              IF(numberOfParameters>0) THEN
                !We do not have any intermediates in the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN                            
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1
                   
                    ASSERT_WITH_CELLML()
                    
#ifdef WITH_CELLML
                                                                    
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parametersData(parameterStartDOF:parameterEndDOF))

#endif
                    
                  ENDIF !modelIdx
                ENDDO !dofIdx
              ELSE
                !We do not have any intermediates or parameters in the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN                    
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    
                    ASSERT_WITH_CELLML()
                   
#ifdef WITH_CELLML
                                                                                         
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,& 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parameters)
                    
#endif
                    
                  ENDIF !modelIdx
                ENDDO !dofIdx
              ENDIF
            ENDIF
          ELSE
            IF(numberOfIntermediates>0) THEN
              IF(numberOfParameters>0) THEN
                !We do not have any states in the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN
                            
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1
                    
                    ASSERT_WITH_CELLML()
                   
#ifdef WITH_CELLML
                                                                                                             
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
                      & intermediateData(intermediateStartDOF:intermediateEndDOF),parametersData( &
                      & parameterStartDOF:parameterEndDOF))

#endif                    
                    
                  ENDIF !modelIdx
                ENDDO !dofIdx
              ELSE
                !We do not have any states or parameters the model
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN                    
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    
                    ASSERT_WITH_CELLML()
                   
#ifdef WITH_CELLML
                                                                                                             
                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
                      & intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)

#endif
                    
                  ENDIF !modelIdx           
                ENDDO !dofIdx
              ENDIF
            ELSE
              CALL FlagError("Invalid CellML model - there are no states or intermediates.",err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
        
    EXITS("SolverCellMLEvaluator_Evaluate")
    RETURN
999 ERRORSEXITS("SolverCellMLEvaluator_Evaluate",err,error)
    RETURN 1
   
  END SUBROUTINE SolverCellMLEvaluator_Evaluate

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a solver 
  SUBROUTINE Solver_CreateFinish(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverIdx

    ENTERS("Solver_CreateFinish",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    
    !Set the finished flag. The final solver finish will be done once the solver equations have been finished.
    DO solverIdx=1,solver%numberOfLinkedSolvers
      solver%linkedSolvers(solverIdx)%ptr%solverFinished=.TRUE.
    ENDDO !solverIdx
    solver%solverFinished=.TRUE.
       
    EXITS("Solver_CreateFinish")
    RETURN
999 ERRORSEXITS("Solver_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise an Adams-Moulton differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_AdamsMoultonFinalise(adamsMoultonSolver,err,error,*)

    !Argument variables
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver !<A pointer the Adams-Moulton differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_AdamsMoultonFinalise",err,error,*999)

    IF(ASSOCIATED(adamsMoultonSolver)) THEN
      DEALLOCATE(adamsMoultonSolver)
    ENDIF
         
    EXITS("SolverDAE_AdamsMoultonFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_AdamsMoultonFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_AdamsMoultonFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an Adams-Moulton solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_AdamsMoultonInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise an Adams-Moulton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_AdamsMoultonInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%adamsMoultonSolver)) &
      & CALL FlagError("Adams-Moulton solver is already associated for this differential-algebraic equation solver.", &
      & err,error,*998)
    
    !Allocate the Adams-Moulton solver
    ALLOCATE(daeSolver%adamsMoultonSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Adams-Moulton solver.",err,error,*999)
    !Initialise
    daeSolver%adamsMoultonSolver%DAESolver=>daeSolver
    daeSolver%adamsMoultonSolver%solverLibrary=0
    !Defaults
         
    EXITS("SolverDAE_AdamsMoultonInitialise")
    RETURN
999 CALL SolverDAE_AdamsMoultonFinalise(daeSolver%adamsMoultonSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_AdamsMoultonInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_AdamsMoultonInitialise

  !
  !================================================================================================================================
  !

  !>Solve using an Adams-Moulton differential-algebraic equation solver.
  SUBROUTINE SolverDAEAdamsMoulton_Solve(adamsMoultonSolver,err,error,*)

    !Argument variables
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver !<A pointer the Adams-Moulton differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAEAdamsMoulton_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(adamsMoultonSolver)) &
      & CALL FlagError("Adams-Moulton differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
          
    EXITS("SolverDAEAdamsMoulton_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEAdamsMoulton_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEAdamsMoulton_Solve

 !
  !================================================================================================================================
  !

  !>Finishes the process of creating a differential-algebraic equation solver 
  SUBROUTINE SolverDAE_CreateFinish(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the differential-algebraic equation solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAE_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
       
    EXITS("SolverDAE_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverDAE_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a backward Euler differential-algebraic equation and deallocate all memory.
  SUBROUTINE SolverDAEEuler_BackwardFinalise(backwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver !<A pointer the backward Euler differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAEEuler_BackwardFinalise",err,error,*999)

    IF(ASSOCIATED(backwardEulerSolver)) THEN
      DEALLOCATE(backwardEulerSolver)
    ENDIF
         
    EXITS("SolverDAEEuler_BackwardFinalise")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_BackwardFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_BackwardFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a backward Euler solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAEEuler_BackwardInitialise(eulerDAESolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer the Euler differential-algebraic equation solver to initialise a backward Euler solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAEEuler_BackwardInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(eulerDAESolver)) &
      & CALL FlagError("Euler differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(eulerDAESolver%backwardEulerSolver)) &
      & CALL FlagError("Backward Euler solver is already associated for this Euler differential-algebraic equation solver.", &
      & err,error,*998)
      
    !Allocate the backward Euler solver
    ALLOCATE(eulerDAESolver%backwardEulerSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate backward Euler solver.",err,error,*999)
    !Initialise
    eulerDAESolver%backwardEulerSolver%eulerDAESolver=>eulerDAESolver
    eulerDAESolver%backwardEulerSolver%solverLibrary=0
    !Defaults
         
    EXITS("SolverDAEEuler_BackwardInitialise")
    RETURN
999 CALL SolverDAEEuler_BackwardFinalise(eulerDAESolver%backwardEulerSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAEEuler_BackwardInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_BackwardInitialise

  !
  !================================================================================================================================
  !

  !>Solve using a backward Euler differential-algebraic equation solver.
  SUBROUTINE SolverDAEEulerBackward_Solve(backwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver !<A pointer the backward Euler differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAEEulerBackward_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(backwardEulerSolver)) &
      & CALL FlagError("Backward Euler differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverDAEEulerBackward_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEEulerBackward_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEulerBackward_Solve

  !
  !================================================================================================================================
  !

  !>Finalise an Euler differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_EulerFinalise(eulerSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer the Euler differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_EulerFinalise",err,error,*999)

    IF(ASSOCIATED(eulerSolver)) THEN
      CALL SolverDAEEuler_ForwardFinalise(eulerSolver%forwardEulerSolver,err,error,*999)
      CALL SolverDAEEuler_BackwardFinalise(eulerSolver%backwardEulerSolver,err,error,*999)
      CALL SolverDAEEuler_ImprovedFinalise(eulerSolver%improvedEulerSolver,err,error,*999)      
      DEALLOCATE(eulerSolver)
    ENDIF
         
    EXITS("SolverDAE_EulerFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_EulerFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_EulerFinalise

  !
  !================================================================================================================================
  !

  !>Finalise a forward Euler differential-algebraic equation and deallocate all memory.
  SUBROUTINE SolverDAEEuler_ForwardFinalise(forwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<A pointer the forward Euler differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAEEuler_ForwardFinalise",err,error,*999)

    IF(ASSOCIATED(forwardEulerSolver)) THEN
      DEALLOCATE(forwardEulerSolver)
    ENDIF
         
    EXITS("SolverDAEEuler_ForwardFinalise")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_ForwardFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_ForwardFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a forward Euler solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAEEuler_ForwardInitialise(eulerDAESolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer the Euler differential-algebraic equation solver to initialise a forward Euler solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAEEuler_ForwardInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(eulerDAESolver)) &
      & CALL FlagError("Euler differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(eulerDAESolver%forwardEulerSolver)) &
      & CALL FlagError("Forward Euler solver is already associated for this Euler differential-algebraic equation solver.", &
      & err,error,*998)
     
    !Allocate the forward Euler solver
    ALLOCATE(eulerDAESolver%forwardEulerSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate forward Euler solver.",err,error,*999)
    !Initialise
    eulerDAESolver%forwardEulerSolver%eulerDAESolver=>eulerDAESolver
    eulerDAESolver%forwardEulerSolver%solverLibrary=SOLVER_CMISS_LIBRARY
    !Defaults
         
    EXITS("SolverDAEEuler_ForwardInitialise")
    RETURN
999 CALL SolverDAEEuler_ForwardFinalise(eulerDAESolver%forwardEulerSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAEEuler_ForwardInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_ForwardInitialise

  !
  !================================================================================================================================
  !

  !>Integrate using a forward Euler differential-algebraic equation solver.
  SUBROUTINE SolverDAEEulerForward_Integrate(forwardEulerSolver,cellML,N,startTime,endTime,timeIncrement, &
    & onlyOneModelIndex,modelsData,maximumNumberOfStates,stateData,maximumNumberOfParameters,parametersData, &
    & maximumNumberOfIntermediates,intermediateData,err,error,*)

    !Argument variables
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<A pointer the forward Euler differential-algebraic equation solver to integrate
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to integrate the equations for.
    INTEGER(INTG), INTENT(IN) :: N !<The number of degrees-of-freedom
    REAL(DP), INTENT(IN) :: startTime !<The start time for the integration
    REAL(DP), INTENT(IN) :: endTime !<The end time for the integration
    REAL(DP), INTENT(INOUT) :: timeIncrement !<The (initial) time increment for the integration
    INTEGER(INTG), INTENT(IN) :: onlyOneModelIndex !<If only one model is used in the models data the index of that model. 0 otherwise.
    INTEGER(INTG), POINTER :: modelsData(:) !<modelsData(dofIdx). The models data for the dofIdx'th dof.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfStates !<The maximum number of state variables per dof
    REAL(DP), POINTER :: stateData(:) !<stateData(stateIdx,dofIdx). The state data for the stateIdx'th state variable of the dofIdx'th dof. stateIdx varies from 1..numberOfStates.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfParameters !<The maximum number of parameter variables per dof.
    REAL(DP), POINTER :: parametersData(:) !<parametersData(parameterIdx,dofIdx). The parameters data for the parameterIdx'th parameter variable of the dofIdx'th dof. parameterIdx varies from 1..numberOfParameters.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIntermediates !<The maximum number of intermediate variables per dof.
    REAL(DP), POINTER :: intermediateData(:) !<intermediateData(intermediateIdx,dofIdx). The intermediate values data for the intermediateIdx'th intermediate variable of the dofIdx'th dof. intermediateIdx varies from 1.NUMBER_INTERMEDIATE    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,dofOrderType,intermediateEndDOF,intermediateIdx,intermediateStartDOF,modelIdx, &
      & numberOfIntermediates,numberOfParameters,numberOfStates,parameterEndDOF,parameterIdx,parameterStartDOF, &
      & stateEndDOF,stateIdx,stateStartDOF
    REAL(DP) :: intermediates(MAX(1,maximumNumberOfIntermediates)),parameters(MAX(1,maximumNumberOfParameters)), &
      & rates(MAX(1,maximumNumberOfStates)),states(MAX(1,maximumNumberOfStates)),time
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(FieldType), POINTER :: modelsField
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverDAEEulerForward_Integrate",err,error,*999)

    IF(.NOT.ASSOCIATED(forwardEulerSolver)) CALL FlagError("Forward Euler solver is not associated.",err,error,*999)
    
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    NULLIFY(modelsField)
    CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
    CALL Field_DOFOrderTypeGet(modelsField,FIELD_U_VARIABLE_TYPE,dofOrderType,err,error,*999)
    IF(dofOrderType==FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
      !Dof components are separated. Will need to copy data to temporary arrays.
      IF(onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !Mulitple models
        DO WHILE(time<=endTime)
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN
              NULLIFY(cellMLModel)
              CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
              CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
              CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
              CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)

              !Copy CellML data to temporary arrays
              DO stateIdx=1,numberOfStates
                states(stateIdx)=stateData((dofIdx-1)*N+stateIdx)
              ENDDO !stateIdx
              DO parameterIdx=1,numberOfParameters
                parameters(parameterIdx)=parametersData((dofIdx-1)*N+parameterIdx)
              ENDDO !parameterIdx

              ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

              CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates,parameters)

#endif

              !Copy temporary data back to CellML arrays
              DO intermediateIdx=1,numberOfIntermediates
                intermediateData((dofIdx-1)*N+intermediateIdx)=intermediates(intermediateIdx)
              ENDDO !intermediateIdx
              DO stateIdx=1,numberOfStates
                stateData((dofIdx-1)*N+stateIdx)=states(stateIdx)+timeIncrement*rates(stateIdx)
              ENDDO !stateIdx

            ENDIF !modelIdx                  
          ENDDO !dofIdx
          time=time+timeIncrement
        ENDDO !time              
      ELSE
        !One one model is used.
        NULLIFY(cellMLModel)
        CALL CellML_CellMLModelGet(cellML,onlyOneModelIndex,cellMLModel,err,error,*999)
        CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
        CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
        CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
        time=startTime
        DO WHILE(time<=endTime)
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN

              !Copy CellML data to temporary arrays
              DO stateIdx=1,numberOfStates
                states(stateIdx)=stateData((dofIdx-1)*N+stateIdx)
              ENDDO !stateIdx
              DO parameterIdx=1,numberOfParameters
                parameters(parameterIdx)=parametersData((dofIdx-1)*N+parameterIdx)
              ENDDO !parameterIdx

              ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

              CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates,parameters)

#endif

              !Copy temporary data back to CellML arrays
              DO intermediateIdx=1,numberOfIntermediates
                intermediateData((dofIdx-1)*N+intermediateIdx)=intermediates(intermediateIdx)
              ENDDO !intermediateIdx
              DO stateIdx=1,numberOfStates
                stateData((dofIdx-1)*N+stateIdx)=states(stateIdx)+timeIncrement*rates(stateIdx)
              ENDDO !stateIdx

            ENDIF !modelIdx
          ENDDO !dofIdx
          time=time+timeIncrement
        ENDDO !time
      ENDIF
    ELSE
      !Dof components are continguous. Can pass data directly.
      IF(onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
        !Mulitple models
        time=startTime
        DO WHILE(time<=endTime)
          DO dofIdx=1,N
            modelIdx=modelsData(dofIdx)
            IF(modelIdx>0) THEN
              NULLIFY(cellMLModel)
              CALL CellML_CellMLModelGet(cellML,modelIdx,cellMLModel,err,error,*999)
              CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
              CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
              CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
              !Call RHS. Note some models might not have state, rates, intermediate or parameter data so call accordingly
              !to avoid referencing null pointers
              IF(numberOfStates>0) THEN
                IF(numberOfIntermediates>0) THEN
                  IF(numberOfParameters>0) THEN
                    !We have states, intermediate and parameters for the model                     
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF: &
                      & stateEndDOF),rates,intermediateData(intermediateStartDOF:intermediateEndDOF), &
                      & parametersData(parameterStartDOF:parameterEndDOF))

#endif                    

                  ELSE
                    !We do not have parameters in the model                    
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF: &
                      & stateEndDOF),rates,intermediateData(intermediateStartDOF:intermediateEndDOF), &
                      & parameters)

#endif

                  ENDIF
                ELSE
                  IF(numberOfParameters>0) THEN
                    !We do not have intermediates in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF: &
                      & stateEndDOF),rates,intermediates,parametersData(parameterStartDOF:parameterEndDOF))

#endif                    

                  ELSE
                    !We do not have intermediates or parameters in the model
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF: &
                      & stateEndDOF),rates,intermediates,parameters)

#endif

                  ENDIF
                ENDIF
              ELSE
                CALL FlagError("Invalid CellML model for integration - there are no states.",err,error,*999)
              ENDIF
              stateData(stateStartDOF:stateEndDOF)=stateData(stateStartDOF:stateEndDOF)+timeIncrement*rates(1:numberOfStates)
            ENDIF
          ENDDO !dofIdx
          time=time+timeIncrement
        ENDDO !time              
      ELSE
        !One one model is used.
        NULLIFY(cellmlModel)
        CALL CellML_CellMLModelGet(cellML,onlyOneModelIndex,cellMLModel,err,error,*999)
        CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
        CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
        CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
        !Call RHS. Note some models might not have state, rates, intermediate or parameter data so call accordingly
        !to avoid referencing null pointers
        IF(numberOfStates>0) THEN
          IF(numberOfIntermediates>0) THEN
            IF(numberOfParameters>0) THEN
              !We have states, intermediate and parameters for the model                      
              time=startTime
              DO WHILE(time<=endTime)
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, &
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF),parametersData( &
                      & parameterStartDOF:parameterEndDOF))

#endif

                    stateData(stateStartDOF:stateEndDOF)=stateData(stateStartDOF:stateEndDOF)+ &
                      & timeIncrement*rates(1:numberOfStates)
                  ENDIF !modelIdx
                ENDDO !dofIdx
                time=time+timeIncrement
              ENDDO !time
            ELSE
              !We do not have parameters in the model
              time=startTime
              DO WHILE(time<=endTime)
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN                      
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    intermediateStartDOF=(dofIdx-1)*maximumNumberOfIntermediates+1
                    intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)

#endif                      

                    stateData(stateStartDOF:stateEndDOF)=stateData(stateStartDOF:stateEndDOF)+ &
                      & timeIncrement*rates(1:numberOfStates)
                  ENDIF !modelIdx
                ENDDO !dofIdx
                time=time+timeIncrement
              ENDDO !time
            ENDIF
          ELSE
            IF(numberOfParameters>0) THEN
              !We do not have intermediates in the model                      
              time=startTime
              DO WHILE(time<=endTime)
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN                      
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1
                    parameterStartDOF=(dofIdx-1)*maximumNumberOfParameters+1
                    parameterEndDOF=parameterStartDOF+numberOfParameters-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parametersData(parameterStartDOF:parameterEndDOF))

#endif                      

                    stateData(stateStartDOF:stateEndDOF)=stateData(stateStartDOF:stateEndDOF)+ &
                      & timeIncrement*rates(1:numberOfStates)
                  ENDIF !modelIdx
                ENDDO !dofIdx
                time=time+timeIncrement
              ENDDO !time
            ELSE
              !We do not have intermediates or parameters in the model
              time=startTime
              DO WHILE(time<=endTime)
                DO dofIdx=1,N
                  modelIdx=modelsData(dofIdx)
                  IF(modelIdx>0) THEN
                    stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
                    stateEndDOF=stateStartDOF+numberOfStates-1

                    ASSERT_WITH_CELLML()

#ifdef WITH_CELLML

                    CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time, & 
                      & stateData(stateStartDOF:stateEndDOF), &
                      & rates,intermediates,parameters)

#endif                      

                    stateData(stateStartDOF:stateEndDOF)=stateData(stateStartDOF:stateEndDOF)+ &
                      & timeIncrement*rates(1:numberOfStates)
                  ENDIF !modelIdx
                ENDDO !dofIdx
                time=time+timeIncrement
              ENDDO !time
            ENDIF
          ENDIF
        ELSE
          CALL FlagError("Invalid CellML model for integration - there are no states.",err,error,*999)
        ENDIF
      ENDIF
    ENDIF
        
    EXITS("SolverDAEEulerForward_Integrate")
    RETURN
999 ERRORSEXITS("SolverDAEEulerForward_Integrate",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEulerForward_Integrate

  !
  !================================================================================================================================
  !

  !>Solve using a forward Euler differential-algebraic equation solver.
  SUBROUTINE SolverDAEEulerForward_Solve(forwardEulerSolver,err,error,*)

    !Argument variables
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver !<A pointer the forward Euler differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx,maximumNumberOfIntermediates,maximumNumberOfParameters,maximumNumberOfStates, &
      & numberOfCellMLEnvironments,onlyOneModelIndex,totalNumberOfDOFs
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP), POINTER :: intermediateData(:),parametersData(:),stateData(:)
    TYPE(CellMLType), POINTER :: cellML
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(EulerDAESolverType), POINTER :: eulerSolver
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(FieldType), POINTER :: modelsField,stateField,parametersField,intermediateField
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverDAEEulerForward_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(forwardEulerSolver)) &
      & CALL FlagError("Forward Euler differential-algebraic equation solver is not associated.",err,error,*999)

    NULLIFY(eulerSolver)
    CALL SolverDAEEulerForward_EulerSolverGet(forwardEulerSolver,eulerSolver,err,error,*999)
    NULLIFY(daeSolver)
    CALL SolverDAEEuler_DAESolverGet(eulerSolver,daeSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverDAE_SolverGet(daeSolver,solver,err,error,*999)
    NULLIFY(cellMLEquations)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    CALL CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*999)
    DO cellmlIdx=1,numberOfCellMLEnvironments
      NULLIFY(cellML)
      CALL CellMLEquations_CellMLGet(cellMLEquations,cellMLIdx,cellML,err,error,*999)
      NULLIFY(cellMLModelsField)
      CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
      NULLIFY(modelsField)
      CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
      CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,totalNumberOfDOFs,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
 
!!TODO: Maybe move this getting of fields earlier up the DAE solver chain? For now keep here.
                      
      !Make sure CellML fields have been updated to the current value of any mapped fields
      CALL CellML_FieldToCellMLUpdate(cellML,err,error,*999)
                      
      !Get the state information if this environment has any.
      NULLIFY(cellMLStateField)
      NULLIFY(stateField)
      NULLIFY(stateData)
      CALL CellML_CellMLStateFieldExists(cellML,cellMLStateField,err,error,*999)
      IF(ASSOCIATED(cellMLStateField)) THEN
        CALL CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*999)
        CALL Field_ParameterSetDataGet(stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & stateData,err,error,*999)
      ENDIF
                      
      !Get the parameters information if this environment has any.
      NULLIFY(cellMLParametersField)
      NULLIFY(parametersField)
      NULLIFY(parametersData)
      CALL CellML_CellMLParametersFieldExists(cellML,cellMLParametersField,err,error,*999)
      IF(ASSOCIATED(cellMLParametersField)) THEN
        CALL CellMLParametersField_ParametersFieldGet(cellMLParametersField,parametersField,err,error,*999)
        CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & parametersData,err,error,*999)
      ENDIF
                      
      !Get the intermediate information if this environment has any.
      NULLIFY(cellMLIntermediateField)
      NULLIFY(intermediateField)
      NULLIFY(intermediateData)
      CALL CellML_CellMLIntermediateFieldExists(cellML,cellMLIntermediateField,err,error,*999)
      IF(ASSOCIATED(cellMLIntermediateField)) THEN
        CALL CellMLIntermediateField_IntermediateFieldGet(cellMLIntermediateField,intermediateField,err,error,*999)
        CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & intermediateData,err,error,*999)                            
      ENDIF

      CALL CellMLModelsField_OnlyOneModelIndexGet(cellMLModelsField,onlyOneModelIndex,err,error,*999)
      CALL CellML_MaximumNumberOfIntermediateGet(cellML,maximumNumberOfIntermediates,err,error,*999)
      CALL CellML_MaximumNumberOfParametersGet(cellML,maximumNumberOfParameters,err,error,*999)
      CALL CellML_MaximumNumberOfStateGet(cellML,maximumNumberOfStates,err,error,*999)
      
      !Integrate these CellML equations
      CALL SolverDAEEulerForward_Integrate(forwardEulerSolver,cellML,totalNumberOfDofs, &
        & daeSolver%startTime,daeSolver%endTime,daeSolver%initialStep, &
        & onlyOneModelIndex,modelsData,maximumNumberOfStates,stateData,maximumNumberOfParameters, &
        & parametersData,maximumNumberOfIntermediates,intermediateData,err,error,*999)
      
      !Restore field data
      CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      IF(ASSOCIATED(stateField)) CALL Field_ParameterSetDataRestore(stateField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,stateData,err,error,*999)                    
      IF(ASSOCIATED(parametersField)) CALL Field_ParameterSetDataRestore(parametersField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)                    
      IF(ASSOCIATED(intermediateField)) CALL Field_ParameterSetDataRestore(intermediateField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData,err,error,*999)
      
      !Make sure fields have been updated to the current value of any mapped CellML fields
      CALL CellML_CellMLToFieldUpdate(cellML,err,error,*999)
                      
    ENDDO !cellmlIdx
         
    EXITS("SolverDAEEulerForward_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEEulerForward_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEulerForward_Solve

  !
  !================================================================================================================================
  !

  !>Finalise an improved Euler differential-algebaic equation and deallocate all memory.
  SUBROUTINE SolverDAEEuler_ImprovedFinalise(improvedEulerSolver,err,error,*)

    !Argument variables
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver !<A pointer the improved Euler differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAEEuler_ImprovedFinalise",err,error,*999)

    IF(ASSOCIATED(improvedEulerSolver)) THEN
      DEALLOCATE(improvedEulerSolver)
    ENDIF
         
    EXITS("SolverDAEEuler_ImprovedFinalise")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_ImprovedFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_ImprovedFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an improved Euler solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAEEuler_ImprovedInitialise(eulerDAESolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerDAESolver !<A pointer the Euler differential-algebraic equation solver to initialise an improved Euler solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAEEuler_ImprovedInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(eulerDAESolver)) &
      & CALL FlagError("Euler differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(eulerDAESolver%improvedEulerSolver)) &
      & CALL FlagError("Improved Euler solver is already associated for this Euler differential-algebraic equation solver.", &
      & err,error,*998)
     
    !Allocate the improved Euler solver
    ALLOCATE(eulerDAESolver%improvedEulerSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate improved Euler solver.",err,error,*999)
    !Initialise
    eulerDAESolver%improvedEulerSolver%eulerDAESolver=>eulerDAESolver
    eulerDAESolver%improvedEulerSolver%solverLibrary=0
    !Defaults
         
    EXITS("SolverDAEEuler_ImprovedInitialise")
    RETURN
999 CALL SolverDAEEuler_ImprovedFinalise(eulerDAESolver%improvedEulerSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAEEuler_ImprovedInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_ImprovedInitialise

  !
  !================================================================================================================================
  !

  !>Solve using an improved Euler differential-algebraic equation solver.
  SUBROUTINE SolverDAEEulerImproved_Solve(improvedEulerSolver,err,error,*)

    !Argument variables
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver !<A pointer the improved Euler differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAEEulerImproved_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(improvedEulerSolver)) &
      & CALL FlagError("Improved Euler differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverDAEEulerImproved_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEEulerImproved_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEulerImproved_Solve

  !
  !================================================================================================================================
  !

  !>Initialise an Euler solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_EulerInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise an Euler solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_EulerInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%eulerSolver)) &
      & CALL FlagError("Euler solver is already associated for this differential-algebraic equation solver.",err,error,*998)
     
    !Allocate the Euler solver
    ALLOCATE(daeSolver%eulerSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Euler solver.",err,error,*999)
    !Initialise
    daeSolver%eulerSolver%DAESolver=>daeSolver
    NULLIFY(daeSolver%eulerSolver%forwardEulerSolver)
    NULLIFY(daeSolver%eulerSolver%backwardEulerSolver)
    NULLIFY(daeSolver%eulerSolver%improvedEulerSolver)
    !Default to a forward Euler solver
    CALL SolverDAEEuler_ForwardInitialise(daeSolver%eulerSolver,err,error,*999)
    daeSolver%eulerSolver%eulerType=SOLVER_DAE_EULER_FORWARD
         
    EXITS("SolverDAE_EulerInitialise")
    RETURN
999 CALL SolverDAE_EulerFinalise(daeSolver%eulerSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_EulerInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_EulerInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for an Euler differential-algebraic equation solver
  SUBROUTINE SolverDAEEuler_LibraryTypeSet(eulerSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer the Euler differential-algebraic equation solver to set the library type for
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the Euler differential-algebraic equation solver to set \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BackwardEulerDAESolverType), POINTER :: backwardEulerSolver
    TYPE(ForwardEulerDAESolverType), POINTER :: forwardEulerSolver
    TYPE(ImprovedEulerDAESolverType), POINTER :: improvedEulerSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAEEuler_LibraryTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(eulerSolver)) &
      & CALL FlagError("The Euler differential-algebraic equation solver is not associated.",err,error,*999)
    
    SELECT CASE(eulerSolver%eulerType)
    CASE(SOLVER_DAE_EULER_BACKWARD)
      NULLIFY(backwardEulerSolver)
      CALL SolverDAEEuler_BackwardEulerSolverGet(eulerSolver,backwardEulerSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a backward Euler DAE solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_EULER_FORWARD)
      NULLIFY(forwardEulerSolver)
      CALL SolverDAEEuler_ForwardEulerSolverGet(eulerSolver,forwardEulerSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        forwardEulerSolver%solverLibrary=SOLVER_CMISS_LIBRARY
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a forward Euler DAE solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_EULER_IMPROVED)
      NULLIFY(improvedEulerSolver)
      CALL SolverDAEEuler_ImprovedEulerSolverGet(eulerSolver,improvedEulerSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for an improved Euler DAE solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The Euler differential-algebraic equations solver type of "// &
        & TRIM(NumberToVString(eulerSolver%eulerType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("SolverDAEEuler_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Solve using an Euler differential-algebraic equation solver.
  SUBROUTINE SolverDAEEuler_Solve(eulerSolver,err,error,*)

    !Argument variables
    TYPE(EulerDAESolverType), POINTER :: eulerSolver !<A pointer the Euler differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAEEuler_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(eulerSolver)) &
      & CALL FlagError("Euler differential-algebraic equation solver is not associated.",err,error,*999)
    
    SELECT CASE(eulerSolver%eulerType)
    CASE(SOLVER_DAE_EULER_FORWARD)
      CALL SolverDAEEulerForward_Solve(eulerSolver%forwardEulerSolver,err,error,*999)
    CASE(SOLVER_DAE_EULER_BACKWARD)
      CALL SolverDAEEulerBackward_Solve(eulerSolver%backwardEulerSolver,err,error,*999)
    CASE(SOLVER_DAE_EULER_IMPROVED)
      CALL SolverDAEEulerImproved_Solve(eulerSolver%improvedEulerSolver,err,error,*999)
    CASE DEFAULT
      localError="The Euler differential-algebraic equation solver type of "// &
        & TRIM(NumberToVString(eulerSolver%eulerType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
          
    EXITS("SolverDAEEuler_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEEuler_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEEuler_Solve

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solve type for an Euler differential-algebraic equation solver.
  SUBROUTINE Solver_DAEEulerSolverTypeSet(solver,daeEulerType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the Euler differential equation solver to set type for 
    INTEGER(INTG), INTENT(IN) :: daeEulerType !<The type of Euler solver for the Euler differential-algebraic equation to set \see SolverRoutines_EulerDAESolverTypes,SolverRoutines.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(EulerDAESolverType), POINTER :: eulerSolver
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Solver_DAEEulerSolverTypeSet",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    CALL SolverDAE_AssertIsEuler(daeSolver,err,error,*999)
    NULLIFY(eulerSolver)
    CALL SolverDAE_EulerSolverGet(daeSolver,eulerSolver,err,error,*999)
    IF(daeEulerType/=eulerSolver%eulerType) THEN
      !Intialise the new Euler differential-algebraic equation solver type
      SELECT CASE(daeEulerType)
      CASE(SOLVER_DAE_EULER_BACKWARD)
        CALL SolverDAEEuler_BackwardInitialise(eulerSolver,err,error,*999)
      CASE(SOLVER_DAE_EULER_FORWARD)
        CALL SolverDAEEuler_ForwardInitialise(eulerSolver,err,error,*999)
      CASE(SOLVER_DAE_EULER_IMPROVED)
        CALL SolverDAEEuler_ImprovedInitialise(eulerSolver,err,error,*999)
      CASE DEFAULT
        localError="The specified Euler differential-algebraic equation solver type of "// &
          & TRIM(NumberToVString(daeEulerType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old Euler differential-algebraic equation solver type
      SELECT CASE(eulerSolver%eulerType)
      CASE(SOLVER_DAE_EULER_BACKWARD)
        CALL SolverDAEEuler_BackwardFinalise(eulerSolver%backwardEulerSolver,err,error,*999)
      CASE(SOLVER_DAE_EULER_FORWARD)
        CALL SolverDAEEuler_ForwardFinalise(eulerSolver%forwardEulerSolver,err,error,*999)
      CASE(SOLVER_DAE_EULER_IMPROVED)
        CALL SolverDAEEuler_ImprovedFinalise(eulerSolver%improvedEulerSolver,err,error,*999)
      CASE DEFAULT
        localError="The Euler differential-algebraic equation solver type of "// &
          & TRIM(NumberToVString(eulerSolver%eulerType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      eulerSolver%eulerType=daeEulerType
    ENDIF
        
    EXITS("Solver_DAEEulerSolverTypeSet")
    RETURN
999 ERRORSEXITS("Solver_DAEEulerSolverTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAEEulerSolverTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise a differential-algebraic equation solver and deallocate all memory
  SUBROUTINE Solver_DAEFinalise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_DAEFinalise",err,error,*999)

    IF(ASSOCIATED(daeSolver)) THEN
      CALL SolverDAE_EulerFinalise(daeSolver%eulerSolver,err,error,*999)
      CALL SolverDAE_CrankNicolsonFinalise(daeSolver%crankNicolsonSolver,err,error,*999)
      CALL SolverDAE_RungeKuttaFinalise(daeSolver%rungeKuttaSolver,err,error,*999)
      CALL SolverDAE_AdamsMoultonFinalise(daeSolver%adamsMoultonSolver,err,error,*999)
      CALL SolverDAE_BDFFinalise(daeSolver%bdfSolver,err,error,*999)
      CALL SolverDAE_RushLarsonFinalise(daeSolver%rushLarsonSolver,err,error,*999)
      CALL SolverDAE_ExternalFinalise(daeSolver%externalSolver,err,error,*999)
      DEALLOCATE(daeSolver)
    ENDIF
         
    EXITS("Solver_DAEFinalise")
    RETURN
999 ERRORSEXITS("Solver_DAEFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAEFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a differential-algebraic equation solver for a solver
  SUBROUTINE Solver_DAEInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the differential-algebraic equation solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_DAEInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%DAESolver)) &
      & CALL FlagError("Differential-algebraic equation solver is already associated for this solver.",err,error,*998)
      
    !Allocate the differential-algebraic equation solver
    ALLOCATE(solver%DAESolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver differential-algebraic equation solver.",err,error,*999)
    !Initialise
    solver%DAESolver%solver=>solver
    solver%DAESolver%daeType=0
    solver%DAESolver%daeSolveType=0
    solver%DAESolver%startTime=0.0_DP
    solver%DAESolver%endTime=0.1_DP
    solver%DAESolver%initialStep=0.1_DP
    NULLIFY(solver%DAESolver%eulerSolver)
    NULLIFY(solver%DAESolver%crankNicolsonSolver)
    NULLIFY(solver%DAESolver%rungeKuttaSolver)
    NULLIFY(solver%DAESolver%adamsMoultonSolver)
    NULLIFY(solver%DAESolver%bdfSolver)
    NULLIFY(solver%DAESolver%rushLarsonSolver)
    NULLIFY(solver%DAESolver%externalSolver)
    !Default to an Euler differential equation solver
    CALL SolverDAE_EulerInitialise(solver%DAESolver,err,error,*999)
    solver%DAESolver%daeSolveType=SOLVER_DAE_EULER
         
    EXITS("Solver_DAEInitialise")
    RETURN
999 CALL Solver_DAEFinalise(solver%DAESolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_DAEInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAEInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_LibraryTypeSet(daeSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to set the library type for
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the differential-algebraic equation solver to set \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(AdamsMoultonDAESolverType), POINTER :: adamsMoultonSolver
    TYPE(BDFDAESolverType), POINTER :: bdfSolver
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver
    TYPE(EulerDAESolverType), POINTER :: eulerSolver
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAE_LibraryTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("DAE solver is not associated.",err,error,*999)
    
    SELECT CASE(daeSolver%daeSolveType)
    CASE(SOLVER_DAE_EULER)
      NULLIFY(eulerSolver)
      CALL SolverDAE_EulerSolverGet(daeSolver,eulerSolver,err,error,*999)
      CALL SolverDAEEuler_LibraryTypeSet(eulerSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DAE_CRANK_NICOLSON)
      NULLIFY(crankNicolsonSolver)
      CALL SolverDAE_CrankNicolsonSolverGet(daeSolver,crankNicolsonSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_RUNGE_KUTTA)
      NULLIFY(rungeKuttaSolver)
      CALL SolverDAE_RungeKuttaSolverGet(daeSolver,rungeKuttaSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_ADAMS_MOULTON)
      NULLIFY(adamsMoultonSolver)
      CALL SolverDAE_AdamsMoultonSolverGet(daeSolver,adamsMoultonSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_BDF)
      NULLIFY(bdfSolver)
      CALL SolverDAE_BDFSolverGet(daeSolver,bdfSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        bdfSolver%solverLibrary = SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_RUSH_LARSON)
      NULLIFY(rushLarsonSolver)
      CALL SolverDAE_RushLarsonSolverGet(daeSolver,rushLarsonSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DAE_EXTERNAL)
      CALL FlagError("Can not set the library type for an external differential-algebraic equation solver is not associated.", &
        & err,error,*999)
    CASE DEFAULT
      localError="The differential-algebraic equations solver type of "// &
        & TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverDAE_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverDAE_LibraryTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise a BDF differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_BDFFinalise(bdfSolver,err,error,*)

    !Argument variables
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<A pointer the BDF differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_BDFFinalise",err,error,*999)

    IF(ASSOCIATED(bdfSolver)) THEN
      DEALLOCATE(bdfSolver)
    ENDIF
         
    EXITS("SolverDAE_BDFFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_BDFFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_BDFFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a BDF solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_BDFInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise a BDF solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_BDFInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%bdfSolver)) &
      & CALL FlagError("BDF solver is already associated for this differential-algebraic equation solver.",err,error,*998)
    
    !Allocate the BDF solver
    ALLOCATE(daeSolver%bdfSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate BDF solver.",err,error,*999)
    !Initialise
    daeSolver%bdfSolver%DAESolver=>daeSolver
    daeSolver%bdfSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    !Defaults
         
    EXITS("SolverDAE_BDFInitialise")
    RETURN
999 CALL SolverDAE_BDFFinalise(daeSolver%bdfSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_BDFInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_BDFInitialise
  !
  !================================================================================================================================
  !

  !>Finalise a CellML PETSc solver context.
  SUBROUTINE Solver_DAECellMLPETScContextFinalise(ctx,err,error,*)

    !Argument variables
    TYPE(CellMLPETScContextType), POINTER :: ctx !<A pointer the CellML-PETSc solver context to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_DAECellMLPETScContextFinalise",err,error,*999)

    IF(ASSOCIATED(ctx)) THEN
      IF(ASSOCIATED(ctx%rates)) DEALLOCATE(ctx%rates)
      IF(ALLOCATED(ctx%ratesIndices)) DEALLOCATE(ctx%ratesIndices)
      DEALLOCATE(ctx)
    ENDIF
         
    EXITS("Solver_DAECellMLPETScContextFinalise")
    RETURN
999 ERRORSEXITS("Solver_DAECellMLPETScContextFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAECellMLPETScContextFinalise


  !
  !================================================================================================================================
  !

  !>Initialise a CellML PETSc context
  SUBROUTINE Solver_DAECellMLPETScContextInitialise(ctx,err,error,*)

    !Argument variables
    TYPE(CellMLPETScContextType), INTENT(OUT), POINTER :: ctx !<A pointer to CellML PETSc context to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("Solver_DAECellMLPETScContextInitialise",err,error,*998)

    IF(ASSOCIATED(ctx)) CALL FlagError("Context is already associated.",err,error,*998)
    
    !Allocate the CTX
    ALLOCATE(ctx,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate context.",err,error,*999)
    !Initialise
    NULLIFY(ctx%solver)
    NULLIFY(ctx%cellml)
    NULLIFY(ctx%rates)
    ctx%dofIdx=0
         
    EXITS("Solver_DAECellMLPETScContextInitialise")
    RETURN
999 CALL Solver_DAECellMLPETScContextFinalise(ctx,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_DAECellMLPETScContextInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DAECellMLPETScContextInitialise
  
  !
  !================================================================================================================================
  !

  !>Set a CellML PETSc context
  SUBROUTINE Solver_DAECellMLPETScContextSet(ctx,solver,cellml,dofIdx,err,error,*)

    !Argument variables
    TYPE(CellMLPETScContextType), INTENT(IN), POINTER :: ctx !<A pointer to initialise a CELLML_PETSC_CONTEXT
    TYPE(SolverType), POINTER, INTENT(IN) :: solver !<A pointer to the solver to set to ctx
    TYPE(CellMLType), POINTER, INTENT(IN) :: cellml !<A pointer to the CellML environment to set to ctx
    INTEGER(INTG), INTENT(IN) :: dofIdx !<The DOF index of the cellml-petsc context
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: arrayIdx,dummyErr
    TYPE(VARYING_STRING) :: dummyError
   
    ENTERS("Solver_DAECellMLPETScContextSet",err,error,*998)

    IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Context is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(.NOT.ASSOCIATED(cellml)) CALL FlagError("CellML environment is not associated.",err,error,*998)
    
    !Set
    ctx%solver=>solver
    ctx%cellml=>cellml
    ctx%dofIdx=dofIdx
    ALLOCATE(ctx%rates(cellml%maximumNumberOfState),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate context rates.",err,error,*999)
    IF (.NOT.ALLOCATED(ctx%ratesIndices)) THEN
      ALLOCATE(ctx%ratesIndices(cellml%maximumNumberOfState),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate context rates indicies.",err,error,*999)
      ctx%ratesIndices=[(arrayIdx,arrayIdx=0,(cellml%maximumNumberOfState-1))]
    ENDIF
         
    EXITS("Solver_DAECellMLPETScContextSet")
    RETURN
999 CALL Solver_DAECellMLPETScContextFinalise(ctx,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_DAECellMLPETScContextSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DAECellMLPETScContextSet

  !
  !================================================================================================================================
  !
  
  !>Integrate using a BDF differential-algebraic equation solver.
  SUBROUTINE SolverDAEBDF_Integrate(bdfSolver,cellML,N,startTime,endTime,timeIncrement,onlyOneModelIndex,modelsData, &
    & maximumNumberOfStates,stateData,maximumNumberOfParameters,parametersData,maximumNumberOfIntermediates, &
    & intermediateData,err,error,*)

    !Argument variables
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<A pointer the BDF differential-algebraic equation solver to integrate
    TYPE(CellMLType), POINTER :: cellML !<A pointer to the CellML environment to integrate the equations for.
    INTEGER(INTG), INTENT(IN) :: N !<The number of degrees-of-freedom
    REAL(DP), INTENT(IN) :: startTime !<The start time for the integration
    REAL(DP), INTENT(IN) :: endTime !<The end time for the integration
    REAL(DP), INTENT(INOUT) :: timeIncrement !<The (initial) time increment for the integration
    INTEGER(INTG), INTENT(IN) :: onlyOneModelIndex !<If only one model is used in the models data the index of that model. 0 otherwise.
    INTEGER(INTG), POINTER, INTENT(IN) :: modelsData(:) !<modelsData(dofIdx). The models data for the dofIdx'th dof.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfStates !<The maximum number of state variables per dof
    REAL(DP), POINTER, INTENT (INOUT) :: stateData(:) !<stateData(stateIdx,dofIdx). The state data for the stateIdx'th state variable of the dofIdx'th dof. stateIdx varies from 1..numberOfStates.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfParameters !<The maximum number of parameter variables per dof.
    REAL(DP), POINTER, INTENT(INOUT) :: parametersData(:) !<parametersData(parameterIdx,dofIdx). The parameters data for the parameterIdx'th parameter variable of the dofIdx'th dof. parameterIdx varies from 1..numberOfParameters.
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIntermediates !<The maximum number of intermediate variables per dof.
    REAL(DP), POINTER, INTENT(INOUT) :: intermediateData(:) !<intermediateData(intermediateIdx,dofIdx). The intermediate values data for the intermediateIdx'th intermediate variable of the dofIdx'th dof. intermediateIdx varies from 1.NUMBER_INTERMEDIATE    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofIdx,dofOrderType,modelIdx,numberOfStates,stateEndDOF,stateIdx,stateStartDOF,arrayIdx
    INTEGER(INTG), ALLOCATABLE :: arrayIndices(:)
    REAL(DP) :: finalSolvedTime,timeStep
    REAL(DP), ALLOCATABLE  :: statesTemp(:),ratesTemp(:)
    TYPE(CellMLModelType), POINTER :: cellMLModel
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLPETScContextType), POINTER :: ctx !<The passed through context
    TYPE(FieldType), POINTER :: modelsField
    TYPE(VARYING_STRING) :: localError
    TYPE(PetscTSType) :: ts !<The PETSc ts type
    TYPE(PetscVecType) :: petscRates
    TYPE(PetscVecType) :: petscCurrentStates !<The initial and final states for the DAE
    EXTERNAL :: Problem_SolverDAECellMLRHSPetsc 
    
    ENTERS("SolverDAEBDF_Integrate",err,error,*999)

    NULLIFY(ctx)
    timeStep=endTime-startTime
    IF(.NOT.ASSOCIATED(bdfSolver)) CALL FlagError("BDF solver is not associated.",err,error,*999)
    
    NULLIFY(cellMLModelsField)
    CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
    NULLIFY(modelsField)
    CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
    CALL Field_DOFOrderTypeGet(modelsField,FIELD_U_VARIABLE_TYPE,dofOrderType,err,error,*999)
    SELECT CASE(bdfSolver%solverLibrary)   
    CASE(SOLVER_PETSC_LIBRARY)
      IF(dofOrderType==FIELD_SEPARATED_COMPONENT_DOF_ORDER) THEN
        CALL FlagError("Not implemented.",err,error,*999)      
      ELSE !dof component order is contiguous
        IF(onlyOneModelIndex==CELLML_MODELS_FIELD_NOT_CONSTANT) THEN
          CALL FlagError("Not implemented.",err,error,*999)
        ELSE !only one model
          NULLIFY(cellMLModel)
          CALL CellML_CellMLModelGet(cellML,onlyOneModelIndex,cellMLModel,err,error,*999)
          !determine no. of states in model and allocate necessary arrays
          CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
          ALLOCATE(statesTemp(0:numberOfStates-1),STAT=err)
          ALLOCATE(ratesTemp(0:numberOfStates-1),STAT=err)
          ALLOCATE(arrayIndices(0:numberOfStates-1),STAT=err)
          arrayIndices = [(arrayIdx,arrayIdx=0,(numberOfStates-1))]          
                  
          !initialize context for petsc solving.
          CALL Solver_DAECellMLPETScContextInitialise(ctx,err,error,*999)
          DO dofIdx=1,N     
            modelIdx = modelsData(dofIdx)
            IF(modelIdx>0) THEN !if model is assigned to dof
              !access the state field data
              stateStartDOF=(dofIdx-1)*maximumNumberOfStates+1
              stateEndDOF=stateStartDOF+numberOfStates-1
              DO stateIdx=1,numberOfStates
                statesTemp(stateIdx-1) = stateData(stateStartDOF+stateIdx-1)
              ENDDO !stateIdx
              
              !create PETSC states vector to initialize solver
              CALL PETSc_VecCreateSeq(PETSC_COMM_SELF,numberOfStates,petscCurrentStates,err,error,*999)
              !CALL PETSc_VecSetSizes(petscCurrentStates,PETSC_DECIDE,(numberOfStates),err,error,*999)
              !CALL PETSc_VecSetFromOptions(petscCurrentStates,err,error,*999)
                      
              !create PETSC rates vector to return values from evaluating rhs routine
              CALL PETSc_VecCreateSeq(PETSC_COMM_SELF,numberOfStates,petscRates,err,error,*999)
              !CALL PETSc_VecSetSizes(petscRates,PETSC_DECIDE,(numberOfStates),err,error,*999)
              !CALL PETSc_VecSetFromOptions(petscRates,err,error,*999)
              
              !Set up PETSC TS context for sundials BDF solver
              CALL PETSc_TSCreate(PETSC_COMM_SELF,ts,err,error,*999)
              CALL PETSc_TSSetProblemType(ts,PETSC_TS_NONLINEAR,err,error,*999)
              CALL PETSc_TSSetType(ts,PETSC_TS_SUNDIALS,err,error,*999)
              CALL PETSc_TSSundialsSetType(ts,PETSC_SUNDIALS_BDF,err,error,*999)
              CALL PETSc_TSSundialsSetTolerance(ts,0.0000001_DP,0.0000001_DP,err,error,*999)
              !set the initial solution to the current state
              CALL PETSc_VecSetValues(petscCurrentStates,(numberOfStates),arrayIndices,statesTemp,PETSC_INSERT_VALUES, &
                & err,error,*999)
              CALL PETSc_VecAssemblyBegin(petscCurrentStates,err,error,*999)
              CALL PETSc_VecAssemblyEnd(petscCurrentStates,err,error,*999)
              CALL PETSc_TSSetSolution(ts,petscCurrentStates,err,error,*999)
              
              !set up the time data
              CALL PETSc_TSSetInitialTimeStep(ts,startTime,timeIncrement,err,error,*999)
              CALL PETSc_TSSetDuration(ts,5000,endTime,err,error,*999)
              CALL PETSc_TSSetExactFinalTime(ts,.TRUE.,err,error,*999)
              
              IF(diagnostics1) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DAE start time = ",startTime,err,error,*999)
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DAE end time   = ",endTime,err,error,*999)
              ENDIF
                      
              !set rhs function and pass through the cellml model context 
              CALL Solver_DAECellMLPETScContextSet(ctx,bdfSolver%DAESolver%solver,cellML,dofIdx,err,error,*999)
              CALL PETSc_TSSetRHSFunction(ts,petscRates,Problem_SolverDAECellMLRHSPetsc,ctx,err,error,*999)
                      
              CALL PETSc_TSSolve(ts,petscCurrentStates,finalSolvedTime,err,error,*999)
              
              IF(diagnostics1) &
                & CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Final solved time = ",finalSolvedTime,err,error,*999)
              
              !update the states to new integrated values
              CALL PETSc_VecAssemblyBegin(petscCurrentStates,err,error,*999)
              CALL PETSc_VecAssemblyEnd(petscCurrentStates,err,error,*999)
              CALL PETSc_VecGetValues(petscCurrentStates,numberOfStates,arrayIndices,statesTemp,err,error,*999)
              
              DO stateIdx=1,numberOfStates                      
                stateData(stateStartDOF+stateIdx-1)=statesTemp(stateIdx-1)
              ENDDO
              CALL PETSc_TSFinalise(ts,err,error,*999) 
            ENDIF !modelIdx
            CALL PETSc_VecDestroy(petscCurrentStates,err,error,*999)
            CALL PETSc_VecDestroy(petscRates,err,error,*999)
          ENDDO !dofIdx
        ENDIF
      ENDIF !dof continguous
    CASE DEFAULT
      localError="The BDF solver library type of  "// &
        & TRIM(NumberToVString(bdfSolver%solverLibrary,"*",err,error))//" is not implemented."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverDAEBDF_Integrate")
    RETURN 
999 ERRORSEXITS("SolverDAEBDF_Integrate",err,error)
    RETURN 1

  END SUBROUTINE SolverDAEBDF_Integrate
  
  !
  !================================================================================================================================
  !

  !>Solve using a BDF differential-algebraic equation solver.
  SUBROUTINE SolverDAEBDF_Solve(bdfSolver,err,error,*)

    !Argument variables
    TYPE(BDFDAESolverType), POINTER :: bdfSolver !<A pointer the BDF differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx,maximumNumberOfIntermediate,maximumNumberOfParameters,maximumNumberOfState, &
      & numberOfCellMLEnvironments,onlyOneModelIndex,totalNumberOfDOFs
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP), POINTER :: intermediateData(:),parametersData(:),stateData(:)
    TYPE(CellMLType), POINTER :: cellML
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(FieldType), POINTER :: modelsField,stateField,parametersField,intermediateField
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAEBDF_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(bdfSolver)) CALL FlagError("BDF differential-algebraic equation solver is not associated.",err,error,*999)
    NULLIFY(daeSolver)
    CALL SolverDAEBDF_DAESolverGet(bdfSolver,daeSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverDAE_SolverGet(daeSolver,solver,err,error,*999)
    NULLIFY(cellMLEquations)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    CALL CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*999)
    DO cellmlIdx=1,numberOfCellMLEnvironments
      NULLIFY(cellML)
      CALL CellMLEquations_CellMLGet(cellMLEquations,cellMLIdx,cellML,err,error,*999)
      NULLIFY(cellMLModelsField)
      CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
      NULLIFY(modelsField)
      CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
      CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,totalNumberOfDOFs,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
 
!!TODO: Maybe move this getting of fields earlier up the DAE solver chain? For now keep here.
      
      !Make sure CellML fields have been updated to the current value of any mapped fields
      CALL CellML_FieldToCellMLUpdate(cellML,err,error,*999)
                      
      !Get the state information if this environment has any.
      NULLIFY(cellMLStateField)
      NULLIFY(stateField)
      NULLIFY(stateData)
      CALL CellML_CellMLStateFieldExists(cellML,cellMLStateField,err,error,*999)
      IF(ASSOCIATED(cellMLStateField)) THEN
        CALL CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*999)
        CALL Field_ParameterSetDataGet(stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,stateData,err,error,*999)
      ENDIF
                      
      !Get the parameters information if this environment has any.
      NULLIFY(cellMLParametersField)
      NULLIFY(parametersField)
      NULLIFY(parametersData)
      CALL CellML_CellMLParametersFieldExists(cellML,cellMLParametersField,err,error,*999)
      IF(ASSOCIATED(cellMLParametersField)) THEN
        CALL CellMLParametersField_ParametersFieldGet(cellMLParametersField,parametersField,err,error,*999)
        CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)
      ENDIF
                      
      !Get the intermediate information if this environment has any.
      NULLIFY(cellMLIntermediateField)
      NULLIFY(intermediateField)
      NULLIFY(intermediateData)
      CALL CellML_CellMLIntermediateFieldExists(cellML,cellMLIntermediateField,err,error,*999)
      IF(ASSOCIATED(cellMLIntermediateField)) THEN
        CALL CellMLIntermediateField_IntermediateFieldGet(cellMLIntermediateField,intermediateField,err,error,*999)
        CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData, &
          & err,error,*999)
      ENDIF
                    
      CALL CellMLModelsField_OnlyOneModelIndexGet(cellMLModelsField,onlyOneModelIndex,err,error,*999)
      CALL CellML_MaximumNumberOfIntermediateGet(cellML,maximumNumberOfIntermediate,err,error,*999)
      CALL CellML_MaximumNumberOfParametersGet(cellML,maximumNumberOfParameters,err,error,*999)
      CALL CellML_MaximumNumberOfStateGet(cellML,maximumNumberOfState,err,error,*999)
      
      !Integrate these CellML equations
      CALL SolverDAEBDF_Integrate(bdfSolver,cellML,totalNumberOfDofs,daeSolver%startTime,daeSolver%endTime,daeSolver%initialStep, &
        & onlyOneModelIndex,modelsData,maximumNumberOfState,stateData,maximumNumberOfParameters,parametersData, &
        & maximumNumberOfIntermediate,intermediateData,err,error,*999)
      
      !Restore field data
      CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
      IF(ASSOCIATED(stateField)) CALL Field_ParameterSetDataRestore(stateField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,stateData,err,error,*999)                    
      IF(ASSOCIATED(parametersField)) CALL Field_ParameterSetDataRestore(parametersField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)                    
      IF(ASSOCIATED(intermediateField)) CALL Field_ParameterSetDataRestore(intermediateField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData,err,error,*999)
      
      !Make sure fields have been updated to the current value of any mapped CellML fields
      CALL CellML_CellMLToFieldUpdate(cellML,err,error,*999)
                                        
    ENDDO !cellmlIdx
         
    EXITS("SolverDAEBDF_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEBDF_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEBDF_Solve
  
  !
  !================================================================================================================================
  !

  !>Finalise a Crank-Nicolson differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_CrankNicolsonFinalise(crankNicolsonSolver,err,error,*)

    !Argument variables
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver !<A pointer the Crank-Nicolson differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_CrankNicolsonFinalise",err,error,*999)

    IF(ASSOCIATED(crankNicolsonSolver)) THEN
      DEALLOCATE(crankNicolsonSolver)
    ENDIF
         
    EXITS("SolverDAE_CrankNicolsonFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_CrankNicolsonFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_CrankNicolsonFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Crank-Nicolson solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_CrankNicolsonInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise a Crank-Nicolson solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_CrankNicolsonInitialise",err,error,*998)

    IF(ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%crankNicolsonSolver)) &
      & CALL FlagError("Crank-Nicolson solver is already associated for this differential-algebraic equation solver.", &
      & err,error,*998)
    
    !Allocate the Crank-Nicholson solver
    ALLOCATE(daeSolver%crankNicolsonSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Crank-Nicolson solver.",err,error,*999)
    !Initialise
    daeSolver%crankNicolsonSolver%DAESolver=>daeSolver
    daeSolver%crankNicolsonSolver%solverLibrary=0
    !Defaults
        
    EXITS("SolverDAE_CrankNicolsonInitialise")
    RETURN
999 CALL SolverDAE_CrankNicolsonFinalise(daeSolver%crankNicolsonSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_CrankNicolsonInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_CrankNicolsonInitialise

  !
  !================================================================================================================================
  !

  !>Solve using a Crank-Nicolson differential-algebraic equation solver.
  SUBROUTINE SolverDAECrankNicolson_Solve(crankNicolsonSolver,err,error,*)

    !Argument variables
    TYPE(CrankNicolsonDAESolverType), POINTER :: crankNicolsonSolver !<A pointer the Crank-Nicolson differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAECrankNicolson_Solve",err,error,*999)

    IF(ASSOCIATED(crankNicolsonSolver)) &
      & CALL FlagError("Crank-Nicolson differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverDAECrankNicolson_Solve")
    RETURN
999 ERRORSEXITS("SolverDAECrankNicolson_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAECrankNicolson_Solve
  
  !
  !================================================================================================================================
  !

  !>Finalise an external differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_ExternalFinalise(externalSolver,err,error,*)

    !Argument variables
    TYPE(ExternalDAESolverType), POINTER :: externalSolver !<A pointer the external differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_ExternalFinalise",err,error,*999)

    IF(ASSOCIATED(externalSolver)) THEN
      DEALLOCATE(externalSolver)
    ENDIF
         
    EXITS("SolverDAE_ExternalFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_ExternalFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_ExternalFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an external solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_ExternalInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise an external solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_ExternalInitialise",err,error,*998)

    IF(ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%externalSolver)) &
      & CALL FlagError("External solver is already associated for this differential-algebraic equation solver.", &
      & err,error,*998)
    
    !Allocate the external solver
    ALLOCATE(daeSolver%externalSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate external solver.",err,error,*999)
    !Initialise
    daeSolver%externalSolver%DAESolver=>daeSolver
    !Defaults
         
    EXITS("SolverDAE_ExternalInitialise")
    RETURN
999 CALL SolverDAE_ExternalFinalise(daeSolver%externalSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_ExternalInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_ExternalInitialise

  !
  !================================================================================================================================
  !

  !>Solve using an external differential-algebraic equation solver.
  SUBROUTINE SolverDAEExternal_Solve(externalSolver,err,error,*)

    !Argument variables
    TYPE(ExternalDAESolverType), POINTER :: externalSolver !<A pointer the external differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx,maximumNumberOfIntermediates,maximumNumberOfParameters,maximumNumberOfStates, &
      & numberOfCellMLEnvironments,onlyOneModelIndex,totalNumberOfDOFs
    INTEGER(INTG), POINTER :: modelsData(:)
    REAL(DP), POINTER :: intermediateData(:),parametersData(:),stateData(:)
    TYPE(CellMLType), POINTER :: cellML
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(CellMLIntermediateFieldType), POINTER :: cellMLIntermediateField
    TYPE(CellMLModelsFieldType), POINTER :: cellMLModelsField
    TYPE(CellMLParametersFieldType), POINTER :: cellMLParametersField
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(FieldVariableType), POINTER :: modelsVariable
    TYPE(FieldType), POINTER :: modelsField,stateField,parametersField,intermediateField
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAEExternal_Solve",err,error,*999)

    IF(ASSOCIATED(externalSolver)) &
      & CALL FlagError("External Euler differential-algebraic equation solver is not associated.",err,error,*999)

    NULLIFY(daeSolver)
    CALL SolverDAEExternal_DAESolverGet(externalSolver,daeSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverDAE_SolverGet(daeSolver,solver,err,error,*999)
    NULLIFY(cellMLEquations)
    CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
    CALL CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*999)
    DO cellmlIdx=1,numberOfCellMLEnvironments
      NULLIFY(cellML)
      CALL CellMLEquations_CellMLGet(cellMLEquations,cellmlIdx,cellML,err,error,*999)
      NULLIFY(cellMLModelsField)
      CALL CellML_CellMLModelsFieldGet(cellML,cellMLModelsField,err,error,*999)
      NULLIFY(modelsField)
      CALL CellMLModelsField_ModelsFieldGet(cellMLModelsField,modelsField,err,error,*999)
      NULLIFY(modelsVariable)
      CALL Field_VariableGet(modelsField,FIELD_U_VARIABLE_TYPE,modelsVariable,err,error,*999)
      CALL FieldVariable_TotalNumberOfDOFsGet(modelsVariable,totalNumberOfDOFs,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)
       
      !Make sure CellML fields have been updated to the current value of any mapped fields
      CALL CellML_FieldToCellMLUpdate(cellML,err,error,*999)
                    
      !Get the state information if this environment has any.
      NULLIFY(cellMLStateField)
      NULLIFY(stateField)
      NULLIFY(stateData)
      CALL CellML_CellMLStateFieldExists(cellML,cellMLStateField,err,error,*999)
      IF(ASSOCIATED(cellMLStateField)) THEN
        CALL CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*999)
        CALL Field_ParameterSetDataGet(stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,stateData,err,error,*999)
      ENDIF
                    
      !Get the parameters information if this environment has any.
      NULLIFY(cellMLParametersField)
      NULLIFY(parametersField)
      NULLIFY(parametersData)
      CALL CellML_CellMLParametersFieldExists(cellML,cellMLParametersField,err,error,*999)
      IF(ASSOCIATED(cellMLParametersField)) THEN
        CALL CellMLParametersField_ParametersFieldGet(cellMLParametersField,parametersField,err,error,*999)
        CALL Field_ParameterSetDataGet(parametersField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & parametersData,err,error,*999)
      ENDIF
                    
      !Get the intermediate information if this environment has any.
      NULLIFY(cellMLIntermediateField)
      NULLIFY(intermediateField)
      NULLIFY(intermediateData)
      CALL CellML_CellMLIntermediateFieldExists(cellML,cellMLIntermediateField,err,error,*999)
      IF(ASSOCIATED(cellMLIntermediateField)) THEN
        CALL CellMLIntermediateField_IntermediateFieldGet(cellMLIntermediateField,intermediateField,err,error,*999)
        CALL Field_ParameterSetDataGet(intermediateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & intermediateData,err,error,*999)                            
      ENDIF

      CALL CellMLModelsField_OnlyOneModelIndexGet(cellMLModelsField,onlyOneModelIndex,err,error,*999)
      CALL CellML_MaximumNumberOfIntermediateGet(cellML,maximumNumberOfIntermediates,err,error,*999)
      CALL CellML_MaximumNumberOfParametersGet(cellML,maximumNumberOfParameters,err,error,*999)
      CALL CellML_MaximumNumberOfStateGet(cellML,maximumNumberOfStates,err,error,*999)
        
      !Call the external solver to integrate these CellML equations
      CALL SolverDAEExternal_Integrate(totalNumberOfDofs,daeSolver%startTime,daeSolver%endTime,daeSolver%initialStep, &
        & onlyOneModelIndex,modelsData,maximumNumberOfStates,stateData,maximumNumberOfParameters,parametersData, &
        & maximumNumberOfIntermediates,intermediateData,err)
      IF(err/=0) CALL FlagError("Error from external solver integrate.",err,error,*999)
      
      !Restore field data
      CALL FieldVariable_ParameterSetDataRestore(modelsVariable,FIELD_VALUES_SET_TYPE,modelsData,err,error,*999)                    
      IF(ASSOCIATED(stateField)) CALL Field_ParameterSetDataRestore(stateField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,stateData,err,error,*999)                    
      IF(ASSOCIATED(parametersField)) CALL Field_ParameterSetDataRestore(parametersField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,parametersData,err,error,*999)                    
      IF(ASSOCIATED(intermediateField)) CALL Field_ParameterSetDataRestore(intermediateField, &
        & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,intermediateData,err,error,*999)                    
      
      !Make sure fields have been updated to the current value of any mapped CellML fields
      CALL CellML_CellMLToFieldUpdate(cellML,err,error,*999)
      
    ENDDO !cellmlIdx
        
    EXITS("SolverDAEExternal_Solve")
    RETURN
999 ERRORSEXITS("SolverDAEExternal_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAEExternal_Solve
  
  !
  !================================================================================================================================
  !

  !>Integrate using a forward Euler differential-algebraic equation solver.
  SUBROUTINE Solver_DAECellMLRHSEvaluate(cellMLModel,time,stateStartIdx,stateDataOffset,stateData,parameterStartIdx, &
    & parameterDataOffset,parameterData,intermediateStartIdx,intermediateDataOffset,intermediateData,rateStartIdx, &
    & rateDataOffset,rateData,err,error,*)

    !Argument variables
    TYPE(CellMLModelType), POINTER :: cellMLModel !<The CellML model to evaluate
    REAL(DP), INTENT(IN) :: time !<The time to evaluate the CellML model at
    INTEGER(INTG), INTENT(IN) :: stateStartIdx !<The state start data offset.
    INTEGER(INTG), INTENT(IN) :: stateDataOffset !<The offset to the next state data
    REAL(DP), POINTER :: stateData(:) !<A pointer to the state data
    INTEGER(INTG), INTENT(IN) :: parameterStartIdx !<The parameter start data offset.
    INTEGER(INTG), INTENT(IN) :: parameterDataOffset !<The offset to the next parameters data
    REAL(DP), POINTER :: parameterData(:) !<A pointer to the parameters data
    INTEGER(INTG), INTENT(IN) :: intermediateStartIdx !<The intermediate start data offset.
    INTEGER(INTG), INTENT(IN) :: intermediateDataOffset !<The offset to the next intermediate data
    REAL(DP), POINTER :: intermediateData(:) !<A pointer to the intermediate data
    INTEGER(INTG), INTENT(IN) :: rateStartIdx !<The rate start data offset.
    INTEGER(INTG), INTENT(IN) :: rateDataOffset !<The offset to the next rates data
    REAL(DP), POINTER :: rateData(:) !<On exit, the rate data.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: intermediateIdx,intermediateEndDOF,intermediateStartDOF,numberOfIntermediates,numberOfParameters, &
      & numberOfStates,parameterIdx,parameterEndDOF,parameterStartDOF,rateIdx,rateEndDOF,rateStartDOF,stateIdx,stateEndDOF, &
      & stateStartDOF
    REAL(DP) :: intermediates(MAX(1,intermediateDataOffset)),parameters(MAX(1,parameterDataOffset)),rates(MAX(1,rateDataOffset)), &
      & states(MAX(1,stateDataOffset))
    
    ENTERS("Solver_DAECellMLRHSEvaluate",err,error,*999)

    CALL CellMLModel_NumberOfIntermediateGet(cellMLModel,numberOfIntermediates,err,error,*999)
    CALL CellMLModel_NumberOfParametersGet(cellMLModel,numberOfParameters,err,error,*999)
    CALL CellMLModel_NumberOfStateGet(cellMLModel,numberOfStates,err,error,*999)
    IF(numberOfStates>0) THEN
      IF(.NOT.ASSOCIATED(stateData)) CALL FlagError("State data is not associated.",err,error,*999)
      IF(.NOT.ASSOCIATED(rateData)) CALL FlagError("Rate data is not associated.",err,error,*999)
    ENDIF
    IF(numberOfParameters>0) THEN
      IF(.NOT.ASSOCIATED(parameterData)) CALL FlagError("Parameter data is not associated.",err,error,*999)
    ENDIF
    IF(numberOfIntermediates>0) THEN
      IF(.NOT.ASSOCIATED(intermediateData)) CALL FlagError("Intermediate data is not associated.",err,error,*999)
    ENDIF
    IF(stateDataOffset>1.OR.numberOfStates==0) THEN
      !State data is not contiguous or there are no states
        
      !Copy state data to temporary array
      DO stateIdx=1,numberOfStates
        states(stateIdx)=stateData((stateStartIdx-1)*stateDataOffset+stateIdx)
      ENDDO !stateIdx
      
      IF(parameterDataOffset>1.OR.numberOfParameters==0) THEN
        !Parameter data is not contiguous or there are no parameters
        
        !Copy parameter data to temporary array
        DO parameterIdx=1,numberOfParameters
          parameters(parameterIdx)=parameterData((parameterStartIdx-1)*parameterDataOffset+parameterIdx)
        ENDDO !parameterIdx
        
        IF(intermediateDataOffset>1.OR.numberOfIntermediates==0) THEN
          !Intermediate data is not contiguous or there are no intermediates
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates,parameters)

#endif            
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rateData(rateStartDOF:rateEndDOF), &
              & intermediates,parameters)

#endif            
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
          ENDIF
            
        ELSE
          !Intermediate data is contiguous
            
          intermediateStartDOF=(intermediateStartIdx-1)*intermediateDataOffset+1
          intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
              & intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)

#endif            
              
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rateData(rateStartDOF:rateEndDOF), &
              & intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)          

#endif
            
          ENDIF
        ENDIF
      ELSE
        !Parameters data is contiguous
        
        parameterStartDOF=(parameterStartIdx-1)*parameterDataOffset+1
        parameterEndDOF=parameterStartDOF+numberOfParameters-1
        
        IF(intermediateDataOffset>1.OR.numberOfIntermediates==0) THEN
          !Intermediate data is not contiguous or there are no intermediates
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates,intermediates, &
              & parameters(parameterStartDOF:parameterEndDOF))

#endif            
              
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
                
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rateData(rateStartDOF:rateEndDOF), &
              & intermediates,parameters(parameterStartDOF:parameterEndDOF))

#endif            
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
          ENDIF
          
        ELSE
          !Intermediate data is contiguous
          
          intermediateStartDOF=(intermediateStartIdx-1)*intermediateDataOffset+1
          intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rates, &
              & intermediateData(intermediateStartDOF:intermediateEndDOF), &
              & parameters(parameterStartDOF:parameterEndDOF))

#endif            
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,states,rateData(rateStartDOF:rateEndDOF), &
              & intermediateData(intermediateStartDOF:intermediateEndDOF), &
              & parameters(parameterStartDOF:parameterEndDOF))

#endif            
              
          ENDIF
        ENDIF
      ENDIF
    ELSE
      !State data is contiguous
      
      stateStartDOF=(stateStartIdx-1)*stateDataOffset+1
      stateEndDOF=stateStartDOF+numberOfStates-1
      
      IF(parameterDataOffset>1.OR.numberOfParameters==0) THEN
        !Parameter data is not contiguous or there are no parameters
        
        !Copy parameter data to temporary array
        DO parameterIdx=1,numberOfParameters
          parameters(parameterIdx)=parameterData((parameterStartIdx-1)*parameterDataOffset+parameterIdx)
        ENDDO !parameterIdx
        
        IF(intermediateDataOffset>1.OR.numberOfIntermediates==0) THEN
          !Intermediate data is not contiguous or there are no intermediates
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rates,intermediates,parameters)

#endif            
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rateData(rateStartDOF:rateEndDOF),intermediates,parameters)          

#endif
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
          ENDIF
            
        ELSE
          !Intermediate data is contiguous
          
          intermediateStartDOF=(intermediateStartIdx-1)*intermediateDataOffset+1
          intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF),rates, &
              & intermediateData(intermediateStartDOF:intermediateEndDOF),parameters)          

#endif
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rateData(rateStartDOF:rateEndDOF),intermediateData(intermediateStartDOF:intermediateEndDOF), &
              & parameters)

#endif            
            
          ENDIF
        ENDIF
      ELSE
        !Parameters data is contiguous
        
        parameterStartDOF=(parameterStartIdx-1)*parameterDataOffset+1
        parameterEndDOF=parameterStartDOF+numberOfParameters-1
        
        IF(intermediateDataOffset>1.OR.numberOfIntermediates==0) THEN
          !Intermediate data is not contiguous or there are no intermediates
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rates,intermediates,parameters(parameterStartDOF:parameterEndDOF))
            
#endif
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rateData(rateStartDOF:rateEndDOF),intermediates,parameters(parameterStartDOF:parameterEndDOF))

#endif            
            
            !Copy intermediate data from temporary array
            DO intermediateIdx=1,numberOfIntermediates
              intermediateData((intermediateStartIdx-1)*intermediateDataOffset+intermediateIdx)=intermediates(intermediateIdx)
            ENDDO !intermediateIdx
              
          ENDIF
          
        ELSE
          !Intermediate data is contiguous
          
          intermediateStartDOF=(intermediateStartIdx-1)*intermediateDataOffset+1
          intermediateEndDOF=intermediateStartDOF+numberOfIntermediates-1
          
          IF(rateDataOffset>1.OR.numberOfStates==0) THEN
            !Rates data is not contiguous or there are no rates

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rates,intermediateData(intermediateStartDOF:intermediateEndDOF), &
              & parameters(parameterStartDOF:parameterEndDOF))

#endif            
            
            !Copy rate data from temporary array
            DO rateIdx=1,numberOfStates
              rateData((rateStartIdx-1)*rateDataOffset+rateIdx)=rates(rateIdx)
            ENDDO !rateIdx
            
          ELSE
            !Rates data is contiguous
            
            rateStartDOF=(rateStartIdx-1)*rateDataOffset+1
            rateEndDOF=rateStartDOF+numberOfStates-1

            ASSERT_WITH_CELLML()

#ifdef WITH_CELLML
            
            CALL CELLML_MODEL_DEFINITION_CALL_RHS_ROUTINE(cellMLModel%ptr,time,stateData(stateStartDOF:stateEndDOF), &
              & rateData(rateStartDOF:rateEndDOF),intermediateData(intermediateStartDOF:intermediateEndDOF), &
              & parameters(parameterStartDOF:parameterEndDOF))

#endif
            
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
    EXITS("Solver_DAECellMLRHSEvaluate")
    RETURN
999 ERRORSEXITS("Solver_DAECellMLRHSEvaluate",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAECellMLRHSEvaluate

  !
  !================================================================================================================================
  !

  !>Finalise a Runge-Kutta differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_RungeKuttaFinalise(rungeKuttaSolver,err,error,*)

    !Argument variables
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver !<A pointer the Runge-Kutta differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_RungeKuttaFinalise",err,error,*999)

    IF(ASSOCIATED(rungeKuttaSolver)) THEN
      DEALLOCATE(rungeKuttaSolver)
    ENDIF
         
    EXITS("SolverDAE_RungeKuttaFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_RungeKuttaFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_RungeKuttaFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Runge-Kutta solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_RungeKuttaInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise a Runge-Kutta solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_RungeKuttaInitialise",err,error,*998)

    IF(ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%rungeKuttaSolver)) &
      & CALL FlagError("Runge-Kutta solver is already associated for this differential-algebraic equation solver.",err,error,*998)
      
    !Allocate the Runge-Kutta solver
    ALLOCATE(daeSolver%rungeKuttaSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Runge-Kutta solver.",err,error,*999)
    !Initialise
    daeSolver%rungeKuttaSolver%DAESolver=>daeSolver
    daeSolver%rungeKuttaSolver%solverLibrary=0
    !Defaults
         
    EXITS("SolverDAE_RungeKuttaInitialise")
    RETURN
999 CALL SolverDAE_RungeKuttaFinalise(daeSolver%rungeKuttaSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_RungeKuttaInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_RungeKuttaInitialise

  !
  !================================================================================================================================
  !

  !>Solve using a Runge-Kutta differential-algebraic equation solver.
  SUBROUTINE SolverDAERungeKutta_Solve(rungeKuttaSolver,err,error,*)

    !Argument variables
    TYPE(RungeKuttaDAESolverType), POINTER :: rungeKuttaSolver !<A pointer the Runge-Kutta differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAERungeKutta_Solve",err,error,*999)

    IF(ASSOCIATED(rungeKuttaSolver)) &
      & CALL FlagError("Runge-Kutta differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverDAERungeKutta_Solve")
    RETURN
999 ERRORSEXITS("SolverDAERungeKutta_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAERungeKutta_Solve
  
  !
  !================================================================================================================================
  !

  !>Finalise a Rush-Larson differential-algebraic equation solver and deallocate all memory.
  SUBROUTINE SolverDAE_RushLarsonFinalise(rushLarsonSolver,err,error,*)

    !Argument variables
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver !<A pointer the Rush-Larson differential-algebraic equation solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
     
    ENTERS("SolverDAE_RushLarsonFinalise",err,error,*999)

    IF(ASSOCIATED(rushLarsonSolver)) THEN
      DEALLOCATE(rushLarsonSolver)
    ENDIF
         
    EXITS("SolverDAE_RushLarsonFinalise")
    RETURN
999 ERRORSEXITS("SolverDAE_RushLarsonFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_RushLarsonFinalise

  !
  !================================================================================================================================
  !

  !>Initialise an Rush-Larson solver for a differential-algebraic equation solver
  SUBROUTINE SolverDAE_RushLarsonInitialise(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to initialise a Rush-Larson solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverDAE_RushLarsonInitialise",err,error,*998)

    IF(ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equation solver is not associated.",err,error,*998)
    IF(ASSOCIATED(daeSolver%rushLarsonSolver)) &
      & CALL FlagError("Rush-Larson solver is already associated for this differential-algebraic equation solver.",err,error,*998)
    
    !Allocate the Rush-Larson solver
    ALLOCATE(daeSolver%rushLarsonSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Rush-Larson solver.",err,error,*999)
    !Initialise
    daeSolver%rushLarsonSolver%DAESolver=>daeSolver
    daeSolver%rushLarsonSolver%solverLibrary=0
    !Defaults
         
    EXITS("SolverDAE_RushLarsonInitialise")
    RETURN
999 CALL SolverDAE_RushLarsonFinalise(daeSolver%rushLarsonSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverDAE_RushLarsonInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_RushLarsonInitialise

  !
  !================================================================================================================================
  !

  !>Solve using a Rush-Larson differential-algebraic equation solver.
  SUBROUTINE SolverDAERushLarson_Solve(rushLarsonSolver,err,error,*)

    !Argument variables
    TYPE(RushLarsonDAESolverType), POINTER :: rushLarsonSolver !<A pointer the Rush-Larson differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverDAERushLarson_Solve",err,error,*999)

    IF(ASSOCIATED(rushLarsonSolver)) &
      & CALL FlagError("Rush-Larson differential-algebraic equation solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
       
    EXITS("SolverDAERushLarson_Solve")
    RETURN
999 ERRORSEXITS("SolverDAERushLarson_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAERushLarson_Solve

  !
  !================================================================================================================================
  !

  !>Solve a differential-algebraic equation solver
  SUBROUTINE SolverDAE_Solve(daeSolver,err,error,*)

    !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer the differential-algebraic equation solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: cellmlIdx,numberOfCellMLEnvironments,solverOutputType
    TYPE(CellMLType), POINTER :: cellML
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(CellMLStateFieldType), POINTER :: cellMLStateField
    TYPE(FieldType), POINTER :: stateField
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDAE_Solve",err,error,*999)

    NULLIFY(solver)
    CALL SolverDAE_SolverGet(daeSolver,solver,err,error,*999)
    CALL Solver_OutputTypeGet(solver,solverOutputType,err,error,*999)
    
    SELECT CASE(daeSolver%daeSolveType)
    CASE(SOLVER_DAE_EULER)
      CALL SolverDAEEuler_Solve(daeSolver%eulerSolver,err,error,*999)
    CASE(SOLVER_DAE_CRANK_NICOLSON)
      CALL SolverDAECrankNicolson_Solve(daeSolver%crankNicolsonSolver,err,error,*999)
    CASE(SOLVER_DAE_RUNGE_KUTTA)
      CALL SolverDAERungeKutta_Solve(daeSolver%rungeKuttaSolver,err,error,*999)
    CASE(SOLVER_DAE_ADAMS_MOULTON)
      CALL SolverDAEAdamsMoulton_Solve(daeSolver%adamsMoultonSolver,err,error,*999)        
    CASE(SOLVER_DAE_BDF)
      CALL SolverDAEBDF_Solve(daeSolver%bdfSolver,err,error,*999)
    CASE(SOLVER_DAE_RUSH_LARSON)
      CALL SolverDAERushLarson_Solve(daeSolver%rushLarsonSolver,err,error,*999)
    CASE(SOLVER_DAE_EXTERNAL)
      CALL SolverDAEExternal_Solve(daeSolver%externalSolver,err,error,*999)
    CASE DEFAULT
      localError="The differential-algebraic equation solver solve type of "// &
        & TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    IF(solverOutputType>SOLVER_SOLVER_OUTPUT) THEN
      
      NULLIFY(cellMLEquations)
      CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
      CALL CellMLEquations_NumberOfCellMLEnvironmentsGet(cellMLEquations,numberOfCellMLEnvironments,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver State vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of CellML environments = ",numberOfCellMLEnvironments,err,error,*999)
      DO cellmlIdx=1,numberOfCellMLEnvironments
        NULLIFY(cellML)
        CALL CellMLEquations_CellMLGet(cellMLEquations,cellMLIdx,cellML,err,error,*999)
        NULLIFY(cellMLStateField)
        CALL CellML_CellMLStateFieldGet(cellML,cellMLStateField,err,error,*999)
        NULLIFY(stateField)
        CALL CellMLStateField_StateFieldGet(cellMLStateField,stateField,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"CellML index : ",cellmlIdx,err,error,*999)
        CALL Field_ParameterSetOutput(GENERAL_OUTPUT_TYPE,stateField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !cellmlIdx
              
    ENDIF
         
    EXITS("SolverDAE_Solve")
    RETURN
999 ERRORSEXITS("SolverDAE_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDAE_Solve

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the solve type for an differential-algebraic equation solver.
  SUBROUTINE Solver_DAESolverTypeSet(solver,daeSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the differential-algebraic equation solver type for 
    INTEGER(INTG), INTENT(IN) :: daeSolveType !<The type of solver for the differential-algebraic equation to set \see SolverRoutines_DAESolverTypes,SolverRoutines.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Solver_DAESolverTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    IF(daeSolveType/=daeSolver%daeSolveType) THEN
      !Intialise the new differential-algebraic equation solver type
      SELECT CASE(daeSolveType)
      CASE(SOLVER_DAE_EULER)
        CALL SolverDAE_EulerInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_CRANK_NICOLSON)
        CALL SolverDAE_CrankNicolsonInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_RUNGE_KUTTA)
        CALL SolverDAE_RungeKuttaInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_ADAMS_MOULTON)
        CALL SolverDAE_AdamsMoultonInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_BDF)
        CALL SolverDAE_BDFInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_RUSH_LARSON)
        CALL SolverDAE_RushLarsonInitialise(daeSolver,err,error,*999)
      CASE(SOLVER_DAE_EXTERNAL)
        CALL SolverDAE_ExternalInitialise(daeSolver,err,error,*999)
      CASE DEFAULT
        localError="The specified differential-algebraic equation solver type of "// &
          & TRIM(NumberToVString(daeSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old differential-algebraic equation solver type
      SELECT CASE(daeSolver%daeSolveType)
      CASE(SOLVER_DAE_EULER)
        CALL SolverDAE_EulerFinalise(daeSolver%eulerSolver,err,error,*999)
      CASE(SOLVER_DAE_CRANK_NICOLSON)
        CALL SolverDAE_CrankNicolsonFinalise(daeSolver%crankNicolsonSolver,err,error,*999)
      CASE(SOLVER_DAE_RUNGE_KUTTA)
        CALL SolverDAE_RungeKuttaFinalise(daeSolver%rungeKuttaSolver,err,error,*999)
      CASE(SOLVER_DAE_ADAMS_MOULTON)
        CALL SolverDAE_AdamsMoultonFinalise(daeSolver%adamsMoultonSolver,err,error,*999)
      CASE(SOLVER_DAE_BDF)
        CALL SolverDAE_BDFFinalise(daeSolver%bdfSolver,err,error,*999)
      CASE(SOLVER_DAE_RUSH_LARSON)
        CALL SolverDAE_RushLarsonFinalise(daeSolver%rushLarsonSolver,err,error,*999)
      CASE(SOLVER_DAE_EXTERNAL)
        CALL SolverDAE_ExternalFinalise(daeSolver%externalSolver,err,error,*999)
      CASE DEFAULT
        localError="The differential-algebraic equation solve type of "// &
          & TRIM(NumberToVString(daeSolver%daeSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      daeSolver%daeSolveType=daeSolveType
    ENDIF
         
    EXITS("Solver_DAESolverTypeSet")
    RETURN
999 ERRORSEXITS("Solver_DAESolverTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAESolverTypeSet

  !
  !================================================================================================================================
  !

  !>Set/change the times for a differential-algebraic equation solver
  SUBROUTINE Solver_DAETimesSet(solver,startTime,endTime,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the differential-algebraic equation solver to set the times for
    REAL(DP), INTENT(IN) :: startTime !<The start time for the differential equation solver
    REAL(DP), INTENT(IN) :: endTime !<The end time for the differential equation solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_DAETimesSet",err,error,*999)

    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    IF(endTime<=startTime) THEN
      localError="The specified end time of "//TRIM(NumberToVString(endTime,"*",err,error))// &
        & " is not > than the specified start time of "//TRIM(NumberToVString(startTime,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    daeSolver%startTime=startTime
    daeSolver%endTime=endTime
         
    EXITS("Solver_DAETimesSet")
    RETURN
999 ERRORSEXITS("Solver_DAETimesSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAETimesSet

  !
  !================================================================================================================================
  !

  !>Set/change the (initial) time step size for a differential-algebraic equation solver
  SUBROUTINE Solver_DAETimeStepSet(solver,timeStep,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the differential-algebraic equation solver to set the times for
    REAL(DP), INTENT(IN) :: timeStep !<The (initial) time step for the differential-algebraic equation solver    
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DAESolverType), POINTER :: daeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_DAETimeStepSet",err,error,*999)

    CALL Solver_AssertIsDAE(solver,err,error,*999)
    NULLIFY(daeSolver)
    CALL Solver_DAESolverGet(solver,daeSolver,err,error,*999)
    IF(ABS(timeStep)<=ZERO_TOLERANCE) THEN
      localError="The specified time step of "//TRIM(NumberToVString(timeStep,"*",err,error))// &
        & " is invalid. The time step must be > zero."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    daeSolver%initialStep=timeStep
         
    EXITS("Solver_DAETimeStepSet")
    RETURN
999 ERRORSEXITS("Solver_DAETimeStepSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DAETimeStepSet

  !
  !================================================================================================================================
  !

  !>Destroys a solver
  SUBROUTINE Solver_Destroy(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
        
    EXITS("Solver_Destroy")
    RETURN
999 ERRORSEXITS("Solver_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_Destroy

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a dynamic solver 
  SUBROUTINE SolverDynamic_CreateFinish(dynamicSolver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dampingMatrixNumber,dependentVariableType,dynamicVariableType,equationsSetIdx,interfaceConditionIdx, &
      & lagrangeVariableType,linearLibraryType,linearMatrixIdx,massMatrixNumber,nonlinearLibraryType,numberOfEquationsSets, &
      & numberOfInterfaceConditions,numberOfLinearMatrices,numberOfResiduals,numberOfResidualVariables,numberOfSources, &
      & residualIdx,residualVariableIdx,rhsVariableType,solverVariableIdx,sourceIdx,sparsityType,symmetryType,variablePositionIdx
    LOGICAL :: lumped
    TYPE(DistributedVectorType), POINTER :: tempDistributedVector
    TYPE(DomainMappingType), POINTER :: colsDomainMapping,rowsDomainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,linearMatrix,massMatrix,stiffnessMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,dynamicVariable,lagrangeVariable,lhsVariable,linearVariable, &
      & residualVariable,rhsVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(SolverType), POINTER :: linearSolver,nonlinearSolver,solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDynamic_CreateFinish",err,error,*999)

    NULLIFY(solver)
    CALL SolverDynamic_SolverGet(dynamicSolver,solver,err,error,*999)
    SELECT CASE(dynamicSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      !Create the parameter sets required for the solver
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(solverMatrixToEquationsMap)
      CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,1,solverMatrixToEquationsMap,err,error,*999)
      NULLIFY(solverVariables)
      CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverVariables,err,error,*999)
      !Initialise for explicit solve
      dynamicSolver%explicit=ABS(dynamicSolver%theta(dynamicSolver%degree))<ZERO_TOLERANCE
      !Loop over the equations set in the solver equations
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(dependentVariable)
        dependentVariableType=0
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        NULLIFY(lhsVariable)
        CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
        NULLIFY(rowsDomainMapping)
        CALL FieldVariable_DomainMappingGet(lhsVariable,rowsDomainMapping,err,error,*999)
        NULLIFY(vectorMatrices)
        CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
        !Dynamic variables and matrices
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(dynamicVariable)
        dynamicVariableType=0
        IF(ASSOCIATED(dynamicMapping)) THEN
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          !Sanity check that the dynamic variable is in the list of solver variables
          NULLIFY(solverVariable)
          CALL SolverMappingVariables_VariableInListCheck(solverVariables,dynamicVariable,solverVariable,err,error,*999)
          IF(.NOT.ASSOCIATED(solverVariable)) THEN
            localError="The dynamic variable for equations set number "// &
              & TRIM(NumberToVString(equationsSet%userNumber,"*",err,error))// &
              & " is not present in the list of solver matrix variables."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          CALL FieldVariable_VariableTypeGet(dynamicVariable,dynamicVariableType,err,error,*999)
          !Set up the parameter sets to hold the required solver parameters
          !1st degree or higher so set up displacement parameter sets
          CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)          
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN                   
            !2nd degree or higher so set up velocity parameter sets
            CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              !3rd degree or higher so set up acceleration parameter sets
              CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(dynamicVariable,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
                & err,error,*999)
            ENDIF
          ENDIF                            
          !Create the dynamic matrices temporary vector for matrix-vector products
          NULLIFY(dynamicMatrices)
          CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
          NULLIFY(tempDistributedVector)
          CALL EquationsMatricesDynamic_TempDistributedVectorExists(dynamicMatrices,tempDistributedVector,err,error,*999)
          IF(.NOT.ASSOCIATED(tempDistributedVector)) THEN
            CALL DistributedVector_CreateStart(rowsDomainMapping,dynamicMatrices%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(dynamicMatrices%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(dynamicMatrices%tempVector,err,error,*999)
          ENDIF
          !Check to see if we have an explicit solve
          IF(dynamicSolver%explicit) THEN
            CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
            IF(dampingMatrixNumber/=0) THEN
              NULLIFY(dampingMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix,err,error,*999)
              CALL EquationsMatrix_LumpedFlagGet(dampingMatrix,lumped,err,error,*999)
              dynamicSolver%explicit=dynamicSolver%explicit.AND.lumped
            ENDIF
            CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
            IF(massMatrixNumber/=0) THEN
              NULLIFY(massMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
              CALL EquationsMatrix_LumpedFlagGet(massMatrix,lumped,err,error,*999)
              dynamicSolver%explicit=dynamicSolver%explicit.AND.lumped
            ENDIF
          ENDIF
        ENDIF
          
        !We now allow for static equation sets for dynamic solvers to be able to couple static eqs - dynamic eqs
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
        IF(ASSOCIATED(nonlinearMapping)) THEN
          IF(dynamicSolver%linearity/=SOLVER_DYNAMIC_NONLINEAR) THEN
            localError="The specified dynamic solver linearity type of "// &
              & TRIM(NumberToVString(dynamicSolver%linearity,"*",err,error))// &
              & " is invalid for a nonlinear equations mapping."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          NULLIFY(nonlinearMatrices)
          CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
          !Loop over the residuals
          CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
          DO residualIdx=1,numberOfResiduals
            NULLIFY(residualMapping)
            CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
            CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
            DO residualVariableIdx=1,numberOfResidualVariables
              NULLIFY(residualVariable)
              CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
              !See if the residual variable is in the list of solver variables
              NULLIFY(solverVariable)
              CALL SolverMappingVariables_VariableInListCheck(solverVariables,residualVariable,solverVariable,err,error,*999)
              IF(ASSOCIATED(solverVariable)) THEN
                CALL FieldVariable_ParameterSetEnsureCreated(residualVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(residualVariable,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE, &
                  & err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(residualVariable,FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(residualVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
              ENDIF
            ENDDO !residualVariableIdx
            NULLIFY(residualVector)
            CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
            CALL EquationsMatricesResidual_DistributedVectorEnsureCreated(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
              & err,error,*999)
            CALL EquationsMatricesResidual_DistributedVectorEnsureCreated(residualVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
              CALL EquationsMatricesResidual_DistributedVectorEnsureCreated(residualVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                & err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                CALL EquationsMatricesResidual_DistributedVectorEnsureCreated(residualVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                  & err,error,*999)
              ENDIF
            ENDIF
          ENDDO !residualIdx
          NULLIFY(tempDistributedVector)
          CALL EquationsMatricesNonlinear_TempDistributedVectorExists(nonlinearMatrices,tempDistributedVector,err,error,*999)
          IF(.NOT.ASSOCIATED(tempDistributedVector)) THEN
            CALL DistributedVector_CreateStart(rowsDomainMapping,nonlinearMatrices%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(nonlinearMatrices%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(nonlinearMatrices%tempVector,err,error,*999)
          ENDIF
        ENDIF
      
        !Check if there are any linear mappings
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
          DO linearMatrixIdx=1,numberOfLinearMatrices
            NULLIFY(linearMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,linearMatrixIdx,linearMatrix,err,error,*999)
            NULLIFY(linearVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixIdx,linearVariable,err,error,*999)
            !See if the linear varible is in the list of solver variables
            NULLIFY(solverVariable)
            CALL SolverMappingVariables_VariableInListCheck(solverVariables,linearVariable,solverVariable,err,error,*999)
            IF(ASSOCIATED(solverVariable)) THEN
              CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE, &
                & err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN                   
                !2nd degree or higher so set up velocity parameter sets
                CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  !3rd degree or higher so set up acceleration parameter sets
                  CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetEnsureCreated(linearVariable,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
                    & err,error,*999)
                ENDIF
              ENDIF
            ENDIF
            !If there are any linear matrices create temporary vector for matrix-vector products
            NULLIFY(tempDistributedVector)
            CALL EquationsMatrix_TempDistributedVectorExists(linearMatrix,tempDistributedVector,err,error,*999)
            IF(.NOT.ASSOCIATED(tempDistributedVector)) THEN
              CALL DistributedVector_CreateStart(rowsDomainMapping,linearMatrix%tempVector,err,error,*999)
              CALL DistributedVector_DataTypeSet(linearMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
              CALL DistributedVector_CreateFinish(linearMatrix%tempVector,err,error,*999)
            ENDIF
          ENDDO !linearMatrixIdx
          NULLIFY(tempDistributedVector)
          CALL EquationsMatricesLinear_TempDistributedVectorExists(linearMatrices,tempDistributedVector,err,error,*999)
          IF(.NOT.ASSOCIATED(tempDistributedVector)) THEN
            CALL DistributedVector_CreateStart(rowsDomainMapping,linearMatrices%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(linearMatrices%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(linearMatrices%tempVector,err,error,*999)
          ENDIF
        ENDIF

        NULLIFY(sourcesMapping)
        CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
        IF(ASSOCIATED(sourcesMapping)) THEN
          CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
          CALL EquationsMappingSources_NumberOfSourcesGet(sourcesMapping,numberOfSources,err,error,*999)
          DO sourceIdx=1,numberOfSources
            NULLIFY(sourceVector)
            CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
            CALL EquationsMatricesSource_DistributedVectorEnsureCreated(sourceVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
              & err,error,*999)
            CALL EquationsMatricesSource_DistributedVectorEnsureCreated(sourceVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
              CALL EquationsMatricesSource_DistributedVectorEnsureCreated(sourceVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                & err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                CALL EquationsMatricesSource_DistributedVectorEnsureCreated(sourceVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                  & err,error,*999)
              ENDIF
            ENDIF
          ENDDO !sourceIdx
          NULLIFY(tempDistributedVector)
          CALL EquationsMatricesSources_TempDistributedVectorExists(sourceVectors,tempDistributedVector,err,error,*999)
          IF(.NOT.ASSOCIATED(tempDistributedVector)) THEN
            CALL DistributedVector_CreateStart(rowsDomainMapping,sourceVectors%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(sourceVectors%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(sourceVectors%tempVector,err,error,*999)
          ENDIF
        ENDIF
        
        NULLIFY(rhsMapping)
        CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
        IF(ASSOCIATED(rhsMapping)) THEN
          NULLIFY(rhsVariable)
          CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
          CALL FieldVariable_VariableTypeGet(rhsVariable,rhsVariableType,err,error,*999)
          NULLIFY(rhsVector)
          CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
          CALL EquationsMatricesRHS_DistributedVectorEnsureCreated(rhsVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
            & err,error,*999)
          CALL EquationsMatricesRHS_DistributedVectorEnsureCreated(rhsVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
            & err,error,*999)
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
            CALL EquationsMatricesRHS_DistributedVectorEnsureCreated(rhsVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL EquationsMatricesRHS_DistributedVectorEnsureCreated(rhsVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                & err,error,*999)
            ENDIF
          ENDIF
        ENDIF
        
      ENDDO !equationsSetIdx
      
      !Loop over any interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
        !See if the linear varible is in the list of solver variables
        NULLIFY(solverVariable)
        CALL SolverMappingVariables_VariableInListCheck(solverVariables,lagrangeVariable,solverVariable,err,error,*999)
        IF(ASSOCIATED(solverVariable)) THEN
          CALL FieldVariable_VariableTypeGet(lagrangeVariable,lagrangeVariableType,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
            !2nd degree or higher so set up velocity parameter sets
            CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              !3rd degree or higher so set up acceleration parameter sets
              CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetEnsureCreated(lagrangeVariable,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
                & err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO !interfaceConditionIdx
                          
      !Create the solver matrices and vectors
      IF(dynamicSolver%linearity==SOLVER_DYNAMIC_LINEAR) THEN
        NULLIFY(linearSolver)
        CALL SolverDynamic_LinkedLinearSolverGet(dynamicSolver,linearSolver,err,error,*999)
        CALL Solver_LibraryTypeGet(linearSolver,linearLibraryType,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
        CALL SolverMatrices_LibraryTypeSet(solverMatrices,linearLibraryType,err,error,*999)
        IF(dynamicSolver%explicit) THEN
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE,err,error,*999)
        ELSE
          CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(SOLVER_SPARSE_MATRICES)
            CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
          CASE(SOLVER_FULL_MATRICES)
            CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
          SELECT CASE(symmetryType)
          CASE(SOLVER_SYMMETRIC_MATRICES)
            CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
          CASE(SOLVER_UNSYMMETRIC_MATRICES)
            CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
        CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
        !Link linear solver
        linearSolver%solverEquations=>solver%solverEquations
        !Finish the creation of the linear solver
        CALL SolverLinear_CreateFinish(linearSolver%linearSolver,err,error,*999)
      ELSE IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
        NULLIFY(nonlinearSolver)
        CALL SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,nonlinearSolver,err,error,*999)
        CALL Solver_LibraryTypeGet(nonlinearSolver,nonlinearLibraryType,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
        CALL SolverMatrices_LibraryTypeSet(solverMatrices,nonlinearLibraryType,err,error,*999)
        IF(dynamicSolver%explicit) THEN
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE,err,error,*999)
        ELSE
          CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(SOLVER_SPARSE_MATRICES)
            CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
          CASE(SOLVER_FULL_MATRICES)
            CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
          SELECT CASE(symmetryType)
          CASE(SOLVER_SYMMETRIC_MATRICES)
            CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
          CASE(SOLVER_UNSYMMETRIC_MATRICES)
            CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
        CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
        !Link nonlinear solver
        nonlinearSolver%solverEquations=>solver%solverEquations
        !Finish the creation of the nonlinear solver
        CALL SolverNonlinear_CreateFinish(nonlinearSolver%nonlinearSolver,err,error,*999)
      ENDIF
    CASE(SOLVER_PETSC_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "//TRIM(NumberToVString(dynamicSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverDynamic_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverDynamic_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDynamic_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the degree of the polynomial used to interpolate time for a dynamic solver.
  SUBROUTINE Solver_DynamicDegreeSet(solver,degree,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    INTEGER(INTG), INTENT(IN) :: degree !<The degree of the polynomial used for time interpolation in a dynamic solver \see SolverRoutines_DynamicDegreeTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: degreeIdx
    REAL(DP), ALLOCATABLE :: newTheta(:)
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_DynamicDegreeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    IF(degree/=dynamicSolver%degree) THEN
      IF(degree>=dynamicSolver%order) THEN
        localError="Invalid dynamic solver setup. The specfied degree of "// &
          & TRIM(NumberToVString(degree,"*",err,error))//" must be >= the current dynamic order of "// &
          & TRIM(NumberToVString(dynamicSolver%order,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      SELECT CASE(degree)
      CASE(SOLVER_DYNAMIC_FIRST_DEGREE,SOLVER_DYNAMIC_SECOND_DEGREE,SOLVER_DYNAMIC_THIRD_DEGREE)
        ALLOCATE(newTheta(degree),STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate new theta.",err,error,*999)
        DO degreeIdx=1,dynamicSolver%degree
          newTheta(degreeIdx)=dynamicSolver%theta(degreeIdx)
        ENDDO !degreeIdx
        IF(degree>dynamicSolver%degree) THEN
          DO degreeIdx=dynamicSolver%degree+1,degree
            newTheta(degreeIdx)=1.0_DP
          ENDDO !degreeIdx
        ENDIF
        CALL MOVE_ALLOC(newTheta,dynamicSolver%theta)
        dynamicSolver%degree=degree
      CASE DEFAULT
        localError="The specified degree of "//TRIM(NumberToVString(degree,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_DynamicDegreeSet")
    RETURN
999 IF(ALLOCATED(newTheta)) DEALLOCATE(newTheta)
    ERRORSEXITS("Solver_DynamicDegreeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicDegreeSet

  !
  !================================================================================================================================
  !

  !>Finalise a dynamic solver and deallocates all memory
  RECURSIVE SUBROUTINE Solver_DynamicFinalise(dynamicSolver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer the dynamic solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_DynamicFinalise",err,error,*999)
    
    IF(ASSOCIATED(dynamicSolver)) THEN
      IF(ALLOCATED(dynamicSolver%theta)) DEALLOCATE(dynamicSolver%theta)
      CALL Solver_Finalise(dynamicSolver%linearSolver,err,error,*999)
      CALL Solver_Finalise(dynamicSolver%nonlinearSolver,err,error,*999)
      DEALLOCATE(dynamicSolver)
    ENDIF
        
    EXITS("Solver_DynamicFinalise")
    RETURN
999 ERRORSEXITS("Solver_DynamicFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DynamicFinalise
 
  !
  !================================================================================================================================
  !

  !>Initialise a dynamic solver for a solver.
  SUBROUTINE Solver_DynamicInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the dynamic solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_DynamicInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%dynamicSolver)) CALL FlagError("Dynamic solver is already associated for this solver.",err,error,*998)
    
    !Allocate memory for dynamic solver and set default values (link solver later on)
    ALLOCATE(solver%dynamicSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver dynamic solver.",err,error,*999)
    solver%dynamicSolver%solver=>solver
    solver%dynamicSolver%solverLibrary=SOLVER_CMISS_LIBRARY
    solver%dynamicSolver%solverInitialised=.FALSE.
    solver%dynamicSolver%numberOfSolves=0
    solver%dynamicSolver%order=SOLVER_DYNAMIC_FIRST_ORDER
    solver%dynamicSolver%degree=SOLVER_DYNAMIC_FIRST_DEGREE
    solver%dynamicSolver%scheme=SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME
    solver%dynamicSolver%startupType=SOLVER_DYNAMIC_PREVIOUS_STARTUP_TYPE
    ALLOCATE(solver%dynamicSolver%theta(1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate theta.",err,error,*999)
    solver%dynamicSolver%theta(1)=1.0_DP/2.0_DP
    solver%dynamicSolver%explicit=.FALSE.
    solver%dynamicSolver%restart=.FALSE.
    solver%dynamicSolver%ale=.TRUE. !this should be .FALSE. eventually and set by the user
    solver%dynamicSolver%currentTime=0.0_DP
    solver%dynamicSolver%timeIncrement=0.01_DP
    NULLIFY(solver%dynamicSolver%linearSolver)
    NULLIFY(solver%dynamicSolver%nonlinearSolver)
    !Make a linear solver by default, and allocate solver%linearSolver
    CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_LINEAR,err,error,*999)
        
    EXITS("Solver_DynamicInitialise")
    RETURN
999 CALL Solver_DynamicFinalise(solver%dynamicSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_DynamicInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_DynamicInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a dynamic solver.
  SUBROUTINE SolverDynamic_LibraryTypeSet(dynamicSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer the dynamic solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the dynamic solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDynamic_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
    
    SELECT CASE(solverLibraryType)
    CASE(SOLVER_CMISS_LIBRARY)
      dynamicSolver%solverLibrary=SOLVER_CMISS_LIBRARY
    CASE DEFAULT
      localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
        & " is invalid for a dynamic solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverDynamic_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverDynamic_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverDynamic_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for the dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicLinearityTypeSet
  SUBROUTINE Solver_DynamicLinearityTypeSet(solver,linearityType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the dynamic solver for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The type of linearity to be set \see SolverRoutines_EquationLinearityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearSolveType,nonlinearSolveType
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(SolverType), POINTER :: linearSolver
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("Solver_DynamicLinearityTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)

    IF(linearityType/=dynamicSolver%linearity) THEN

      CALL Solver_LinkedSolverRemove(solver,SOLVER_LINEAR_TYPE,err,error,*999)
      CALL Solver_Finalise(dynamicSolver%linearSolver,err,error,*999)
      CALL Solver_Finalise(dynamicSolver%nonlinearSolver,err,error,*999)
      
      SELECT CASE(linearityType)
      CASE(SOLVER_DYNAMIC_LINEAR)
        ALLOCATE(dynamicSolver%linearSolver,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver linear solver.",err,error,*999)
        NULLIFY(dynamicSolver%linearSolver%solvers)
        dynamicSolver%linearity=SOLVER_DYNAMIC_LINEAR
        CALL Solver_Initialise(dynamicSolver%linearSolver,err,error,*999)
        CALL Solver_LinearInitialise(dynamicSolver%linearSolver,err,error,*999)
        CALL Solver_LinkedSolverAdd(solver,dynamicSolver%linearSolver,SOLVER_LINEAR_TYPE,err,error,*999)
        CALL SolverLinear_SolverTypeGet(dynamicSolver%linearSolver%linearSolver,linearSolveType,err,error,*999)
        IF(linearSolveType==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) &
          & CALL Solver_LinearIterativeSolutionInitialiseTypeSet(dynamicSolver%linearSolver,SOLVER_SOLUTION_INITIALISE_ZERO, &
          & err,error,*999)
      CASE(SOLVER_DYNAMIC_NONLINEAR)
        ALLOCATE(dynamicSolver%nonlinearSolver,STAT=err)
        IF(err/=0) CALL FlagError("Could not allocate solver nonlinear solver.",err,error,*999)
        NULLIFY(dynamicSolver%nonlinearSolver%solvers)
        dynamicSolver%linearity=SOLVER_DYNAMIC_NONLINEAR
        CALL Solver_Initialise(dynamicSolver%nonlinearSolver,err,error,*999)
        CALL Solver_NonlinearInitialise(dynamicSolver%nonlinearSolver,err,error,*999)
        CALL Solver_LinkedSolverAdd(solver,dynamicSolver%nonlinearSolver,SOLVER_NONLINEAR_TYPE,err,error,*999)
        CALL SolverNonlinear_SolverTypeGet(dynamicSolver%nonlinearSolver%nonlinearSolver,nonlinearSolveType,err,error,*999)
        IF(nonlinearSolveType==SOLVER_NONLINEAR_NEWTON) THEN
          CALL Solver_NewtonSolutionInitialiseTypeSet(dynamicSolver%nonlinearSolver,SOLVER_SOLUTION_INITIALISE_ZERO,err,error,*999)
          NULLIFY(newtonSolver)
          CALL SolverNonlinear_NewtonSolverGet(dynamicSolver%nonlinearSolver%nonlinearSolver,newtonSolver,err,error,*999)
          NULLIFY(linearSolver)
          CALL SolverNonlinearNewton_LinkedLinearSolverGet(newtonSolver,linearSolver,err,error,*999)
          CALL SolverLinear_SolverTypeGet(linearSolver%linearSolver,linearSolveType,err,error,*999)
          IF(linearSolveType==SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE) &
            & CALL Solver_LinearIterativeSolutionInitialiseTypeSet(linearSolver,SOLVER_SOLUTION_INITIALISE_ZERO,err,error,*999)
        ENDIF        
      CASE DEFAULT
        localError="The specified solver equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_DynamicLinearityTypeSet")
    RETURN
999 ERRORSEXITS("Solver_DynamicLinearityTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicLinearityTypeSet

  !
  !================================================================================================================================
  ! 

  !>Copies the current to previous time-step, calculates mean predicted values, predicted values and previous residual values.
  SUBROUTINE Solver_DynamicMeanPredictedCalculate(solver,err,error,*)

    !Argument variableg
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dynamicVariableType,equationsSetIdx,linearMatrixIdx,numberOfEquationsSets,numberOfLinearMatrices &
      & ,numberOfResiduals,numberOfResidualVariables,residualIdx,residualVariableIdx,variablePositionIdx
    REAL(DP) :: deltaT,firstMeanPredictionFactor,secondMeanPredictionFactor,thirdMeanPredictionFactor
    REAL(DP) :: firstPredictionFactor,secondPredictionFactor,thirdPredictionFactor
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dynamicVariable,linearVariable,residualVariable,solverVariable
    TYPE(FieldVariablesListType), POINTER :: processedVariablesList
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Solver_DynamicMeanPredictedCalculate",err,error,*999)

    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    IF(dynamicSolver%solverInitialised) THEN
      deltaT=dynamicSolver%timeIncrement
      SELECT CASE(dynamicSolver%degree)
      CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
        firstMeanPredictionFactor=1.0_DP
        firstPredictionFactor=1.0_DP
      CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
        firstMeanPredictionFactor=1.0_DP
        secondMeanPredictionFactor=dynamicSolver%theta(1)*deltaT
        firstPredictionFactor=1.0_DP
        secondPredictionFactor=deltaT
      CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
        firstMeanPredictionFactor=1.0_DP
        secondMeanPredictionFactor=dynamicSolver%theta(1)*deltaT
        thirdMeanPredictionFactor=dynamicSolver%theta(2)*deltaT*deltaT/2.0_DP
        firstPredictionFactor=1.0_DP
        secondPredictionFactor=deltaT
        thirdPredictionFactor=deltaT*deltaT/2.0_DP
      CASE DEFAULT
        localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    IF(dynamicSolver%solverInitialised.OR.(.NOT.dynamicSolver%solverInitialised.AND. &
      & ((dynamicSolver%order==SOLVER_DYNAMIC_FIRST_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
      & (dynamicSolver%order==SOLVER_DYNAMIC_SECOND_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE)))) &
      & THEN
      !We now need to calculate the the mean predicited displacement, velocity and accelerations for any field variables in
      !the solver matrix. We only want to do this to a field variable once and so we need to be careful that we only process
      !variables that are shared between equations sets or dynamic, nonlinear and linear matrices once.
      NULLIFY(processedVariablesList)
      CALL FieldVariablesList_CreateStart(processedVariablesList,err,error,*999)
      CALL FieldVariablesList_CreateFinish(processedVariablesList,err,error,*999)
      !Loop over the equations sets
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMatrices)
        CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(dynamicVariable)
        IF(ASSOCIATED(dynamicMapping)) THEN
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          CALL FieldVariable_VariableTypeGet(dynamicVariable,dynamicVariableType,err,error,*999)
          IF(dynamicSolver%solverInitialised) THEN
            !Check if we have already processed this variable
            CALL FieldVariablesList_VariableInListCheck(processedVariablesList,dynamicVariable,variablePositionIdx,err,error,*999)
            IF(variablePositionIdx/=0) THEN
              !CPB: THIS MESSES UP DYNANIC BC'S. HOWEVER, IT CLEARLY WAS ADDED FOR A REASON (I THINK TO DO
              !     WITH MONODOMAIN AND SPLITTING. COMMENT FOR NOW.
              !
!!As the dynamic solver may be part of a workflow of solvers within a control loop it is possible
!!that the current dependent field values are not equal to the current previous values that were set
!!at the beginning of the control loop. 
!!Copy the current field values to the previous values
              !CALL Field_ParameterSetsCopy(dependentField,dynamicVariableType,FIELD_VALUES_SET_TYPE, &
              !  & FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,err,error,*999)
              IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) &
                & CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_RESIDUAL_SET_TYPE,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
                & 1.0_DP,err,error,*999)
              !Calculate the mean predicted and predicted values for this dependent field.
              SELECT CASE(dynamicSolver%degree)
              CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                !The mean predicited displacement is the current displacement
                CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                  & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                  !The predicted displacement is just the current displacement
                  CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                    & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                ENDIF
              CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                !The mean predicted displacement comes from the previous displacement and the previous velocity
                CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor], &
                  & [FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                  & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                !The mean predicted velocity is the current velocity
                CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                  & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
                IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                  !The predicted displacement comes from the previous displacement and the previous velocity
                  CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstPredictionFactor,secondPredictionFactor], &
                    & [FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                    & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                ENDIF
              CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                !The mean predicted displacement comes from the previous displacement and the previous
                !velocity and acceleration
                CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor, &
                  & thirdMeanPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                  & FIELD_PREVIOUS_ACCELERATION_SET_TYPE],FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                !The mean predicted velocity comes from the previous velocity and acceleration
                CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor], &
                  & [FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE], &
                  & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
                !The mean predicted acceleration is the current acceleration
                CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE, &
                  & FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
                IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                  !The predicted displacement comes from the previous displacement and the previous
                  !velocity and acceleration
                  CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstPredictionFactor,secondPredictionFactor, &
                    & thirdPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                    & FIELD_PREVIOUS_ACCELERATION_SET_TYPE],FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="The dynamic solver degree of "// &
                  & TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)                        
              END SELECT
              CALL FieldVariablesList_VariableAdd(processedVariablesList,dynamicVariable,err,error,*999)
            ENDIF
          ENDIF
        ENDIF !dynamic mapping
        
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
        IF(ASSOCIATED(nonlinearMapping)) THEN
          IF(dynamicSolver%solverInitialised) THEN
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables, &
                & err,error,*999)
              DO residualVariableIdx=1,numberOfResidualVariables
                NULLIFY(residualVariable)
                CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
                !Check if the residual variable is in the solver variables
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,residualVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) THEN
                  !Check if the variable has already been processed
                  NULLIFY(solverVariable)
                  CALL FieldVariablesList_VariableInListCheck(processedVariablesList,residualVariable,variablePositionIdx, &
                    & err,error,*999)
                  IF(variablePositionIdx==0) THEN
                    !CPB: THIS MESSES UP DYNANIC BC'S. HOWEVER, IT CLEARLY WAS ADDED FOR A REASON (I THINK TO DO
                    !     WITH MONODOMAIN AND SPLITTING. COMMENT FOR NOW.
                    !
                    !!As the dynamic solver may be part of a workflow of solvers within a control loop it is possible
                    !!that the current dependent field values are not equal to the current previous values that were set
                    !!at the beginning of the control loop. 
                    !!Copy the current field values to the previous values
                    !CALL FieldVariable_ParameterSetsCopy(residualVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                    !  & 1.0_DP,err,error,*999)
                    !IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                    !  CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_RESIDUAL_SET_TYPE, &
                    !    & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
                    !ENDIF
                    !Calculate the mean predicted and predicted values for this dependent field.
                    SELECT CASE(dynamicSolver%degree)
                    CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                      !The mean predicited displacement is the current displacement
                      CALL FieldVariable_ParameterSetsCopy(residualVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                        & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                      IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                        !The predicted displacement is just the current displacement
                        CALL FieldVariable_ParameterSetsCopy(residualVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                          & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                      ENDIF
                    CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                      !The mean predicted displacement comes from the previous displacement and the previous velocity
                      CALL FieldVariable_ParameterSetsAdd(residualVariable,[firstMeanPredictionFactor, &
                        & secondMeanPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                        & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                      !The mean predicted velocity is the current velocity
                      CALL FieldVariable_ParameterSetsCopy(residualVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                        & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
                      IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                        !The predicted displacement comes from the previous displacement and the previous velocity
                        CALL FieldVariable_ParameterSetsAdd(residualVariable,[firstPredictionFactor,secondPredictionFactor], &
                          & [FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                          & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                      END IF
                    CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                      !The mean predicted displacement comes from the previous displacement and the previous
                      !velocity and acceleration
                      CALL FieldVariable_ParameterSetsAdd(residualVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor, &
                        & thirdMeanPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE, &
                        & FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE], &
                        & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                      !The mean predicted velocity comes from the previous velocity and acceleration
                      CALL FieldVariable_ParameterSetsAdd(residualVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor], &
                        & [FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE], &
                        & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
                      !The mean predicted acceleration is the current acceleration
                      CALL FieldVariable_ParameterSetsCopy(residualVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE, &
                        & FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
                      IF(dynamicSolver%LINEARITY==SOLVER_DYNAMIC_NONLINEAR) THEN
                        !The predicted displacement comes from the previous displacement and the previous
                        !velocity and acceleration
                        CALL FieldVariable_ParameterSetsAdd(residualVariable,[firstPredictionFactor,secondPredictionFactor, &
                          & thirdPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE, &
                          & FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE], &
                          & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                      END IF
                    CASE DEFAULT
                      localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                        & " is invalid."
                      CALL FlagError(localError,err,error,*999)                        
                    END SELECT
                  ENDIF
                ENDIF
              ENDDO !residualVariableIdx
            ENDDO !residualIdx
          ENDIF!initialised
        ENDIF !nonlinear mapping
        
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          IF(dynamicSolver%solverInitialised) THEN
            NULLIFY(linearMatrices)
            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
            CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
            DO linearMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(linearMatrix)
              CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,linearMatrixIdx,linearMatrix,err,error,*999)
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixIdx,linearVariable,err,error,*999)
              !See if the linear varible is in the list of solver variables
              NULLIFY(solverMappingVariable)
              CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,linearVariable,solverMappingVariable, &
                & err,error,*999)
              IF(ASSOCIATED(solverMappingVariable)) THEN
                !Check if the variable has already been processed
                NULLIFY(solverVariable)
                CALL FieldVariablesList_VariableInListCheck(processedVariablesList,residualVariable,variablePositionIdx, &
                  & err,error,*999)
                IF(variablePositionIdx==0) THEN
                  !CPB: THIS MESSES UP DYNANIC BC'S. HOWEVER, IT CLEARLY WAS ADDED FOR A REASON (I THINK TO DO
                  !     WITH MONODOMAIN AND SPLITTING. COMMENT FOR NOW.
                  !
!!As the dynamic solver may be part of a workflow of solvers within a control loop it is possible
!!that the current dependent field values are not equal to the current previous values that were set
!!at the beginning of the control loop. 
!!Copy the current field values to the previous values
                  !CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                  !  & 1.0_DP,err,error,*999)
                  !IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                  !  CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_RESIDUAL_SET_TYPE, &
                  !    & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
                  !ENDIF
                  !Calculate the mean predicted and predicted values for this dependent field.
                  SELECT CASE(dynamicSolver%degree)
                  CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                    !The mean predicited displacement is the current displacement
                    CALL FieldVariable_ParameterSetsCopy(linearVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                      & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                    IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                      !The predicted displacement is just the current displacement
                      CALL FieldVariable_ParameterSetsCopy(linearVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                        & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,1.0_DP,err,error,*999)
                    ENDIF
                  CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                    !The mean predicted displacement comes from the previous displacement and the previous velocity
                    CALL FieldVariable_ParameterSetsAdd(linearVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor], &
                      & [FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                      & FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                    !The mean predicted velocity is the current velocity
                    CALL FieldVariable_ParameterSetsCopy(linearVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                      & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,1.0_DP,err,error,*999)
                    IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
                      !The predicted displacement comes from the previous displacement and the previous velocity
                      CALL FieldVariable_ParameterSetsAdd(linearVariable,[firstPredictionFactor,secondPredictionFactor], &
                        & [FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE], &
                        & FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                    ENDIF
                  CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                    !The mean predicted displacement comes from the previous displacement and the previous
                    !velocity and acceleration
                    CALL FieldVariable_ParameterSetsAdd(linearVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor, &
                      & thirdMeanPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                      & FIELD_PREVIOUS_ACCELERATION_SET_TYPE],FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                    !The mean predicted velocity comes from the previous velocity and acceleration
                    CALL FieldVariable_ParameterSetsAdd(linearVariable,[firstMeanPredictionFactor,secondMeanPredictionFactor], &
                      & [FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_ACCELERATION_SET_TYPE], &
                      & FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE,err,error,*999)
                    !The mean predicted acceleration is the current acceleration
                    CALL FieldVariable_ParameterSetsCopy(linearVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE, &
                      & FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
                    IF(dynamicSolver%LINEARITY==SOLVER_DYNAMIC_NONLINEAR) THEN
                      !The predicted displacement comes from the previous displacement and the previous
                      !velocity and acceleration
                      CALL FieldVariable_ParameterSetsAdd(linearVariable,[firstPredictionFactor,secondPredictionFactor, &
                        & thirdPredictionFactor],[FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                        & FIELD_PREVIOUS_ACCELERATION_SET_TYPE],FIELD_PREDICTED_DISPLACEMENT_SET_TYPE,err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)                        
                  END SELECT
                ENDIF
              ENDIF
            ENDDO !linearMatrixIdx
          ENDIF
        ENDIF !linear mapping
      ENDDO !equationsSetIdx
    ENDIF

    EXITS("Solver_DynamicMeanPredictedCalculate")
    RETURN
999 ERRORSEXITS("Solver_DynamicMeanPredictedCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicMeanPredictedCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the restart value  for a dynamic solver.
  SUBROUTINE Solver_DynamicRestartSet(solver,restart,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    LOGICAL, INTENT(IN) :: restart !<The restart value to be set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    
    ENTERS("Solver_DynamicRestartSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    dynamicSolver%restart=restart
   
    EXITS("Solver_DynamicRestartSet")
    RETURN
999 ERRORSEXITS("Solver_DynamicRestartSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicRestartSet

  !
  !================================================================================================================================
  !

  !>Monitors the differential-algebraic equations solve.
  SUBROUTINE SolverDAE_TimeSteppingMonitor(daeSolver,steps,time,err,error,*)

   !Argument variables
    TYPE(DAESolverType), POINTER :: daeSolver !<A pointer to the differential-algebraic equations solver to monitor
    INTEGER(INTG), INTENT(IN) :: steps !<The number of iterations
    REAL(DP), INTENT(IN) :: time !<The current time
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverDAE_TimeSteppingMonitor",err,error,*999)

    IF(.NOT.ASSOCIATED(daeSolver)) CALL FlagError("Differential-algebraic equations solver is not associated.",err,error,*999)
        
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"Differential-algebraic equations solve monitor: ",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
    CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of steps = ",steps,err,error,*999)
    CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Current time    = ",time,err,error,*999)
    
    EXITS("SolverDAE_TimeSteppingMonitor")
    RETURN
999 ERRORSEXITS("SolverDAE_TimeSteppingMonitor",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDAE_TimeSteppingMonitor

  !
  !================================================================================================================================
  !

  !>Sets/changes the order for a dynamic solver.
  SUBROUTINE Solver_DynamicOrderSet(solver,order,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    INTEGER(INTG), INTENT(IN) :: order !<The order of the dynamic solver \see SolverRoutines_DynamicOrderTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_DynamicOrderSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    IF(order==SOLVER_DYNAMIC_SECOND_ORDER.AND.dynamicSolver%degree==SOLVER_DYNAMIC_FIRST_DEGREE) THEN
      localError="Invalid dynamic solver degree. You must have at least a second degree polynomial "// &
        & "interpolation for a second order dynamic solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(order)
    CASE(SOLVER_DYNAMIC_FIRST_ORDER)
      dynamicSolver%order=SOLVER_DYNAMIC_FIRST_ORDER
    CASE(SOLVER_DYNAMIC_SECOND_ORDER)
      dynamicSolver%order=SOLVER_DYNAMIC_SECOND_ORDER
    CASE DEFAULT
      localError="The specified order of "//TRIM(NumberToVString(order,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_DynamicOrderSet")
    RETURN
999 ERRORSEXITS("Solver_DynamicOrderSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicOrderSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the scheme for a dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicSchemeSet
  SUBROUTINE Solver_DynamicSchemeSet(solver,scheme,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the scheme for
    INTEGER(INTG), INTENT(IN) :: scheme !<The scheme used for a dynamic solver \see SolverRoutines_DynamicSchemeTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: alpha,beta,gamma,theta
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_DynamicSchemeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    
    SELECT CASE(scheme)
    CASE(SOLVER_DYNAMIC_EULER_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_EULER_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,0.0_DP,err,error,*999)
    CASE(SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,1.0_DP,err,error,*999)
    CASE(SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,1.0_DP/2.0_DP,err,error,*999)
    CASE(SOLVER_DYNAMIC_GALERKIN_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_GALERKIN_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,2.0_DP/3.0_DP,err,error,*999)
    CASE(SOLVER_DYNAMIC_ZLAMAL_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_ZLAMAL_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[5.0_DP/6.0_DP,2.0_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_SECOND_DEGREE_GEAR_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[3.0_DP/2.0_DP,2.0_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER1_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[1.0848_DP,1.0_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_SECOND_DEGREE_LINIGER2_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[1.2184_DP,1.292_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_NEWMARK1_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_NEWMARK1_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      beta=0.5_DP
      gamma=2.0_DP
      CALL Solver_DynamicThetaSet(solver,[gamma,2.0_DP*beta],err,error,*999)
    CASE(SOLVER_DYNAMIC_NEWMARK2_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_NEWMARK2_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      beta=0.3025_DP
      gamma=0.6_DP
      CALL Solver_DynamicThetaSet(solver,[gamma,2.0_DP*beta],err,error,*999)
    CASE(SOLVER_DYNAMIC_NEWMARK3_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_NEWMARK3_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_SECOND_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      beta=0.25_DP
      gamma=0.5_DP
      CALL Solver_DynamicThetaSet(solver,[gamma,2.0_DP*beta],err,error,*999)
    CASE(SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_THIRD_DEGREE_GEAR_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[2.0_DP,11.0_DP/3.0_DP,6.0_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER1_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[1.84_DP,3.07_DP,4.5_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_THIRD_DEGREE_LINIGER2_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[0.80_DP,1.03_DP,1.29_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_HOUBOLT_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_HOUBOLT_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      CALL Solver_DynamicThetaSet(solver,[2.0_DP,11.0_DP/3.0_DP,6.0_DP],err,error,*999)
    CASE(SOLVER_DYNAMIC_WILSON_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_WILSON_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      theta=1.4_DP
      CALL Solver_DynamicThetaSet(solver,[theta,theta**2,theta**3],err,error,*999)
    CASE(SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_BOSSAK_NEWMARK1_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      alpha=-0.1_DP
      beta=0.3025_DP
      gamma=0.5_DP-alpha
      CALL Solver_DynamicThetaSet(solver,[1.0_DP-alpha,2.0_DP/3.0_DP-alpha+2.0_DP*beta,6.0_DP*beta],err,error,*999)
    CASE(SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_BOSSAK_NEWMARK2_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      alpha=-0.1_DP
      beta=1.0_DP/6.0_DP-1.0_DP/2.0_DP*alpha
      gamma=1.0_DP/2.0_DP-alpha
      CALL Solver_DynamicThetaSet(solver,[1.0_DP-alpha,1.0_DP-2.0_DP*alpha,1.0_DP-3.0_DP*alpha],err,error,*999)
    CASE(SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR1_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      alpha=-0.1_DP
      beta=0.3025_DP
      gamma=0.5_DP-alpha
      CALL Solver_DynamicThetaSet(solver,[1.0_DP,2.0_DP/3.0_DP+2.0_DP*beta-2.0_DP*alpha**2, &
        & 6.0_DP*beta*(1.0_DP+alpha)],err,error,*999)
    CASE(SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_HILBERT_HUGHES_TAYLOR2_SCHEME
      CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_THIRD_DEGREE,err,error,*999)
      CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_SECOND_ORDER,err,error,*999)
      alpha=-0.3_DP
      beta=0.3025_DP
      gamma=0.5_DP-alpha
      CALL Solver_DynamicThetaSet(solver,[1.0_DP,2.0_DP/3.0_DP+2.0_DP*beta-2.0_DP*alpha**2, &
        & 6.0_DP*beta*(1.0_DP+alpha)],err,error,*999)
    CASE(SOLVER_DYNAMIC_USER_DEFINED_SCHEME)
      dynamicSolver%scheme=SOLVER_DYNAMIC_USER_DEFINED_SCHEME
    CASE DEFAULT
      localError="The specified scheme of "//TRIM(NumberToVString(scheme,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
   
    EXITS("Solver_DynamicSchemeSet")
    RETURN
999 ERRORSEXITS("Solver_DynamicSchemeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicSchemeSet

  !
  !================================================================================================================================
  !
  
  !>Solve a dynamic solver 
  SUBROUTINE SolverDynamic_Solve(dynamicSolver,err,error,*)

    !Argument variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver !<A pointer to the dynamic solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfSolverMatrices,solverMatrixIdx
    TYPE(DistributedVectorType), POINTER :: distributedSolverVector
    TYPE(SolverType), POINTER :: linearSolver,solver,nonlinearSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverDynamic_Solve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(dynamicSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
    
    SELECT CASE(dynamicSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      NULLIFY(solver)
      CALL SolverDynamic_SolverGet(dynamicSolver,solver,err,error,*999)
      SELECT CASE(dynamicSolver%linearity)
      CASE(SOLVER_DYNAMIC_LINEAR)
        !Solve the linear dynamic problem
        NULLIFY(linearSolver)
        CALL SolverDynamic_LinkedLinearSolverGet(dynamicSolver,linearSolver,err,error,*999)
        IF(dynamicSolver%solverInitialised) THEN
          !Assemble the solver equations
          CALL Solver_DynamicMeanPredictedCalculate(solver,err,error,*999)
          CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_LINEAR_ONLY,err,error,*999)
          !Solve the linear system
          CALL Solver_Solve(linearSolver,err,error,*999)
          !Update dependent field with solution
          CALL Solver_VariablesDynamicFieldUpdate(solver,err,error,*999)
        ELSE
          !If we need to initialise the solver
          IF((dynamicSolver%order==SOLVER_DYNAMIC_FIRST_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
            & (dynamicSolver%order==SOLVER_DYNAMIC_SECOND_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE)) THEN
            !Assemble the solver equations
            CALL Solver_DynamicMeanPredictedCalculate(solver,err,error,*999)
            CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_LINEAR_ONLY,err,error,*999)
            !Solve the linear system
            CALL Solver_Solve(linearSolver,err,error,*999)
            !Update dependent field with solution
            CALL Solver_VariablesDynamicFieldUpdate(solver,err,error,*999)
          ENDIF
          !Set initialised flag
          dynamicSolver%solverInitialised=.TRUE.
        ENDIF
      CASE(SOLVER_DYNAMIC_NONLINEAR) 
        !Solve the nonlinear dynamic problem
        NULLIFY(nonlinearSolver)
        CALL SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,nonlinearSolver,err,error,*999)
        IF(dynamicSolver%solverInitialised) THEN
          !Calculate predicted values
          CALL Solver_DynamicMeanPredictedCalculate(solver,err,error,*999)
          !Solve the nonlinear system
          CALL Solver_Solve(nonlinearSolver,err,error,*999)
          !Update dependent field with solution
          CALL Solver_VariablesDynamicFieldUpdate(solver,err,error,*999)
        ELSE
          !If we need to initialise the solver
          IF((dynamicSolver%order==SOLVER_DYNAMIC_FIRST_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
            & (dynamicSolver%order==SOLVER_DYNAMIC_SECOND_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE)) THEN
            !No nonlinear solver for the first (starting) time step.
            !Use the nonlinear solvers linear solver to find starting velocities and accelerations etc.
            NULLIFY(linearSolver)
            CALL SolverNonlinear_LinearSolverGet(nonlinearSolver,linearSolver,err,error,*999)
            !Assemble the solver equations
            CALL Solver_DynamicMeanPredictedCalculate(solver,err,error,*999)
            CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY,err,error,*999)
            !Solve
            CALL Solver_Solve(linearSolver,err,error,*999)
          ELSE
            !Assemble the solver equations for the intial values of the residual and RHS.
            CALL Solver_DynamicMeanPredictedCalculate(solver,err,error,*999)
            CALL Solver_DynamicAssemble(solver,SOLVER_MATRICES_RHS_RESIDUAL_ONLY,err,error,*999)
          ENDIF
          !Update dependent field with solution
          CALL Solver_VariablesDynamicFieldUpdate(solver,err,error,*999)
          !Set initialised flag
          dynamicSolver%solverInitialised=.TRUE.
        ENDIF
      CASE DEFAULT
        localError="The dynamic solver linearity type of "// &
          & TRIM(NumberToVString(dynamicSolver%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT

      IF(dynamicSolver%solverInitialised) dynamicSolver%numberOfSolves=dynamicSolver%numberOfSolves+1
          
      IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver solution vectors:",err,error,*999)
        CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of solution vectors = ",numberOfSolverMatrices,err,error,*999)
        DO solverMatrixIdx=1,numberOfSolverMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solution vector for solver matrix : ",solverMatrixIdx,err,error,*999)
          NULLIFY(solverMatrix)
          CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
          NULLIFY(distributedSolverVector)
          CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,distributedSolverVector,err,error,*999)
          CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,distributedSolverVector,err,error,*999)
        ENDDO !solverMatrixIdx
      ENDIF
            
    CASE(SOLVER_PETSC_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "//TRIM(NumberToVString(dynamicSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverDynamic_Solve")
    RETURN
999 ERRORSEXITS("SolverDynamic_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverDynamic_Solve
        
  !
  !================================================================================================================================
  !

  !>Sets/changes a single theta value for a dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicThetaSet
  SUBROUTINE Solver_DynamicThetaSetDP0(solver,theta,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    REAL(DP), INTENT(IN) :: theta !<The theta value to set for the first degree polynomial
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   
    ENTERS("Solver_DynamicThetaSetDP0",err,error,*999)

    CALL Solver_DynamicThetaSetDP1(solver,[theta],err,error,*999)
    
    EXITS("Solver_DynamicThetaSetDP0")
    RETURN
999 ERRORSEXITS("Solver_DynamicThetaSetDP0",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicThetaSetDP0

  !
  !================================================================================================================================
  !

  !>Sets/changes the theta value for a dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicThetaSet
  SUBROUTINE Solver_DynamicThetaSetDP1(solver,theta,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    REAL(DP), INTENT(IN) :: theta(:) !<theta(degreeIdx). The theta value to set for the degreeIdx-1'th polynomial
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: degreeIdx
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_DynamicThetaSetDP1",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)    
    IF(SIZE(theta,1)<dynamicSolver%degree) THEN
      localError="Invalid number of the thetas. The supplied number of thetas ("// &
        & TRIM(NumberToVString(SIZE(theta,1),"*",err,error))//") must be equal to the interpolation degree ("// &
        & TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))//")."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    DO degreeIdx=1,dynamicSolver%degree
      IF(theta(degreeIdx)<0.0_DP) THEN
        localError="The specified theta "//TRIM(NumberToVString(degreeIdx,"*",err,error))// &
          & " value of "//TRIM(NumberToVString(theta(degreeIdx),"*",err,error))// &
          & " is invalid. The theta value must be >= 0.0."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      dynamicSolver%theta(degreeIdx)=theta(degreeIdx)
    ENDDO !degreeIdx
    
    EXITS("Solver_DynamicThetaSetDP1")
    RETURN
999 ERRORSEXITS("Solver_DynamicThetaSetDP1",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicThetaSetDP1

  !
  !================================================================================================================================
  !

  !>Sets/changes the ALE flag for a dynamic solver.
  SUBROUTINE Solver_DynamicALESet(solver,ale,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the theta value for
    LOGICAL :: ale !<The ALE flag for a dynamic solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    
    ENTERS("Solver_DynamicALESet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)    
    dynamicSolver%ale=ale
    
    EXITS("Solver_DynamicALESet")
    RETURN
999 ERRORSEXITS("Solver_DynamicALESet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicALESet

  !
  !================================================================================================================================
  !

  !>Sets/changes the dynamic times for a dynamic solver. \see OpenCMISS::Iron::cmfe_Solver_DynamicTimesSet
  SUBROUTINE Solver_DynamicTimesSet(solver,currentTime,timeIncrement,err,error,*)

   !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the dynamic solver to set the times for
    REAL(DP), INTENT(IN) :: currentTime !<The current time to set
    REAL(DP), INTENT(IN) :: timeIncrement !<The time increment to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_DynamicTimesSet",err,error,*999)

    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)    
    IF(ABS(timeIncrement)<=ZERO_TOLERANCE) THEN
      localError="The specified time increment of "//TRIM(NumberToVString(timeIncrement,"*",err,error))// &
        & " is invalid. The time increment must not be zero."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    dynamicSolver%currentTime=currentTime
    dynamicSolver%timeIncrement=timeIncrement
     
    EXITS("Solver_DynamicTimesSet")
    RETURN
999 ERRORSEXITS("Solver_DynamicTimesSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicTimesSet

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a eigenproblem solver 
  SUBROUTINE SolverEigenproblem_CreateFinish(eigenproblemSolver,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer to the eigenproblem solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEigenproblem_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(eigenproblemSolver)) CALL FlagError("Eigenproblem solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
          
    EXITS("SolverEigenproblem_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverEigenproblem_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEigenproblem_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a eigenproblem solver for a solver.
  SUBROUTINE Solver_EigenproblemFinalise(eigenproblemSolver,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer the eigenproblem solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_EigenproblemFinalise",err,error,*999)

    IF(ASSOCIATED(eigenproblemSolver)) THEN        
      DEALLOCATE(eigenproblemSolver)
    ENDIF
         
    EXITS("Solver_EigenproblemFinalise")
    RETURN
999 ERRORSEXITS("Solver_EigenproblemFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_EigenproblemFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a eigenproblem solver for a solver.
  SUBROUTINE Solver_EigenproblemInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the eigenproblem solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_EigenproblemInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%eigenproblemSolver)) &
      & CALL FlagError("Eigenproblem solver is already associated for this solver.",err,error,*998)
     
    ALLOCATE(solver%eigenproblemSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver eigenproblem solver.",err,error,*999)
    solver%eigenproblemSolver%solver=>solver
    solver%eigenproblemSolver%solverLibrary=SOLVER_PETSC_LIBRARY
         
    EXITS("Solver_EigenproblemInitialise")
    RETURN
999 CALL Solver_EigenproblemFinalise(solver%eigenproblemSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_EigenproblemInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_EigenproblemInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for an eigenproblem solver.
  SUBROUTINE SolverEigenproblem_LibraryTypeSet(eigenproblemSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer the eigenproblem solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the eigenproblem solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverEigenproblem_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(eigenproblemSolver)) CALL FlagError("Dynamic solver is not associated.",err,error,*999)
    
    SELECT CASE(solverLibraryType)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The specified solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
        & " is invalid for an eigenproblem solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("SolverEigenproblem_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverEigenproblem_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEigenproblem_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Solve a eigenproblem solver
  SUBROUTINE SolverEigenproblem_Solve(eigenproblemSolver,err,error,*)

    !Argument variables
    TYPE(EigenproblemSolverType), POINTER :: eigenproblemSolver !<A pointer the eigenproblem solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEigenproblem_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(eigenproblemSolver)) CALL FlagError("Eigenproblem solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverEigenproblem_Solve")
    RETURN
999 ERRORSEXITS("SolverEigenproblem_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEigenproblem_Solve

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating solver equations
  SUBROUTINE SolverEquations_CreateFinish(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver

    ENTERS("SolverEquations_CreateFinish",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not finish solver equations creation for a solver that has been linked.",err,error,*999)
      
    solverEquations%solverEquationsFinished=.TRUE.
         
    EXITS("SolverEquations_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverEquations_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating solver equations
  SUBROUTINE SolverEquations_CreateStart(solver,solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to start the creation of solver equations on
    TYPE(SolverEquationsType), POINTER :: solverEquations !<On return, A pointer the solver equations. Must not be associated on entry
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverEquations_CreateStart",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not start solver equations creation for a solver that has been linked.",err,error,*999)       
    IF(ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is already associated.",err,error,*999)
    
    CALL SolverEquations_Initialise(solver,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverMapping_CreateStart(solver%solverEquations,solverMapping,err,error,*999)
    SELECT CASE(solver%solveType)
    CASE(SOLVER_LINEAR_TYPE)
      CALL SolverMapping_NumberOfSolverMatricesSet(solverMapping,1,err,error,*999)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL SolverMapping_NumberOfSolverMatricesSet(solverMapping,1,err,error,*999)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL SolverMapping_NumberOfSolverMatricesSet(solverMapping,1,err,error,*999)
    CASE(SOLVER_DAE_TYPE)
      CALL SolverMapping_NumberOfSolverMatricesSet(solverMapping,0,err,error,*999)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL SolverMapping_NumberOfSolverMatricesSet(solverMapping,2,err,error,*999)
    CASE DEFAULT
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    solverEquations=>solver%solverEquations
         
    EXITS("SolverEquations_CreateStart")
    RETURN
999 ERRORSEXITS("SolverEquations_CreateStart",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_CreateStart
        
  !
  !================================================================================================================================
  !

  !>Destroys the solver equations
  SUBROUTINE SolverEquations_Destroy(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to destroy.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(solverEquations)) CALL FlagError("Solver equations is not associated.",err,error,*999)
    
    CALL SolverEquations_Finalise(solverEquations,err,error,*999)
        
    EXITS("SolverEquations_Destroy")
    RETURN
999 ERRORSEXITS("SolverEquations_Destroy",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_Destroy
        
  !
  !================================================================================================================================
  !

  !>Adds equations sets to solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_EquationsSetAdd
  SUBROUTINE SolverEquations_EquationsSetAdd(solverEquations,equationsSet,equationsSetIndex,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to add the equations set to.
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to add
    INTEGER(INTG), INTENT(OUT) :: equationsSetIndex !<On exit, the index of the equations set that has been added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsLinearity,equationsTimeDependence,solverLinearity,solverTimeDependence
    LOGICAL :: timeCompatible,linearityCompatible
    TYPE(EquationsType), POINTER :: equations
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverEquations_EquationsSetAdd",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not add an equations set for a solver that has been linked.",err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    timeCompatible=.TRUE.
    linearityCompatible=.TRUE.
    !Check solver equations and equations set time dependence is compatible
    CALL SolverEquations_TimeDependenceTypeGet(solverEquations,solverTimeDependence,err,error,*999)
    CALL SolverEquations_LinearityTypeGet(solverEquations,solverLinearity,err,error,*999)
    CALL Equations_TimeDependenceTypeGet(equations,equationsTimeDependence,err,error,*999)
    CALL Equations_LinearityTypeGet(equations,equationsLinearity,err,error,*999)
    SELECT CASE(solverTimeDependence)
    CASE(SOLVER_EQUATIONS_STATIC,SOLVER_EQUATIONS_QUASISTATIC)
      SELECT CASE(equationsTimeDependence)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)
        !OK
      CASE DEFAULT
        timeCompatible=.FALSE.
      END SELECT
    CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC)
      SELECT CASE(equationsTimeDependence)                  
      CASE(EQUATIONS_STATIC)
        !OK for now, just to test!!!
      CASE(EQUATIONS_QUASISTATIC)
        !Not yet implemented, this needs to be checked to see that it works
        timeCompatible=.FALSE.
        CALL FlagError("Static equations set equations with dynamic solver equations is not yet implemented.",err,error,*999)
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC)
        !OK
      CASE DEFAULT
        timeCompatible=.FALSE.
      END SELECT
    CASE(SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equationsTimeDependence)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC,EQUATIONS_FIRST_ORDER_DYNAMIC)
        !Not implemented, this needs to be checked to see that it works
        !timeCompatible=.FALSE.
        !localError="Static or first order dynamic equations set equations with a second order dynamic "// &
        !  & "solver equations is not yet implemented."
        !CALL FlagError(localError,err,error,*999)
      CASE(EQUATIONS_SECOND_ORDER_DYNAMIC)
        !OK
      CASE DEFAULT
        timeCompatible=.FALSE.
      END SELECT
    CASE DEFAULT
      timeCompatible=.FALSE.
      localError="Invalid time dependence for solver equations, "//TRIM(NumberToVString(solverTimeDependence,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(.NOT.timeCompatible) THEN
      localError="Invalid equations set up. The time dependence of the equations set to add of "// &
        & TRIM(NumberToVString(equationsTimeDependence,"*",err,error))// &
        & " is not compatible with the solver equations time dependence of "// &
        & TRIM(NumberToVString(solverTimeDependence,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !Check solver equations and equations set linearity is compatible
    SELECT CASE(solverLinearity)
    CASE(SOLVER_EQUATIONS_LINEAR)
      SELECT CASE(equationsLinearity)
      CASE(EQUATIONS_LINEAR)
        !OK
      CASE DEFAULT
        linearityCompatible=.FALSE.
      END SELECT
    CASE(SOLVER_EQUATIONS_NONLINEAR)
      SELECT CASE(equationsLinearity)
      CASE(EQUATIONS_LINEAR,EQUATIONS_NONLINEAR)
        !OK
      CASE DEFAULT
        linearityCompatible=.FALSE.
      END SELECT
    CASE DEFAULT
      linearityCompatible=.FALSE.
      localError="The solver equations linearity of "//TRIM(NumberToVString(solverLinearity,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF (.NOT.linearityCompatible) THEN
      localError="Invalid equations set up. The linearity of the equations set to add of "// &
        & TRIM(NumberToVString(equationsLinearity,"*",err,error))// &
        & " is not compatible with the solver equations linearity of "// &
        & TRIM(NumberToVString(solverLinearity,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(timeCompatible.AND.linearityCompatible) &
      & CALL SolverMapping_EquationsSetAdd(solverMapping,equationsSet,equationsSetIndex,err,error,*999)
        
    EXITS("SolverEquations_EquationsSetAdd")
    RETURN
999 ERRORSEXITS("SolverEquations_EquationsSetAdd",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_EquationsSetAdd
        
  !
  !================================================================================================================================
  !

  !>Finalises the solver equations and deallocates all memory.
  SUBROUTINE SolverEquations_Finalise(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_Finalise",err,error,*999)

    IF(ASSOCIATED(solverEquations)) THEN
      IF(ASSOCIATED(solverEquations%solverMapping)) CALL SolverMapping_Destroy(solverEquations%solverMapping,err,error,*999)
      IF(ASSOCIATED(solverEquations%solverMatrices)) CALL SolverMatrices_Destroy(solverEquations%solverMatrices,err,error,*999)
      IF(ASSOCIATED(solverEquations%boundaryConditions)) CALL BoundaryConditions_Destroy(solverEquations%boundaryConditions, &
        & err,error,*999)
    ENDIF

    EXITS("SolverEquations_Finalise")
    RETURN
999 ERRORSEXITS("SolverEquations_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_Finalise
        
  !
  !================================================================================================================================
  !

  !>Initialises the solver equations for a solver.
  SUBROUTINE SolverEquations_Initialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverEquations_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%solverEquations)) CALL FlagError("Solver equations is already associated for this solver.",err,error,*998)
     
    ALLOCATE(solver%solverEquations,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver equations.",err,error,*999)
    solver%solverEquations%solver=>solver
    solver%solverEquations%solverEquationsFinished=.FALSE.
    solver%solverEquations%sparsityType=SOLVER_SPARSE_MATRICES
    solver%solverEquations%symmetryType=SOLVER_UNSYMMETRIC_MATRICES
    NULLIFY(solver%solverEquations%solverMapping)
    NULLIFY(solver%solverEquations%solverMatrices)
    NULLIFY(solver%solverEquations%boundaryConditions)
        
    EXITS("SolverEquations_Initialise")
    RETURN
999 CALL SolverEquations_Finalise(solver%solverEquations,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverEquations_Initialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_Initialise
        
  !
  !================================================================================================================================
  !

  !>Adds an interface condition to the solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_InterfaceConditionAdd
  SUBROUTINE SolverEquations_InterfaceConditionAdd(solverEquations,interfaceCondition,interfaceConditionIndex,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to add the interface condition to.
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition to add
    INTEGER(INTG), INTENT(OUT) :: interfaceConditionIndex !<On exit, the index of the interface condition that has been added
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverMappingType), POINTER :: solverMapping
   
    ENTERS("SolverEquations_InterfaceConditionAdd",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not add an equations set for a solver that has been linked.",err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    CALL SolverMapping_InterfaceConditionAdd(solverMapping,interfaceCondition,interfaceConditionIndex,err,error,*999)
        
    EXITS("SolverEquations_InterfaceConditionAdd")
    RETURN
999 ERRORSEXITS("SolverEquations_InterfaceConditionAdd",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_InterfaceConditionAdd
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the linearity type for solver equations
  SUBROUTINE SolverEquations_LinearityTypeSet(solverEquations,linearityType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to set the linearity type for
    INTEGER(INTG), INTENT(IN) :: linearityType !<The type of linearity to be set \see SolverRoutines_EquationLinearityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverEquations_LinearityTypeSet",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not set equations linearity for a solver that has been linked.",err,error,*999)
      
    SELECT CASE(linearityType)
    CASE(SOLVER_EQUATIONS_LINEAR)
      solverEquations%linearity=SOLVER_EQUATIONS_LINEAR
    CASE(SOLVER_EQUATIONS_NONLINEAR)
      solverEquations%linearity=SOLVER_EQUATIONS_NONLINEAR
    CASE DEFAULT
      localError="The specified solver equations linearity type of "//TRIM(NumberToVString(linearityType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverEquations_LinearityTypeSet")
    RETURN
999 ERRORSEXITS("SolverEquations_LinearityTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_LinearityTypeSet
        
  !
  !================================================================================================================================
  !

  !>Finishes the creation of boundary conditions for the given solver equations
  SUBROUTINE SolverEquations_BoundaryConditionsCreateFinish(solverEquations,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER, INTENT(IN) :: solverEquations !<A pointer to the solver equations to create boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverEquations_BoundaryConditionsCreateFinish",err,error,*999)

    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not finish solver equations creation for a solver that has been linked.",err,error,*999)
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    CALL BoundaryConditions_CreateFinish(boundaryConditions,err,error,*999)        
    !Finish of the solver mapping
    CALL SolverMapping_CreateFinish(solverEquations%solverMapping,err,error,*999)
    !Now finish off with the solver specific actions
    SELECT CASE(solver%solveType)
    CASE(SOLVER_LINEAR_TYPE)
      CALL SolverLinear_CreateFinish(solver%linearSolver,err,error,*999)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL SolverNonlinear_CreateFinish(solver%nonlinearSolver,err,error,*999)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL SolverDynamic_CreateFinish(solver%dynamicSolver,err,error,*999)
    CASE(SOLVER_DAE_TYPE)
      CALL SolverDAE_CreateFinish(solver%DAESolver,err,error,*999)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL SolverEigenproblem_CreateFinish(solver%eigenproblemSolver,err,error,*999)
    CASE DEFAULT
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("SolverEquations_BoundaryConditionsCreateFinish")
    RETURN
999 ERRORS("SolverEquations_BoundaryConditionsCreateFinish",err,error)
    EXITS("SolverEquations_BoundaryConditionsCreateFinish")
    RETURN 1

  END SUBROUTINE SolverEquations_BoundaryConditionsCreateFinish

  !
  !================================================================================================================================
  !

  !>Starts the creation of boundary conditions for the given solver equations, and returns a pointer to the boundary conditions
  SUBROUTINE SolverEquations_BoundaryConditionsCreateStart(solverEquations,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER, INTENT(IN) :: solverEquations !<A pointer to the solver equations to create boundary conditions for
    TYPE(BoundaryConditionsType), POINTER, INTENT(OUT) :: boundaryConditions !<On return, a pointer the boundary conditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverEquations_BoundaryConditionsCreateStart",err,error,*999)

    CALL SolverEquations_AssertIsFinished(solverEquations,err,error,*999)
    IF(ASSOCIATED(solverEquations%boundaryConditions)) &
      & CALL FlagError("Solver equations boundary conditions is already associated.",err,error,*999)
    
    CALL BoundaryConditions_CreateStart(solverEquations,boundaryConditions,err,error,*999)
    
    EXITS("SolverEquations_BoundaryConditionsCreateStart")
    RETURN
999 ERRORS("SolverEquations_BoundaryConditionsCreateStart",err,error)
    EXITS("SolverEquations_BoundaryConditionsCreateStart")
    RETURN 1

  END SUBROUTINE SolverEquations_BoundaryConditionsCreateStart

  !
  !================================================================================================================================
  !

  !>Sets/changes the sparsity type for solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_SparsityTypeSet
  SUBROUTINE SolverEquations_SparsityTypeSet(solverEquations,sparsityType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: sparsityType !<The type of solver equations sparsity to be set \see SolverRoutines_EquationsSparsityTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverEquations_SparsityTypeSet",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not set equations sparsity for a solver that has been linked.",err,error,*999)
     
!!TODO: Maybe set the sparsity in the different types of solvers. e.g., a sparse integrator doesn't mean much.
    SELECT CASE(sparsityType)
    CASE(SOLVER_SPARSE_MATRICES)
      solverEquations%sparsityType=SOLVER_SPARSE_MATRICES
    CASE(SOLVER_FULL_MATRICES)
      solverEquations%sparsityType=SOLVER_FULL_MATRICES
    CASE DEFAULT
      localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverEquations_SparsityTypeSet")
    RETURN
999 ERRORSEXITS("SolverEquations_SparsityTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_SparsityTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the symmetry type for solver equations. \see OpenCMISS::Iron::cmfe_SolverEquations_SymmetryTypeSet
  SUBROUTINE SolverEquations_SymmetryTypeSet(solverEquations,symmetryType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to set the symmetry type for
    INTEGER(INTG), INTENT(IN) :: symmetryType !<The type of solver equations symmetry to be set \see SolverRoutines_SymmetryTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverEquations_SymmetryTypeSet",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)

    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) CALL FlagError("Can not set equations symmetry for a solver that has been linked.", &
      & err,error,*999)
    
    SELECT CASE(symmetryType)
    CASE(SOLVER_SYMMETRIC_MATRICES)
      solverEquations%symmetryType=SOLVER_SYMMETRIC_MATRICES
    CASE(SOLVER_FULL_MATRICES)
      solverEquations%symmetryType=SOLVER_UNSYMMETRIC_MATRICES
    CASE DEFAULT
      localError="The specified solver equations symmetry type of "// &
        & TRIM(NumberToVString(symmetryType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverEquations_SymmetryTypeSet")
    RETURN
999 ERRORSEXITS("SolverEquations_SymmetryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_SymmetryTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the time dependence type for solver equations
  SUBROUTINE SolverEquations_TimeDependenceTypeSet(solverEquations,timeDependenceType,err,error,*)

    !Argument variables
    TYPE(SolverEquationsType), POINTER :: solverEquations !<A pointer the solver equations to set the sparsity type for
    INTEGER(INTG), INTENT(IN) :: timeDependenceType !<The type of time dependence to be set \see SolverRoutines_EquationTimeDependenceTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverEquations_TimeDependenceTypeSet",err,error,*999)

    CALL SolverEquations_AssertNotFinished(solverEquations,err,error,*999)
    NULLIFY(solver)
    CALL SolverEquations_SolverGet(solverEquations,solver,err,error,*999)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not set equations time dependence for a solver that has been linked.",err,error,*999)
      
    SELECT CASE(timeDependenceType)
    CASE(SOLVER_EQUATIONS_STATIC)
      solverEquations%timeDependence=SOLVER_EQUATIONS_STATIC
    CASE(SOLVER_EQUATIONS_QUASISTATIC)
      solverEquations%timeDependence=SOLVER_EQUATIONS_QUASISTATIC
    CASE(SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC)
      solverEquations%timeDependence=SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC
    CASE(SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC)
      solverEquations%timeDependence=SOLVER_EQUATIONS_SECOND_ORDER_DYNAMIC
    CASE DEFAULT
      localError="The specified solver equations time dependence type of "// &
        & TRIM(NumberToVString(timeDependenceType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverEquations_TimeDependenceTypeSet")
    RETURN
999 ERRORSEXITS("SolverEquations_TimeDependenceTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverEquations_TimeDependenceTypeSet
        
  !
  !================================================================================================================================
  !

  !>Finalises a solver and deallocates all memory.
  RECURSIVE SUBROUTINE Solver_Finalise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_Finalise",err,error,*999)

    IF(ASSOCIATED(solver)) THEN
      solver%label=""
      CALL Solver_LinearFinalise(solver%linearSolver,err,error,*999)
      CALL Solver_NonlinearFinalise(solver%nonlinearSolver,err,error,*999)
      CALL Solver_DynamicFinalise(solver%dynamicSolver,err,error,*999)        
      CALL Solver_DAEFinalise(solver%DAESolver,err,error,*999)        
      CALL Solver_EigenproblemFinalise(solver%eigenproblemSolver,err,error,*999)
      CALL Solver_OptimiserFinalise(solver%optimiserSolver,err,error,*999)
      CALL Solver_CellMLEvaluatorFinalise(solver%cellMLEvaluatorSolver,err,error,*999)
      CALL Solver_GeometricTransformationFinalise(solver%geometricTransformationSolver,err,error,*999)
      IF(.NOT.ASSOCIATED(solver%linkingSolver)) CALL SolverEquations_Finalise(solver%solverEquations,err,error,*999)
      IF(ALLOCATED(solver%linkedSolverTypeMap)) DEALLOCATE(solver%linkedSolverTypeMap)
      IF(ALLOCATED(solver%linkedSolvers)) DEALLOCATE(solver%linkedSolvers)
      DEALLOCATE(solver)
    ENDIF 
        
    EXITS("Solver_Finalise")
    RETURN
999 ERRORSEXITS("Solver_Finalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_Finalise
  
  !
  !================================================================================================================================
  !

  !>Set the arbitrary path logical for geometric transformation solver 
  SUBROUTINE Solver_GeometricTransformationArbitraryPathSet(solver,arbitraryPath,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    LOGICAL, INTENT(IN) :: arbitraryPath !<.TRUE. if the the transformation has an arbitrary path, .FALSE. if the path is uni-directional
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver

    ENTERS("Solver_GeometricTransformationArbitraryPathSet",err,error,*999)

    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)
    
    geometricTransformationSolver%arbitraryPath=arbitraryPath
         
    EXITS("Solver_GeometricTransformationArbitraryPathSet")
    RETURN
999 ERRORS("Solver_GeometricTransformationArbitraryPathSet",err,error)
    EXITS("Solver_GeometricTransformationArbitraryPathSet")
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationArbitraryPathSet
  
  !
  !================================================================================================================================
  !

  !>Clear transformation for a geometric transformation solver 
  SUBROUTINE Solver_GeometricTransformationClear(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: incrementIdx,i
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    
    ENTERS("Solver_GeometricTransformationClear",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)
    
    geometricTransformationSolver%transformationMatrices=0.0_DP
    DO incrementIdx=1,geometricTransformationSolver%numberOfIncrements
      DO i=1,SIZE(geometricTransformationSolver%transformationMatrices,1)
        geometricTransformationSolver%transformationMatrices(i,i,incrementIdx)=1.0_DP
      ENDDO !i
    ENDDO !incrementIdx
    IF(ALLOCATED(geometricTransformationSolver%scalings)) DEALLOCATE(geometricTransformationSolver%scalings)
        
    EXITS("Solver_GeometricTransformationClear")
    RETURN    
999 ERRORSEXITS("Solver_GeometricTransformationClear",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationClear
  
  !
  !================================================================================================================================
  !

  !>Set the field and field variable type for geometric transformation solver 
  SUBROUTINE Solver_GeometricTransformationFieldSet(solver,field,variableType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    TYPE(FieldType), POINTER :: field !<A pointer to the field to transformed
    INTEGER(INTG), INTENT(IN) :: variableType !<The variable type of the field to be transformed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfGeometricComponents,i,j
    TYPE(FieldType), POINTER :: geometricField
    TYPE(FieldVariableType), POINTER :: fieldVariable,geometricVariable
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver

    ENTERS("Solver_GeometricTransformationFieldSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(field,variableType,fieldVariable,err,error,*999)
    NULLIFY(geometricField)
    CALL Field_GeometricFieldGet(field,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfGeometricComponents,err,error,*999)
    IF(ALLOCATED(geometricTransformationSolver%transformationMatrices)) &
      & CALL FlagError("Transformation matrices are already allocated for the geometric transformation solver.", &
      & err,error,*999)
    IF(geometricTransformationSolver%arbitraryPath) THEN !Allocate memory for transformation matrix at each load increment if the transformation is arbitrary at each step
      ALLOCATE(geometricTransformationSolver%transformationMatrices(numberOfGeometricComponents+1, &
        & numberOfGeometricComponents+1,geometricTransformationSolver%numberOfIncrements),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate transformation matrices for geometric transformation solver.",err,error,*999)
    ELSE !Only allocate 1 matrix if the transformation is uni-directional.
      ALLOCATE(geometricTransformationSolver%transformationMatrices(numberOfGeometricComponents+1, &
        & numberOfGeometricComponents+1,1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate transformation matrices for geometric transformation solver.",err,error,*999)
    ENDIF
    geometricTransformationSolver%transformationMatrices=0.0_DP
    ! Set all transformation matrices to be identity matrices
    DO i=1,SIZE(geometricTransformationSolver%transformationMatrices,3)
      DO j=1,numberOfGeometricComponents+1
        geometricTransformationSolver%transformationMatrices(j,j,i)=1.0_DP
      ENDDO !j
    ENDDO !i
    geometricTransformationSolver%field=>field
    geometricTransformationSolver%fieldVariableType=variableType
    geometricTransformationSolver%fieldVariable=>fieldVariable
       
    EXITS("Solver_GeometricTransformationFieldSet")
    RETURN    
999 ERRORSEXITS("Solver_GeometricTransformationFieldSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationFieldSet
  
  !
  !================================================================================================================================
  !

  !>Set the full transformation matrix for a geometric transformation at a load increment 
  SUBROUTINE Solver_GeometricTransformationMatrixSet(solver,incrementIdx,matrix,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    INTEGER(INTG), INTENT(IN) :: incrementIdx !<The load increment index
    REAL(DP), INTENT(IN) :: matrix(:,:) !<The full transformation matrix to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: size1,size2
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_GeometricTransformationMatrixSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)    
    NULLIFY(fieldVariable)
    CALL SolverGeometricTransformation_FieldVariableGet(geometricTransformationSolver,fieldVariable,err,error,*999)
    IF(incrementIdx<1.OR.incrementIdx>geometricTransformationSolver%numberOfIncrements) THEN
      localError="The specified increment index of "//TRIM(NumberToVString(incrementIdx,"*",err,error))// &
        & " is invalid. The increment index should be >= 1 and <= "// &
        & TRIM(NumberToVString(geometricTransformationSolver%numberOfIncrements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    size1=SIZE(geometricTransformationSolver%transformationMatrices(:,:,incrementIdx),1)
    size2=SIZE(geometricTransformationSolver%transformationMatrices(:,:,incrementIdx),2)
    IF(SIZE(matrix,1)<size1.OR.SIZE(matrix,2)<size2) THEN
      localError="The size of the specified transformation matrix of "// &
        & TRIM(NumberToVString(SIZE(matrix,1),"*",err,error))//"x"//TRIM(NumberToVString(SIZE(matrix,1),"*",err,error))// &
        & " is invalid. The size should be >= "// &
        & TRIM(NumberToVString(size1,"*",err,error))//"x"//TRIM(NumberToVString(size2,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    geometricTransformationSolver%transformationMatrices(1:size1,1:size2,incrementIdx)=matrix(1:size1,1:size2)
    
    EXITS("Solver_GeometricTransformationMatrixSet")
    RETURN
999 ERRORSEXITS("Solver_GeometricTransformationMatrixSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationMatrixSet
  
  !
  !================================================================================================================================
  !

  !>Set the number of load increments for geometric transformation solver 
  SUBROUTINE Solver_GeometricTransformationNumberOfLoadIncrementsSet(solver,numberOfIncrements,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    INTEGER(INTG), INTENT(IN) :: numberOfIncrements !<The number of load increments to apply the transformation
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_GeometricTransformationNumberOfLoadIncrementsSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    IF(numberOfIncrements<1) THEN
      localError="The specified number of increments of "//TRIM(NumberToVString(numberOfIncrements,"*",err,error))// &
        & " is invalid. The number of increments should be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)
    
    geometricTransformationSolver%numberOfIncrements=numberOfIncrements
        
    EXITS("Solver_GeometricTransformationNumberOfLoadIncrementsSet")
    RETURN   
999 ERRORS("Solver_GeometricTransformationNumberOfLoadIncrementsSet",err,error)
    EXITS("Solver_GeometricTransformationNumberOfLoadIncrementsSet")
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationNumberOfLoadIncrementsSet
  
  !
  !================================================================================================================================
  !

  !>Set the rotation for a geometric transformation 
  SUBROUTINE Solver_GeometricTransformationRotationSet(solver,incrementIdx,pivotPoint,axis,theta,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    INTEGER(INTG), INTENT(IN) :: incrementIdx !<The load increment index
    REAL(DP), INTENT(IN) :: pivotPoint(:) !<The pivot point to rotate about
    REAL(DP), INTENT(IN) :: axis(:) !<The axis to  to rotate around
    REAL(DP), INTENT(IN) :: theta !<The angle to rotate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfGeometricComponents,variableType
    REAL(DP) :: u,v,w,vectorLength,rotationMatrix(4,4),transformationMatrix(4,4)
    INTEGER(INTG) :: size1,size2
    TYPE(FieldType), POINTER :: field
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_GeometricTransformationRotationSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)    
    NULLIFY(field)
    CALL SolverGeometricTransformation_FieldGet(geometricTransformationSolver,field,variableType,err,error,*999)
    IF(incrementIdx<1.OR.incrementIdx>geometricTransformationSolver%numberOfIncrements) THEN
      localError="The specified increment index of "//TRIM(NumberToVString(incrementIdx,"*",err,error))// &
        & " is invalid. The increment index should be >= 1 and <= "// &
        & TRIM(NumberToVString(geometricTransformationSolver%numberOfIncrements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(incrementIdx>1.AND..NOT.geometricTransformationSolver%arbitraryPath) &
      & CALL FlagError("Rotating a field through multiple load increments must be specified through an arbitrary path.", &
      & err,error,*999) ! Due to difficulty to scale rotation
    numberOfGeometricComponents=SIZE(geometricTransformationSolver%transformationMatrices,1)-1
    !Add rotation to matrix at a specific step
    IF(SIZE(pivotPoint,1)<numberOfGeometricComponents) THEN
      localError="The size of the specified pivot point of "//TRIM(NumberToVString(SIZE(pivotPoint,1),"*",err,error))// &
        & " is invalid. The size should be >= the number of geometric components of "// &
        & TRIM(NumberToVString(numberOfGeometricComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(SIZE(axis,1)<numberOfGeometricComponents) THEN
      localError="The size of the specified axis of "//TRIM(NumberToVString(SIZE(axis,1),"*",err,error))// &
        & " is invalid. The size should be >= the number of geometric components of "// &
        & TRIM(NumberToVString(numberOfGeometricComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    SELECT CASE(numberOfGeometricComponents)
    CASE(1)
      CALL FlagError("Cannot set a rotation matrix with one geometric coordinate.",err,error,*999)
    CASE(2)
      !2D rotation
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(3)
      !3D rotation
      vectorLength=SQRT(axis(1)**2+axis(2)**2+axis(3)**2)
      u=axis(1)/vectorLength
      v=axis(2)/vectorLength
      w=axis(3)/vectorLength
      rotationMatrix=0.0_DP
      rotationMatrix(1,1)=u**2+(v**2+w**2)*COS(theta)
      rotationMatrix(1,2)=u*v*(1.0_DP-COS(theta))-w*SIN(theta)
      rotationMatrix(1,3)=u*w*(1-COS(theta))+v*SIN(theta)
      rotationMatrix(2,1)=u*v*(1-COS(theta))+w*SIN(theta)
      rotationMatrix(2,2)=v**2+(u**2+w**2)*COS(theta)
      rotationMatrix(2,3)=v*w*(1-COS(theta))-u*SIN(theta)
      rotationMatrix(3,1)=u*w*(1-COS(theta))-v*SIN(theta)
      rotationMatrix(3,2)=v*w*(1-COS(theta))+u*SIN(theta)
      rotationMatrix(3,3)=w**2+(u**2+v**2)*COS(theta)
      rotationMatrix(1,4)=(pivotPoint(1)*(v**2+w**2)-u*(pivotPoint(2)*v+pivotPoint(3)*w))*(1-COS(theta))+ &
        & (pivotPoint(2)*w-pivotPoint(3)*v)*SIN(theta)
      rotationMatrix(2,4)=(pivotPoint(2)*(u**2+w**2)-v*(pivotPoint(1)*u+pivotPoint(3)*w))*(1-COS(theta))+ &
        & (pivotPoint(3)*u-pivotPoint(1)*w)*SIN(theta)
      rotationMatrix(3,4)=(pivotPoint(3)*(u**2+v**2)-w*(pivotPoint(1)*u+pivotPoint(2)*v))*(1-COS(theta))+ &
        & (pivotPoint(1)*v-pivotPoint(2)*u)*SIN(theta)
      rotationMatrix(4,4)=1.0_DP
    CASE DEFAULT
      localError="The number of geometric components of "//TRIM(NumberToVString(numberOfGeometricComponents,"*",err,error))// &
        & " is invalid. The number of geometric components should be >=1 and <= 3."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Calculate new transformation matrix by multiplying the old matrix stored with the new rotation matrix
    transformationMatrix(1:numberOfGeometricComponents+1,1:numberOfGeometricComponents+1)= &
      & MATMUL(solver%geometricTransformationSolver%transformationMatrices(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1,incrementIdx),rotationMatrix(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1))
    ! Store the new transformation matrix
    solver%geometricTransformationSolver%transformationMatrices(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1,incrementIdx)=transformationMatrix(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1)
    
    EXITS("Solver_GeometricTransformationRotationSet")
    RETURN
999 ERRORSEXITS("Solver_GeometricTransformationRotationSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationRotationSet
  
  !
  !================================================================================================================================
  !

  !>Set the scalings for geometric transformation solver 
  SUBROUTINE Solver_GeometricTransformationScalingsSet(solver,scalings,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    REAL(DP), INTENT(IN) :: scalings(:) !<The scalings vector to set for uni-directional transformation
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_GeometricTransformationScalingsSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)    
    IF(geometricTransformationSolver%arbitraryPath) &
      & CALL FlagError("Transformation with arbitrary path can not have uni-directional scalings.",err,error,*999)
    IF(SIZE(scalings,1)<geometricTransformationSolver%numberOfIncrements) THEN
      localError="The size of the specified scalings of "//TRIM(NumberToVString(SIZE(scalings,1),"*",err,error))// &
        & " is invalid. The size should be >= the number of increments of "// &
        & TRIM(NumberToVString(geometricTransformationSolver%numberOfIncrements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    IF(ALLOCATED(geometricTransformationSolver%scalings)) DEALLOCATE(geometricTransformationSolver%scalings)
    ALLOCATE(geometricTransformationSolver%scalings(geometricTransformationSolver%numberOfIncrements),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate scalings for geometric transformation solver.",err,error,*999)
    geometricTransformationSolver%scalings(1:geometricTransformationSolver%numberOfIncrements)= &
      & scalings(1:geometricTransformationSolver%numberOfIncrements)
        
    EXITS("Solver_GeometricTransformationScalingsSet")
    RETURN   
999 ERRORSEXITS("Solver_GeometricTransformationScalingsSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationScalingsSet
  
  !
  !================================================================================================================================
  !

  !>Set the translation for a geometric transformation 
  SUBROUTINE Solver_GeometricTransformationTranslationSet(solver,incrementIdx,translation,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the field for
    INTEGER(INTG), INTENT(IN) :: incrementIdx !<The load increment index
    REAL(DP), INTENT(IN) :: translation(:) !<The translation vector to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfGeometricComponents,i,variableType
    REAL(DP) :: transformationMatrix(4,4),translationMatrix(4,4)
    TYPE(FieldType), POINTER :: field
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_GeometricTransformationTranslationSet",err,error,*999)
    
    CALL Solver_AssertIsGeometricTransformation(solver,err,error,*999)
    NULLIFY(geometricTransformationSolver)
    CALL Solver_GeometricTransformationSolverGet(solver,geometricTransformationSolver,err,error,*999)    
    NULLIFY(field)
    CALL SolverGeometricTransformation_FieldGet(geometricTransformationSolver,field,variableType,err,error,*999)
    IF(incrementIdx<1.OR.incrementIdx>geometricTransformationSolver%numberOfIncrements) THEN
      localError="The specified increment index of "//TRIM(NumberToVString(incrementIdx,"*",err,error))// &
        & " is invalid. The increment index should be >= 1 and <= "// &
        & TRIM(NumberToVString(geometricTransformationSolver%numberOfIncrements,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    numberOfGeometricComponents=SIZE(solver%geometricTransformationSolver%transformationMatrices,incrementIdx)-1
    IF(SIZE(translation,1)<numberOfGeometricComponents) THEN
      localError="The size of the specified translation vector of "//TRIM(NumberToVString(SIZE(translation,1),"*",err,error))// &
        & " is invalid. The size should be >= the number of geometric components of "// &
        & TRIM(NumberToVString(numberOfGeometricComponents,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Add translation to matrix at a specific step
    translationMatrix=0.0_DP
    transformationMatrix=0.0_DP
    DO i=1,4
      translationMatrix(i,i)=1.0_DP
    ENDDO
    translationMatrix(1:numberOfGeometricComponents,numberOfGeometricComponents+1)=translation(1:numberOfGeometricComponents)
    !Calculate the new transformation matrix by multiplying the old matrix with the new translation matrix
    transformationMatrix=MATMUL(geometricTransformationSolver%transformationMatrices(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1,incrementIdx),translationMatrix(1:numberOfGeometricComponents+1, &
      & 1:numberOfGeometricComponents+1))
    !Store the new transformation matrix
    geometricTransformationSolver%transformationMatrices(1:numberOfGeometricComponents+1,1:numberOfGeometricComponents+1, &
      & incrementIdx)=transformationMatrix(1:numberOfGeometricComponents+1,1:numberOfGeometricComponents+1)
   
    EXITS("Solver_GeometricTransformationTranslationSet")
    RETURN
999 ERRORS("Solver_GeometricTransformationTranslationSet",err,error)
    EXITS("Solver_GeometricTransformationTranslationSet")
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationTranslationSet
  
  !
  !================================================================================================================================
  !

  !>Finalise a geometric transformation solver for a solver.
  SUBROUTINE Solver_GeometricTransformationFinalise(geometricTransformationSolver,err,error,*)

    !Argument variables
    TYPE(GeometricTransformationSolverType), POINTER :: geometricTransformationSolver !<A pointer the linear solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_GeometricTransformationFinalise",err,error,*999)

    IF(ASSOCIATED(geometricTransformationSolver)) THEN   
      NULLIFY(geometricTransformationSolver%solver)
      NULLIFY(geometricTransformationSolver%field)
      geometricTransformationSolver%arbitraryPath=.FALSE.
      IF(ALLOCATED(geometricTransformationSolver%scalings)) DEALLOCATE(geometricTransformationSolver%scalings)
      IF(ALLOCATED(geometricTransformationSolver%transformationMatrices))  &
        & DEALLOCATE(geometricTransformationSolver%transformationMatrices)
      geometricTransformationSolver%numberOfIncrements=0
      NULLIFY(geometricTransformationSolver%field)
      geometricTransformationSolver%fieldVariableType=0
      NULLIFY(geometricTransformationSolver%fieldVariable)
      DEALLOCATE(geometricTransformationSolver)
    ENDIF
        
    EXITS("Solver_GeometricTransformationFinalise")
    RETURN
999 ERRORSEXITS("Solver_GeometricTransformationFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialise a geometric transformation solver for a solver.
  SUBROUTINE Solver_GeometricTransformationInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the geometric transformation solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,loopType,maximumNumberOfIterations
    TYPE(SolversType), POINTER :: solvers
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ControlLoopWhileType), POINTER :: whileLoop
    TYPE(ControlLoopLoadIncrementType), POINTER :: loadIncrementLoop
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_GeometricTransformationInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%geometricTransformationSolver)) &
      & CALL FlagError("Geometric transformation solver is already associated for this solver.",err,error,*998)
      
    !Allocate and initialise a geometric transformation solver
    ALLOCATE(solver%geometricTransformationSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver geometric transformation solver.",err,error,*999)
    solver%geometricTransformationSolver%solver=>solver
    solver%geometricTransformationSolver%arbitraryPath=.FALSE.
    !Set default number of load increment
    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    IF(loopType==CONTROL_WHILE_LOOP_TYPE) THEN
      CALL ControlLoop_MaximumNumberOfIterationsGet(controlLoop,maximumNumberOfIterations,err,error,*999)
      solver%geometricTransformationSolver%numberOfIncrements=maximumNumberOfIterations
    ELSE IF(loopType==CONTROL_LOAD_INCREMENT_LOOP_TYPE) THEN
      CALL ControlLoop_MaximumNumberOfIterationsGet(controlLoop,maximumNumberOfIterations,err,error,*999)
      solver%geometricTransformationSolver%numberOfIncrements=maximumNumberOfIterations
    ELSE
      !For other loop types set number of increment to be 1  
      solver%geometricTransformationSolver%numberOfIncrements=1
    ENDIF
    !Nullify field
    NULLIFY(solver%geometricTransformationSolver%field)
    solver%geometricTransformationSolver%fieldVariableType=0
    NULLIFY(solver%geometricTransformationSolver%fieldVariable)
       
    EXITS("Solver_GeometricTransformationInitialise")
    RETURN
999 CALL Solver_GeometricTransformationFinalise(solver%geometricTransformationSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_GeometricTransformationInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_GeometricTransformationInitialise

  !
  !================================================================================================================================
  !
  
  !>Create a CellML evaluator solver for the Newton solver
  SUBROUTINE Solver_NewtonCellMLEvaluatorCreate(solver,cellMLSolver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to create the CellML evaluator solver for
    TYPE(SolverType), POINTER :: cellMLSolver !<On return, a pointer to the created CellML evaluator solver. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(SolverType), POINTER :: linkedNonlinearSolver
            
    ENTERS("Solver_NewtonCellMLEvaluatorCreate",err,error,*999)

    IF(ASSOCIATED(cellMLSolver)) CALL FlagError("CellML solver is already associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(nonlinearSolver)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverExists(solver,dynamicSolver,err,error,*999)
    IF(ASSOCIATED(dynamicSolver)) THEN
      NULLIFY(linkedNonlinearSolver)
      CALL SolverDynamic_LinkedNonlinearSolverGet(dynamicSolver,linkedNonlinearSolver,err,error,*999)
      CALL Solver_NonlinearSolverGet(linkedNonlinearSolver,nonlinearSolver,err,error,*999)
    ELSE
      CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    ENDIF
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(ASSOCIATED(newtonSolver%cellMLEvaluatorSolver)) &
      & CALL FlagError("Newton solver CellML evaluator solver is already associated.",err,error,*999)
    ALLOCATE(newtonSolver%cellMLEvaluatorSolver,STAT=err)
    IF(err/=0) CALL FlagError("Cannot allocate CellML evaluator solver.",err,error,*999)
    NULLIFY(cellMLSolver)
    CALL SolverNonlinearNewton_LinkedCellMLSolverGet(newtonSolver,cellMLSolver,err,error,*999)
    NULLIFY(cellMLSolver%solvers)
    CALL Solver_Initialise(cellMLSolver,err,error,*999)
    CALL Solver_CellMLEvaluatorInitialise(cellMLSolver,err,error,*999)
        
    EXITS("Solver_NewtonCellMLEvaluatorCreate")
    RETURN
999 ERRORSEXITS("Solver_NewtonCellMLEvaluatorCreate",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonCellMLEvaluatorCreate

  !
  !================================================================================================================================
  !

  !>Initialise a solver for solvers
  SUBROUTINE Solvers_SolverInitialise(solvers,solverIndex,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer the solvers to initialise the solver for
    INTEGER(INTG), INTENT(IN) :: solverIndex !<The solver index in solvers to initialise the solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solvers_SolverInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*998)
    IF(solverIndex<1.OR.solverIndex>solvers%numberOfSolvers) THEN
      localError="The solver index of "//TRIM(NumberToVString(solverIndex,"*",err,error))// &
        & " is invalid. The solver index must be >= 1 and <= "// &
        & TRIM(NumberToVString(solvers%numberOfSolvers,"*",err,error))//"."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    IF(.NOT.ALLOCATED(solvers%solvers)) CALL FlagError("Solvers solvers is not allocated.",err,error,*998)    
    IF(ASSOCIATED(solvers%solvers(solverIndex)%ptr)) &
      & CALL FlagError("Solver pointer is already associated for this solver index.",err,error,*998)
          
    ALLOCATE(solvers%solvers(solverIndex)%ptr,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver.",err,error,*999)
    solvers%solvers(solverIndex)%ptr%solvers=>solvers
    CALL Solver_Initialise(solvers%solvers(solverIndex)%ptr,err,error,*999)
    solvers%solvers(solverIndex)%ptr%globalNumber=solverIndex
    !Default to a linear solver and initialise
    CALL Solver_LinearInitialise(solvers%solvers(solverIndex)%ptr,err,error,*999)
    solvers%solvers(solverIndex)%ptr%solveType=SOLVER_LINEAR_TYPE
        
    EXITS("Solvers_SolverInitialise")
    RETURN
999 CALL Solver_Finalise(solvers%solvers(solverIndex)%ptr,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solvers_SolverInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solvers_SolverInitialise

  !
  !================================================================================================================================
  !

  !>Initialise a solver 
  SUBROUTINE Solver_Initialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverIdx
    
    ENTERS("Solver_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(solver%linkingSolver)
    ALLOCATE(solver%linkedSolverTypeMap(SOLVER_NUMBER_OF_SOLVER_TYPES),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate linked solver type map.",err,error,*999)
    DO solverIdx=1,SOLVER_NUMBER_OF_SOLVER_TYPES
      NULLIFY(solver%linkedSolverTypeMap(solverIdx)%ptr)
    ENDDO !solverIdx
    solver%numberOfLinkedSolvers=0
    solver%solverFinished=.FALSE.
    solver%label=""
    solver%outputType=SOLVER_NO_OUTPUT
    NULLIFY(solver%linearSolver)
    NULLIFY(solver%nonlinearSolver)
    NULLIFY(solver%dynamicSolver)
    NULLIFY(solver%DAESolver)
    NULLIFY(solver%eigenproblemSolver)
    NULLIFY(solver%optimiserSolver)
    NULLIFY(solver%cellMLEvaluatorSolver)
    NULLIFY(solver%solverEquations)
    NULLIFY(solver%cellMLEquations)
    NULLIFY(solver%geometricTransformationSolver)
    NULLIFY(solver%workGroup)
    
    EXITS("Solver_Initialise")
    RETURN
999 ERRORSEXITS("Solver_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_Initialise

  !
  !================================================================================================================================
  !

  !>Sets the label of a solver. \see OpenCMISS::Iron::cmfe_Solver_LabelSet
  SUBROUTINE Solver_LabelSetC(solver,label,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to set the label for 
    CHARACTER(LEN=*), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_LabelSetC",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    solver%label=label
     
    EXITS("Solver_LabelSetC")
    RETURN
999 ERRORSEXITS("Solver_LabelSetC",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LabelSetC

  !
  !================================================================================================================================
  !

  !>Sets the label of a solver. \see OpenCMISS::Iron::cmfe_Solver_LabelSet
  SUBROUTINE Solver_LabelSetVS(solver,label,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to set the label for 
    TYPE(VARYING_STRING), INTENT(IN) :: label !<The label to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_LabelSetVS",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    solver%label=label
    
    EXITS("Solver_LabelSetVS")
    RETURN
999 ERRORSEXITS("Solver_LabelSetVS",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LabelSetVS

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library type to use for the solver. \see OpenCMISS::Iron::cmfe_Solver_LibraryTypeSet
  SUBROUTINE Solver_LibraryTypeSet(solver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the type of
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library to use for the solver \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_LibraryTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    
    SELECT CASE(solver%solveType)
    CASE(SOLVER_LINEAR_TYPE)
      CALL SolverLinear_LibraryTypeSet(solver%linearSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL SolverNonlinear_LibraryTypeSet(solver%nonlinearSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL SolverDynamic_LibraryTypeSet(solver%dynamicSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_DAE_TYPE)
      CALL SolverDAE_LibraryTypeSet(solver%daeSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL SolverEigenproblem_LibraryTypeSet(solver%eigenproblemSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_OPTIMISER_TYPE)
      CALL SolverOptimiser_LibraryTypeSet(solver%optimiserSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_CELLML_EVALUATOR_TYPE)
      CALL SolverCellMLEvaluator_LibraryTypeSet(solver%cellMLEvaluatorSolver,solverLibraryType,err,error,*999)
    CASE DEFAULT
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("Solver_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("Solver_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LibraryTypeSet
  
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear solver 
  SUBROUTINE SolverLinear_CreateFinish(linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: nonlinearLibraryType
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: linkingSolver,solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinear_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    linkingSolver=>solver%linkingSolver
    IF(ASSOCIATED(linkingSolver)) THEN
      IF(linkingSolver%solveType==SOLVER_NONLINEAR_TYPE) THEN
        NULLIFY(nonlinearSolver)
        CALL Solver_NonlinearSolverGet(linkingSolver,nonlinearSolver,err,error,*999)
        CALL SolverNonlinear_LibraryTypeGet(nonlinearSolver,nonlinearLibraryType,err,error,*999)
        linearSolver%linkedNewtonPetSCSolver=(nonlinearLibraryTYpe==SOLVER_PETSC_LIBRARY)
      ENDIF
    ENDIF
    SELECT CASE(linearSolver%linearSolveType)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      CALL SolverLinearDirect_CreateFinish(linearSolver%directSolver,err,error,*999)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      CALL SolverLinearIterative_CreateFinish(linearSolver%iterativeSolver,err,error,*999)
    CASE DEFAULT
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverLinear_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverLinear_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a Cholesky direct linear solver and deallocate all memory.
  SUBROUTINE SolverLinearDirect_CholeskyFinalise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the linear direct solver to finalise the Cholesky solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinearDirect_CholeskyFinalise",err,error,*999)

    IF(ASSOCIATED(directSolver)) THEN
      CALL FlagError("Not implemented.",err,error,*999)
    ENDIF

    EXITS("SolverLinearDirect_CholeskyFinalise")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_CholeskyFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_CholeskyFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Cholesky direct linear solver for a direct linear solver.
  SUBROUTINE SolverLinearDirect_CholeskyInitialise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer the direct linear solver to initialise the Cholesky direct linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverLinearDirect_CholeskyInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverLinearDirect_CholeskyInitialise")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_CholeskyInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_CholeskyInitialise

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear direct solver 
  SUBROUTINE SolverLinearDirect_CreateFinish(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the linear direct solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: groupCommunicator,numberOfMatrices,sparsityType,symmetryType
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(SolverType), POINTER :: linkingSolver,solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinearDirect_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Linear direct solver is not associated.",err,error,*999)
    NULLIFY(linearSolver)
    CALL SolverLinearDirect_LinearSolverGet(directSolver,linearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    SELECT CASE(directSolver%directSolverType)
    CASE(SOLVER_DIRECT_LU)
      linkingSolver=>solver%linkingSolver
      IF(ASSOCIATED(linkingSolver)) THEN
        !Matrices have already been set up by linking solver
        SELECT CASE(directSolver%solverLibrary)
        CASE(SOLVER_CMISS_LIBRARY) !All non-PETSc libraries
          CALL FlagError("Non-PETSc linear solver cannot be linked to PETSc nonlinear solver.",err,error,*999)
        END SELECT
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(linkingSolver,solverEquations,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
        solverEquations=>solver%linkingSolver%solverEquations
      ELSE
        !Set up solver matrices
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Create the solver matrices
        NULLIFY(solverMatrices)
        CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
        CALL SolverMatrices_LibraryTypeSet(solverMatrices,directSolver%solverLibrary,err,error,*999)
        CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
        SELECT CASE(sparsityType)
        CASE(SOLVER_SPARSE_MATRICES)
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
        CASE(SOLVER_FULL_MATRICES)
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
        CASE DEFAULT
          localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
        SELECT CASE(symmetryType)
        CASE(SOLVER_SYMMETRIC_MATRICES)
          CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
        CASE(SOLVER_UNSYMMETRIC_MATRICES)
          CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
        CASE DEFAULT
          localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      ENDIF
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
      IF(numberOfMatrices/=1) THEN
        localError="The given number of solver matrices of "//TRIM(NumberToVString(numberOfMatrices,"*",err,error))// &
          & " is invalid. There should only be one solver matrix for a linear direct solver."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
      NULLIFY(distributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)

      !Set up direct solver
      SELECT CASE(directSolver%solverLibrary)
      CASE(SOLVER_CMISS_LIBRARY)
        !Nothing else to do
      CASE(SOLVER_MUMPS_LIBRARY,SOLVER_SUPERLU_LIBRARY,SOLVER_PASTIX_LIBRARY,SOLVER_LAPACK_LIBRARY)
        !Set up solver through PETSc
        CALL PETSc_KSPCreate(groupCommunicator,directSolver%ksp,err,error,*999)
        
        !Set any further KSP options from the command line options
        CALL PETSc_KSPSetFromOptions(directSolver%ksp,err,error,*999)
        !Set the solver matrix to be the KSP matrix
        NULLIFY(petscMatrix)
        CALL DistributedMatrix_PetSCMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
        CALL PETSc_KSPSetOperators(directSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
        !Check that the solver supports the matrix sparsity type
        CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
        SELECT CASE(sparsityType)
        CASE(SOLVER_FULL_MATRICES)
          SELECT CASE(directSolver%solverLibrary)
          CASE(SOLVER_MUMPS_LIBRARY,SOLVER_SUPERLU_LIBRARY,SOLVER_PASTIX_LIBRARY)
            CALL FlagError("Solver library does not support full matrices. Please use sparse matrices "// &
              & "or select the LAPACK library type for the linear direct solver.",err,error,*999)
          END SELECT
        CASE(SOLVER_SPARSE_MATRICES)
          SELECT CASE(directSolver%solverLibrary)
          CASE(SOLVER_LAPACK_LIBRARY)
            CALL FlagError("Solver library does not support sparse matrices. Please use full matrices "// &
              & "or select another solver library type for the linear direct solver.",err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Set the KSP type to preonly
        CALL PETSc_KSPSetType(directSolver%ksp,PETSC_KSPPREONLY,err,error,*999)
        !Get the pre-conditioner
        CALL PETSc_KSPGetPC(directSolver%ksp,directSolver%pc,err,error,*999)
        !Set the PC type to LU
        CALL PETSc_PCSetType(directSolver%pc,PETSC_PCLU,err,error,*999)
        SELECT CASE(directSolver%solverLibrary)
        CASE(SOLVER_MUMPS_LIBRARY)
          !Set the PC factorisation package to MUMPS
          CALL PETSc_PCFactorSetMatSolverPackage(directSolver%pc,PETSC_MAT_SOLVER_MUMPS,err,error,*999)
        CASE(SOLVER_SUPERLU_LIBRARY)
          !Set the PC factorisation package to SuperLU_DIST
          CALL PETSc_PCFactorSetMatSolverPackage(directSolver%pc,PETSC_MAT_SOLVER_SUPERLU_DIST,err,error,*999)
        CASE(SOLVER_LAPACK_LIBRARY)
          CALL FlagError("LAPACK not available in this version of PETSc.",err,error,*999)
        CASE(SOLVER_PASTIX_LIBRARY)
          !Set the PC factorisation package to PaStiX
          CALL PETSc_PCFactorSetMatSolverPackage(directSolver%pc,PETSC_MAT_SOLVER_PASTIX,err,error,*999)
        END SELECT
      CASE(SOLVER_SPOOLES_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_UMFPACK_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_LUSOL_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_ESSL_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(directSolver%solverLibrary,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DIRECT_CHOLESKY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_DIRECT_SVD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The direct solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverLinearDirect_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinearDirect_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a direct linear solver for a linear solver and deallocate all memory.
  SUBROUTINE SolverLinear_DirectFinalise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the lienar direct solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver

    ENTERS("SolverLinear_DirectFinalise",err,error,*999)

    IF(ASSOCIATED(directSolver)) THEN
      linearSolver=>directSolver%linearSolver
      IF(ASSOCIATED(linearSolver)) THEN
        IF(.NOT.linearSolver%linkedNewtonPetSCSolver) CALL SolverLinearDirect_LUFinalise(directSolver,err,error,*999)
      ENDIF
      DEALLOCATE(directSolver)
    ENDIF

    EXITS("SolverLinear_DirectFinalise")
    RETURN
999 ERRORSEXITS("SolverLinear_DirectFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_DirectFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a direct linear solver for a lienar solver
  SUBROUTINE SolverLinear_DirectInitialise(linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer the linear solver to initialise the direct solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverLinear_DirectInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*998)
    IF(ASSOCIATED(linearSolver%directSolver)) &
      & CALL FlagError("Direct solver is already associated for this linear solver.",err,error,*998)
      
    ALLOCATE(linearSolver%directSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate linear solver direct solver.",err,error,*999)
    linearSolver%directSolver%linearSolver=>linearSolver
    !Default to an LU direct linear solver
    linearSolver%directSolver%directSolverType=SOLVER_DIRECT_LU
    CALL SolverLinearDirect_LUInitialise(linearSolver%directSolver,err,error,*999)
       
    EXITS("SolverLinear_DirectInitialise")
    RETURN
999 CALL SolverLinear_DirectFinalise(linearSolver%directSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverLinear_DirectInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_DirectInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a direct linear solver.
  SUBROUTINE SolverLinearDirect_LibraryTypeSet(directSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer the direct linear solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the direct linear solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinearDirect_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*999)
    
    SELECT CASE(directSolver%directSolverType)
    CASE(SOLVER_DIRECT_LU)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemeted.",err,error,*999)
      CASE(SOLVER_MUMPS_LIBRARY)
        directSolver%solverLibrary=SOLVER_MUMPS_LIBRARY
      CASE(SOLVER_SUPERLU_LIBRARY)
        directSolver%solverLibrary=SOLVER_SUPERLU_LIBRARY
      CASE(SOLVER_SPOOLES_LIBRARY)
        CALL FlagError("Not implemeted.",err,error,*999)
      CASE(SOLVER_LUSOL_LIBRARY)
        CALL FlagError("Not implemeted.",err,error,*999)
      CASE(SOLVER_ESSL_LIBRARY)
        CALL FlagError("Not implemeted.",err,error,*999)
      CASE(SOLVER_LAPACK_LIBRARY)
        directSolver%solverLibrary=SOLVER_LAPACK_LIBRARY
      CASE(SOLVER_PASTIX_LIBRARY)
        directSolver%solverLibrary=SOLVER_PASTIX_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a LU direct linear solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DIRECT_CHOLESKY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_DIRECT_SVD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The direct linear solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverLinearDirect_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Finalise a LU direct linear solver and deallocate all memory.
  SUBROUTINE SolverLinearDirect_LUFinalise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the linear direct solver to finalise the LU solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinearDirect_LUFinalise",err,error,*999)

    IF(ASSOCIATED(directSolver)) THEN
      SELECT CASE(directSolver%solverLibrary)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_MUMPS_LIBRARY)
        !Call MUMPS through PETSc
        CALL PETSc_PCFinalise(directSolver%pc,err,error,*999)
        CALL PETSc_KSPFinalise(directSolver%ksp,err,error,*999)
      CASE(SOLVER_SUPERLU_LIBRARY)
        !Call SuperLU through PETSc
        CALL PETSc_PCFinalise(directSolver%pc,err,error,*999)
        CALL PETSc_KSPFinalise(directSolver%ksp,err,error,*999)
      CASE(SOLVER_SPOOLES_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_UMFPACK_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_LUSOL_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_ESSL_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_LAPACK_LIBRARY)
        !Call SuperLU through PETSc
        CALL PETSc_PCFinalise(directSolver%pc,err,error,*999)
        CALL PETSc_KSPFinalise(directSolver%ksp,err,error,*999)
      CASE(SOLVER_PASTIX_LIBRARY)
        !Call PaStiX through PETSc
        CALL PETSc_PCFinalise(directSolver%pc,err,error,*999)
        CALL PETSc_KSPFinalise(directSolver%ksp,err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(directSolver%solverLibrary,"*",err,error))// &
          & " is invalid for a LU direct linear solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF

    EXITS("SolverLinearDirect_LUFinalise")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_LUFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_LUFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a LU direct linear solver for a direct linear solver.
  SUBROUTINE SolverLinearDirect_LUInitialise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer the direct linear solver to initialise the LU direct linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("SolverLinearDirect_LUInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*998)
    
    !Default to MUMPS library
    directSolver%solverLibrary=SOLVER_MUMPS_LIBRARY
    !Call MUMPS through PETSc
    CALL PETSc_PCInitialise(directSolver%pc,err,error,*999)
    CALL PETSc_KSPInitialise(directSolver%ksp,err,error,*999)
       
    EXITS("SolverLinearDirect_LUInitialise")
    RETURN
999 CALL SolverLinearDirect_LUFinalise(directSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverLinearDirect_LUInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_LUInitialise

  !
  !================================================================================================================================
  !

!!\todo Allow for the mumps parameters to be set during the solver creation (i.e., cache and defer setting until we have PETSc matrix)
  
  !>Sets MUMPS ICNTL(icntl)=ivalue through PETSc Mat API (see MUMPS user guide for more info). Must be called after the boundary conditions have been set up.
  SUBROUTINE Solver_MumpsSetIcntl(solver,icntl,ivalue,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS ICNTL integer control parameter 
    INTEGER(INTG), INTENT(IN) :: ivalue !<The MUMPS ICNTL integer value to set: ICNTL(icntl)=ivalue
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfMatrices
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(SolverType), POINTER :: linkingSolver
    TYPE(SolverEquationsType), POINTER :: linkingSolverEquations,solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(PetscMatType) :: petscFactoredMatrix !<The factored matrix obtained by calling MatGetFactor() from PETSc-MUMPS interface 
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_MumpsSetIcntl",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsDirect(linearSolver,err,error,*999)
    NULLIFY(directSolver)
    CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
    SELECT CASE(directSolver%directSolverType)
    CASE(SOLVER_DIRECT_LU)
      SELECT CASE(directSolver%solverLibrary)
      CASE(SOLVER_MUMPS_LIBRARY)        
        solverEquations=>solver%solverEquations
        NULLIFY(solverMatrices)
        IF(ASSOCIATED(solverEquations)) THEN
          !Solver equations for this solver.
          CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
        ELSE
          !No solver equations. See if there are solver equations in the linking solver.
          linkingSolver=>solver%linkingSolver
          IF(.NOT.ASSOCIATED(linkingSolver)) &
            & CALL FlagError("Solver matrices could not be found for solver or linking solver.",err,error,*999)
          NULLIFY(linkingSolverEquations)
          CALL Solver_SolverEquationsGet(linkingSolver,linkingSolverEquations,err,error,*999)
          CALL SolverEquations_SolverMatricesGet(linkingSolverEquations,solverMatrices,err,error,*999)
        ENDIF
        CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
        IF(numberOfMatrices/=1) THEN
          localError="The given number of solver matrices of "//TRIM(NumberToVstring(numberOfMatrices,"*",err,error))// &
            & " is invalid. There should only be one solver matrix for a linear direct solver."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
        NULLIFY(distributedMatrix)
        CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
        NULLIFY(petscMatrix)
        CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
        !Call MatGetFactor to create matrix petscFactoredMatrix from preconditioner context
        CALL PETSc_PCFactorSetUpMatSolverPackage(directSolver%pc,err,error,*999)
        CALL PETSc_PCFactorGetMatrix(directSolver%pc,petscFactoredMatrix,err,error,*999)
        !Set ICNTL(icntl)=ivalue
        CALL PETSc_MatMumpsSetIcntl(petscFactoredMatrix,icntl,ivalue,err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVstring(directSolver%solverLibrary,"*",err,error))// &
          & " is invalid. Use the MUMPS library when calling Solver_MumpsSetIcntl"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DIRECT_CHOLESKY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_DIRECT_SVD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The direct solver type of "//TRIM(NumberToVstring(directSolver%directSolverType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_MumpsSetIcntl")
    RETURN
999 ERRORSEXITS("Solver_MumpsSetIcntl",err,error)    
    RETURN 1
    
  END SUBROUTINE Solver_MumpsSetIcntl

  !
  !================================================================================================================================
  !

!!\todo Allow for the mumps parameters to be set during the solver creation (i.e., cache and defer setting until we have PETSc matrix)

  !>Sets MUMPS CNTL(icntl)=val through PETSc Mat API (see MUMPS user guide for more info). Must be called after the boundary conditions have been set up.
  SUBROUTINE Solver_MumpsSetCntl(solver,icntl,val,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(IN) :: icntl !<The MUMPS ICNTL integer control parameter 
    REAL(DP), INTENT(IN) :: val !<The MUMPS CNTL real value to set: CNTL(icntl)=val
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfMatrices
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(PetscMatType) :: petscFactoredMatrix !<The factored matrix obtained by calling MatGetFactor() from PETSc-MUMPS interface 
    TYPE(SolverType), POINTER :: linkingSolver
    TYPE(SolverEquationsType), POINTER :: linkingSolverEquations,solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_MumpsSetCntl",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsDirect(linearSolver,err,error,*999)
    NULLIFY(directSolver)
    CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
    SELECT CASE(directSolver%directSolverType)
    CASE(SOLVER_DIRECT_LU)
      SELECT CASE(directSolver%solverLibrary)
      CASE(SOLVER_MUMPS_LIBRARY)
        solverEquations=>solver%solverEquations
        NULLIFY(solverMatrices)
        IF(ASSOCIATED(solverEquations)) THEN
          !Solver equations for this solver.
          CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
        ELSE
          !No solver equations. See if there are solver equations in the linking solver.
          linkingSolver=>solver%linkingSolver
          IF(.NOT.ASSOCIATED(linkingSolver)) &
            & CALL FlagError("Solver matrices could not be found for solver or linking solver.",err,error,*999)
          NULLIFY(linkingSolverEquations)
          CALL Solver_SolverEquationsGet(linkingSolver,linkingSolverEquations,err,error,*999)
          CALL SolverEquations_SolverMatricesGet(linkingSolverEquations,solverMatrices,err,error,*999)
        ENDIF
        CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
        IF(numberOfMatrices/=1) THEN
          localError="The given number of solver matrices of "//TRIM(NumberToVstring(numberOfMatrices,"*",err,error))// &
            & " is invalid. There should only be one solver matrix for a linear direct solver."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
        NULLIFY(distributedMatrix)
        CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
        NULLIFY(petscMatrix)
        CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
        !Call MatGetFactor to create matrix petscFactoredMatrix from preconditioner context
        CALL PETSc_PCFactorSetUpMatSolverPackage(directSolver%pc,err,error,*999)
        CALL PETSc_PCFactorGetMatrix(directSolver%pc,petscFactoredMatrix,err,error,*999)
        !Set CNTL(icntl)=val
        CALL PETSc_MatMumpsSetCntl(petscFactoredMatrix,icntl,val,err,error,*999)
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(directSolver%solverLibrary,"*",err,error))// &
          & " is invalid. Use the MUMPS library when calling Solver_MumpsSetCntl"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_DIRECT_CHOLESKY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_DIRECT_SVD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The direct solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_MumpsSetCntl")
    RETURN
999 ERRORSEXITS("Solver_MumpsSetCntl",err,error)    
    RETURN 1
    
  END SUBROUTINE Solver_MumpsSetCntl

  !
  !================================================================================================================================
  !

  !>Solve a linear direct solver 
  SUBROUTINE SolverLinearDirect_Solve(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the linear direct solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: globalRow,localRowIdx,numberOfMatrices,numberOfRows,storageType
    REAL(DP) :: solverValue,matrixValue
    REAL(DP), POINTER :: rhsData(:)
    LOGICAL :: updateMatrix
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(DistributedVectorType), POINTER :: rhsVector,solverVector
    TYPE(DistributedVectorPETScType), POINTER :: rhsPETScVector,solverPETScVector
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinearDirect_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Linear direct solver is not associated.",err,error,*999)
    NULLIFY(linearSolver)
    CALL SolverLinearDirect_LinearSolverGet(directSolver,linearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
    IF(numberOfMatrices/=1) THEN
      localError="The number of solver matrices of "//TRIM(NumberToVString(numberOfMatrices,"*",err,error))// &
        & " is invalid. There should only be one solver matrix for a linear direct solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(solverMatrix)
    CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
    CALL SolverMatrix_UpdateMatrixGet(solverMatrix,updateMatrix,err,error,*999)
    NULLIFY(distributedMatrix)
    CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
    NULLIFY(rhsVector)
    CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,rhsVector,err,error,*999)
    NULLIFY(solverVector)
    CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
    CALL DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*999)
    IF(storageType==DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(rowDOFsMapping)
      CALL SolverMapping_RowDOFSMappingGet(solverMapping,rowDOFsMapping,err,error,*999)
      CALL DistributedVector_DataGet(rhsVector,rhsData,err,error,*999)
      CALL SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*999)
      DO localRowIdx=1,numberOfRows
        CALL DomainMapping_LocalToGlobalGet(rowDOFsMapping,localRowIdx,globalRow,err,error,*999)
        CALL DistributedMatrix_ValuesGet(distributedMatrix,localRowIdx,globalRow,matrixValue,err,error,*999)
        IF(ABS(matrixValue)<=ZERO_TOLERANCE) THEN
          localError="The linear solver matrix has a zero pivot at local row "// &
            & TRIM(NumberToVString(localRowIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        solverValue=rhsData(localRowIdx)/matrixValue
        CALL DistributedVector_ValuesSet(solverVector,localRowIdx,solverValue,err,error,*999)
      ENDDO !localRowIdx
      CALL DistributedVector_DataRestore(rhsVector,rhsData,err,error,*999)
    ELSE
      SELECT CASE(directSolver%directSolverType)
      CASE(SOLVER_DIRECT_LU)
        SELECT CASE(directSolver%solverLibrary)
        CASE(SOLVER_CMISS_LIBRARY)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(SOLVER_MUMPS_LIBRARY)
          !Call MUMPS through PETSc
          NULLIFY(rhsPETScVector)
          CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
          NULLIFY(solverPETScVector)
          CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
          NULLIFY(petscMatrix)
          CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
          IF(updateMatrix) THEN
            CALL PETSc_KSPSetOperators(directSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
          ELSE
            CALL PETSc_PCSetReusePreconditioner(directSolver%pc,.TRUE.,err,error,*999)
          ENDIF
          !Solve the linear system
          CALL PETSc_KSPSolve(directSolver%ksp,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999) 
        CASE(SOLVER_SUPERLU_LIBRARY)
          !Call SuperLU through PETSc
          NULLIFY(rhsPETScVector)
          CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
          NULLIFY(solverPETScVector)
          CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
          NULLIFY(petscMatrix)
          CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
          IF(updateMatrix) THEN
            CALL PETSc_KSPSetOperators(directSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
          ELSE
            CALL PETSc_PCSetReusePreconditioner(directSolver%pc,.TRUE.,err,error,*999)
          ENDIF
          !Solve the linear system
          CALL PETSc_KSPSolve(directSolver%ksp,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999) 
        CASE(SOLVER_SPOOLES_LIBRARY)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(SOLVER_UMFPACK_LIBRARY)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(SOLVER_LUSOL_LIBRARY)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(SOLVER_ESSL_LIBRARY)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(SOLVER_LAPACK_LIBRARY)
          !Call LAPACK through PETSc
          NULLIFY(rhsPETScVector)
          CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
          NULLIFY(solverPETScVector)
          CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
          NULLIFY(petscMatrix)
          IF(updateMatrix) THEN
            CALL PETSc_KSPSetOperators(directSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
          ELSE
            CALL PETSc_PCSetReusePreconditioner(directSolver%pc,.TRUE.,err,error,*999)
          ENDIF
          !Solve the linear system
          CALL PETSc_KSPSolve(directSolver%ksp,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999) 
        CASE(SOLVER_PASTIX_LIBRARY)
          !Call PASTIX through PETSc
          NULLIFY(rhsPETScVector)
          CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
          NULLIFY(solverPETScVector)
          CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
          NULLIFY(petscMatrix)
          CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
          IF(updateMatrix) THEN
            CALL PETSc_KSPSetOperators(directSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
          ELSE
            CALL PETSc_PCSetReusePreconditioner(directSolver%pc,.TRUE.,err,error,*999)
          ENDIF
          !Solve the linear system
          CALL PETSc_KSPSolve(directSolver%ksp,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999) 
        CASE DEFAULT
          localError="The solver library type of "//TRIM(NumberToVString(directSolver%solverLibrary,"*",err,error))// &
            & " is invalid for a LU direct linear solver."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(SOLVER_DIRECT_CHOLESKY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_DIRECT_SVD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The direct linear solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
        
    EXITS("SolverLinearDirect_Solve")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinearDirect_Solve
        
  !
  !================================================================================================================================
  !

  !>Finalise a SVD direct linear solver and deallocate all memory.
  SUBROUTINE SolverLinearDirect_SVDFinalise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer to the linear direct solver to finalise the SVD solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverLinearDirect_SVDFinalise",err,error,*999)

    IF(ASSOCIATED(directSolver)) THEN
      CALL FlagError("Not implemented.",err,error,*999)
    ENDIF

    EXITS("SolverLinearDirect_SVDFinalise")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_SVDFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_SVDFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a SVD direct linear solver for a direct linear solver.
  SUBROUTINE SolverLinearDirect_SVDInitialise(directSolver,err,error,*)

    !Argument variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver !<A pointer the direct linear solver to initialise the SVD direct linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverLinearDirect_SVDInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(directSolver)) CALL FlagError("Direct linear solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
        
    EXITS("SolverLinearDirect_SVDInitialise")
    RETURN
999 ERRORSEXITS("SolverLinearDirect_SVDInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearDirect_SVDInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of direct linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearDirectTypeSet
  SUBROUTINE Solver_LinearDirectTypeSet(solver,directSolverType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the direct linear solver type for.
    INTEGER(INTG), INTENT(IN) :: directSolverType !<The type of direct linear solver to set \see SolverRoutines_DirectLinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearDirectTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsDirect(linearSolver,err,error,*999)
    NULLIFY(directSolver)
    CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
    IF(directSolverType/=directSolver%directSolverType) THEN
      !Finalise the old direct solver
      SELECT CASE(directSolver%solverLibrary)
      CASE(SOLVER_DIRECT_LU)
        CALL SolverLinearDirect_LUFinalise(directSolver,err,error,*999)
      CASE(SOLVER_DIRECT_CHOLESKY)
        CALL SolverLinearDirect_CholeskyFinalise(directSolver,err,error,*999)
      CASE(SOLVER_DIRECT_SVD)
        CALL SolverLinearDirect_SVDFinalise(directSolver,err,error,*999)
      CASE DEFAULT
        localError="The direct solver type of "//TRIM(NumberToVString(directSolver%directSolverType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Initialise the new library
      SELECT CASE(directSolverType)
      CASE(SOLVER_DIRECT_LU)
        CALL SolverLinearDirect_LUInitialise(directSolver,err,error,*999)
      CASE(SOLVER_DIRECT_CHOLESKY)
        CALL SolverLinearDirect_CholeskyInitialise(directSolver,err,error,*999)
      CASE(SOLVER_DIRECT_SVD)
        CALL SolverLinearDirect_SVDInitialise(directSolver,err,error,*999)
      CASE DEFAULT
        localError="The direct solver type of "//TRIM(NumberToVString(directSolverType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_LinearDirectTypeSet")
    RETURN
999 ERRORSEXITS("Solver_LinearDirectTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearDirectTypeSet
        
  !
  !================================================================================================================================
  !

  !>Finalise a linear solver for a solver.
  SUBROUTINE Solver_LinearFinalise(linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer the linear solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_LinearFinalise",err,error,*999)

    IF(ASSOCIATED(linearSolver)) THEN
      CALL SolverLinear_DirectFinalise(linearSolver%directSolver,err,error,*999)
      CALL SolverLinear_IterativeFinalise(linearSolver%iterativeSolver,err,error,*999)
      DEALLOCATE(linearSolver)
    ENDIF
        
    EXITS("Solver_LinearFinalise")
    RETURN
999 ERRORSEXITS("Solver_LinearFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a linear solver for a solver.
  SUBROUTINE Solver_LinearInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the linear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_LinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%linearSolver)) CALL FlagError("Linear solver is already associated for this solver.",err,error,*998)
      
    !Allocate and initialise a linear solver
    ALLOCATE(solver%linearSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver linear solver.",err,error,*999)
    solver%linearSolver%solver=>solver
    solver%linearSolver%linkedNewtonPetSCSolver=.FALSE.
    NULLIFY(solver%linearSolver%directSolver)
    NULLIFY(solver%linearSolver%iterativeSolver)
    !Default to an iterative solver
    solver%linearSolver%linearSolveType=SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
    CALL SolverLinear_IterativeInitialise(solver%linearSolver,err,error,*999)
         
    EXITS("Solver_LinearInitialise")
    RETURN
999 CALL Solver_LinearFinalise(solver%linearSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_LinearInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for an iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeAbsoluteToleranceSet
  SUBROUTINE Solver_LinearIterativeAbsoluteToleranceSet(solver,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set 
    REAL(DP), INTENT(IN) :: absoluteTolerance !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeAbsoluteToleranceSet",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(absoluteTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified absolute tolerance of "//TRIM(NumberToVString(absoluteTolerance,"*",err,error))// &
        & " is invalid. The absolute tolerance must be > 0.0"
      CALL FlagError(localError,err,error,*999)
    ENDIF
    iterativeSolver%absoluteTolerance=absoluteTolerance
    
    EXITS("Solver_LinearIterativeAbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeAbsoluteToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeAbsoluteToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a linear iterative solver 
  SUBROUTINE SolverLinearIterative_CreateFinish(iterativeSolver,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer to the linear iterative solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: groupCommunicator,sparsityType,symmetryType
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: petscMatrix
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: newtonLinesearchSolver
    TYPE(NewtonTrustregionSolverType), POINTER :: newtonTrustregionSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: quasiNewtonLinesearchSolver
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: quasiNewtonTrustregionSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: linkingSolver,solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinearIterative_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Linear iterative solver is not associated.",err,error,*999)
    NULLIFY(linearSolver)
    CALL SolverLinearIterative_LinearSolverGet(iterativeSolver,linearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    !Should really check iterative types here and then the solver library but as they are all PETSc for now hold off.
    SELECT CASE(iterativeSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      NULLIFY(solverEquations)
      NULLIFY(linkingSolver)
      CALL Solver_LinkingSolverExists(solver,linkingSolver,err,error,*999)
      IF(ASSOCIATED(linkingSolver)) THEN
        CALL Solver_SolverEquationsGet(linkingSolver,solverEquations,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      ELSE
        !Create the solver matrices and vectors
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMatrices)
        CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
        CALL SolverMatrices_LibraryTypeSet(solverMatrices,SOLVER_PETSC_LIBRARY,err,error,*999)
        CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
        SELECT CASE(sparsityType)
        CASE(SOLVER_SPARSE_MATRICES)
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
        CASE(SOLVER_FULL_MATRICES)
          CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
        CASE DEFAULT
          localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
        SELECT CASE(symmetryType)
        CASE(SOLVER_SYMMETRIC_MATRICES)
          CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
        CASE(SOLVER_UNSYMMETRIC_MATRICES)
          CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
        CASE DEFAULT
          localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      ENDIF
      !Create the PETSc KSP solver
      IF(linearSolver%linkedNewtonPetSCSolver) THEN
        IF(.NOT.ASSOCIATED(linkingSolver)) CALL FlagError("Linking solver is not associated.",err,error,*999)
        NULLIFY(nonlinearSolver)
        CALL Solver_NonlinearSolverGet(linkingSolver,nonlinearSolver,err,error,*999)
        IF(nonlinearSolver%nonlinearSolveType==SOLVER_NONLINEAR_NEWTON) THEN
          NULLIFY(newtonSolver)
          CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
          SELECT CASE(newtonSolver%newtonSolveType)
          CASE(SOLVER_NEWTON_LINESEARCH)
            NULLIFY(newtonLinesearchSolver)
            CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,newtonLinesearchSolver,err,error,*999)
            CALL PETSc_SNESGetKSP(newtonLinesearchSolver%snes,iterativeSolver%ksp,err,error,*999)
          CASE(SOLVER_NEWTON_TRUSTREGION)
            NULLIFY(newtonTrustregionSolver)
            CALL SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,newtonTrustregionSolver,err,error,*999)
            CALL PETSc_SNESGetKSP(newtonTrustregionSolver%snes,iterativeSolver%ksp,err,error,*999)
          CASE DEFAULT
            localError="The Newton solve type of "// &
              & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//"is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE IF(nonlinearSolver%nonlinearSolveType==SOLVER_NONLINEAR_QUASI_NEWTON) THEN
          NULLIFY(quasiNewtonSolver)
          CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
          SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
          CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
            NULLIFY(quasiNewtonLinesearchSolver)
            CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,quasiNewtonLinesearchSolver,err,error,*999)
            CALL PETSc_SNESGetKSP(quasiNewtonLinesearchSolver%snes,iterativeSolver%ksp,err,error,*999)
          CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
            NULLIFY(quasiNewtonTrustregionSolver)
            CALL SolverNonlinearQuasiNewton_TrustregionSolverGet(quasiNewtonSolver,quasiNewtonTrustregionSolver,err,error,*999)
            CALL PETSc_SNESGetKSP(quasiNewtonTrustregionSolver%snes,iterativeSolver%ksp,err,error,*999)
          CASE DEFAULT
            localError="The Quasi-Newton solve type of "// &
              & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//"is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ELSE
        CALL PETSc_KSPCreate(groupCommunicator,iterativeSolver%ksp,err,error,*999)
      ENDIF
      !Set the iterative solver type
      SELECT CASE(iterativeSolver%iterativeSolverType)
      CASE(SOLVER_ITERATIVE_RICHARDSON)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPRICHARDSON,err,error,*999)
      CASE(SOLVER_ITERATIVE_CHEBYSHEV)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPCHEBYSHEV,err,error,*999)
      CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPCG,err,error,*999)
      CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPBICG,err,error,*999)
      CASE(SOLVER_ITERATIVE_GMRES)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPGMRES,err,error,*999)
        CALL PETSc_KSPGMRESSetRestart(iterativeSolver%ksp,iterativeSolver%gmresRestart,err,error,*999)
      CASE(SOLVER_ITERATIVE_BiCGSTAB)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPBCGS,err,error,*999)
      CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
        CALL PETSc_KSPSetType(iterativeSolver%ksp,PETSC_KSPCGS,err,error,*999)
      CASE DEFAULT
        localError="The iterative solver type of "// &
          & TRIM(NumberToVString(iterativeSolver%iterativeSolverType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Get the pre-conditioner
      CALL PETSc_KSPGetPC(iterativeSolver%ksp,iterativeSolver%pc,err,error,*999)
      !Set the pre-conditioner type
      SELECT CASE(iterativeSolver%iterativePreconditionerType)
      CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCNONE,err,error,*999)
      CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCJACOBI,err,error,*999)
      CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCBJACOBI,err,error,*999)
      CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCSOR,err,error,*999)
      CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCICC,err,error,*999)
      CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCILU,err,error,*999)
      CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
        CALL PETSc_PCSetType(iterativeSolver%pc,PETSC_PCASM,err,error,*999)
      CASE DEFAULT
        localError="The iterative preconditioner type of "// &
          & TRIM(NumberToVString(iterativeSolver%iterativePreconditionerType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the tolerances for the KSP solver
      CALL PETSc_KSPSetTolerances(iterativeSolver%ksp,iterativeSolver%relativeTolerance, &
        & iterativeSolver%absoluteTolerance,iterativeSolver%divergenceTolerance, &
        & iterativeSolver%maximumNumberOfIterations,err,error,*999)
      !Set any further KSP options from the command line options
      CALL PETSc_KSPSetFromOptions(iterativeSolver%ksp,err,error,*999)
      !Set the solver matrix to be the KSP matrix
      IF(solverMatrices%numberOfMatrices/=1) THEN
        localError="The number of solver matrices of "// &
          & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))// &
          & " is invalid. There should only be one solver matrix for a linear iterative solver."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
      NULLIFY(distributedMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
      NULLIFY(petscMatrix)
      CALL DistributedMatrix_PETScMatrixGet(distributedMatrix,petscMatrix,err,error,*999)
      CALL PETSc_KSPSetOperators(iterativeSolver%ksp,petscMatrix%matrix,petscMatrix%matrix,err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "// &
        & TRIM(NumberToVString(iterativeSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverLinearIterative_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverLinearIterative_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverLinearIterative_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum divergence tolerance for an iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeDivergenceToleranceSet
  SUBROUTINE Solver_LinearIterativeDivergenceToleranceSet(solver,divergenceTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set 
    REAL(DP), INTENT(IN) :: divergenceTolerance !<The divergence tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeDivergenceToleranceSet",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(divergenceTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified divergence tolerance of "// &
        & TRIM(NumberToVString(divergenceTolerance,"*",err,error))// &
        & " is invalid. The divergence tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    iterativeSolver%divergenceTolerance=divergenceTolerance
       
    EXITS("Solver_LinearIterativeDivergenceToleranceSet")
    RETURN
999 ERRORS("Solver_LinearIterativeDivergenceToleranceSet",err,error)
    EXITS("Solver_LinearIterativeDivergenceToleranceSet")
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeDivergenceToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Finalise an iterative linear solver for a linear solver and deallocate all memory.
  SUBROUTINE SolverLinear_IterativeFinalise(iterativeSolver,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer the linear iterative solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver

    ENTERS("SolverLinear_IterativeFinalise",err,error,*999)

    IF(ASSOCIATED(iterativeSolver)) THEN
      linearSolver=>iterativeSolver%linearSolver
      IF(ASSOCIATED(linearSolver)) THEN
        IF(.NOT.linearSolver%linkedNewtonPetSCSolver) THEN
          CALL PETSc_PCFinalise(iterativeSolver%pc,err,error,*999)
          CALL PETSc_KSPFinalise(iterativeSolver%ksp,err,error,*999)
        ENDIF
      ENDIF
      DEALLOCATE(iterativeSolver)
    ENDIF
        
    EXITS("SolverLinear_IterativeFinalise")
    RETURN
999 ERRORSEXITS("SolverLinear_IterativeFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_IterativeFinalise

  !
  !================================================================================================================================
  !

  !>Sets/changes the GMRES restart value for a GMRES iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeGMRESRestartSet
  SUBROUTINE Solver_LinearIterativeGMRESRestartSet(solver,gmresRestart,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the GMRES restart value
    INTEGER(INTG), INTENT(IN) :: gmresRestart !<The GMRES restart value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeGMRESRestartSet",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(iterativeSolver%iterativeSolverType/=SOLVER_ITERATIVE_GMRES) &
      & CALL FlagError("The linear iterative solver is not a GMRES linear iterative solver.",err,error,*999)
    IF(gmresRestart<1) THEN    
      localError="The specified GMRES restart value of "//TRIM(NumberToVString(gmresRestart,"*",err,error))// &
        & " is invalid. The GMRES restart value must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    iterativeSolver%gmresRestart=gmresRestart
    
    EXITS("Solver_LinearIterativeGMRESRestartSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeGMRESRestartSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeGMRESRestartSet
        
  !
  !================================================================================================================================
  !

  !>Initialise an iterative linear solver for a linear solver
  SUBROUTINE SolverLinear_IterativeInitialise(linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer the linear solver to initialise the iterative solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("SolverLinear_IterativeInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*998)
    IF(ASSOCIATED(linearSolver%iterativeSolver)) &
      & CALL FlagError("Iterative solver is already associated for this linear solver.",err,error,*998)
      
    !Allocate and initialise an iterative solver
    ALLOCATE(linearSolver%iterativeSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate linear solver iterative solver.",err,error,*999)
    linearSolver%iterativeSolver%linearSolver=>linearSolver
    linearSolver%iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    linearSolver%iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_GMRES
    linearSolver%iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_JACOBI_PRECONDITIONER
    linearSolver%iterativeSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    linearSolver%iterativeSolver%maximumNumberOfIterations=100000
    linearSolver%iterativeSolver%relativeTolerance=1.0E-05_DP
    linearSolver%iterativeSolver%absoluteTolerance=1.0E-10_DP
    linearSolver%iterativeSolver%divergenceTolerance=1.0E5_DP
    linearSolver%iterativeSolver%gmresRestart=30
    CALL PETSc_PCInitialise(linearSolver%iterativeSolver%pc,err,error,*999)
    CALL PETSc_KSPInitialise(linearSolver%iterativeSolver%ksp,err,error,*999)
        
    EXITS("SolverLinear_IterativeInitialise")
    RETURN
999 CALL SolverLinear_IterativeFinalise(linearSolver%iterativeSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverLinear_IterativeInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_IterativeInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for an iterative linear solver.
  SUBROUTINE SolverLinearIterative_LibraryTypeSet(iterativeSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer the iterative linear solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the iterative linear solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverLinearIterative_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Iterative linear solver is not associated.",err,error,*999)
    
    SELECT CASE(iterativeSolver%iterativeSolverType)
    CASE(SOLVER_ITERATIVE_RICHARDSON)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Richardson iterative linear solver."
      END SELECT
    CASE(SOLVER_ITERATIVE_CHEBYSHEV)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Chebychev iterative linear solver."
      END SELECT
    CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Conjugate gradient iterative linear solver."
      END SELECT
    CASE(SOLVER_ITERATIVE_GMRES)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a GMRES iterative linear solver."
      END SELECT
    CASE(SOLVER_ITERATIVE_BiCGSTAB)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a BiCGSTAB iterative linear solver."
      END SELECT
    CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        iterativeSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The specified solver library type of "// &
          & TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Conjugate gradient squared iterative linear solver."
      END SELECT
    CASE DEFAULT
      localError="The iterative linear solver type of "// &
        & TRIM(NumberToVString(iterativeSolver%iterativeSolverType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("SolverLinearIterative_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverLinearIterative_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearIterative_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of iterations for an iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeMaximumIterationsSet
  SUBROUTINE Solver_LinearIterativeMaximumIterationsSet(solver,maximumIterations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the maximum number of iterations
    INTEGER(INTG), INTENT(IN) :: maximumIterations !<The maximum number of iterations
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeMaximumIterationsSet",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(maximumIterations<=0) THEN
      localError="The specified maximum iterations of "//TRIM(NumberToVString(maximumIterations,"*",err,error))// &
        & " is invalid. The maximum number of iterations must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    iterativeSolver%maximumNumberOfIterations=maximumIterations
    
    EXITS("Solver_LinearIterativeMaximumIterationsSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeMaximumIterationsSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeMaximumIterationsSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of preconditioner for an iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativePreconditionerTypeSet
  SUBROUTINE Solver_LinearIterativePreconditionerTypeSet(solver,iterativePreconditionerType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: iterativePreconditionerType !<The type of iterative preconditioner to set \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativePreconditionerTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(iterativePreconditionerType/=iterativeSolver%iterativePreconditionerType) THEN
      !Intialise the new preconditioner type
      SELECT CASE(iterativeSolver%solverLibrary)
      CASE(SOLVER_PETSC_LIBRARY)
        SELECT CASE(iterativePreconditionerType)
        CASE(SOLVER_ITERATIVE_NO_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_NO_PRECONDITIONER
        CASE(SOLVER_ITERATIVE_JACOBI_PRECONDITIONER)
          CALL FlagError("Iterative Jacobi preconditioning is not implemented for a PETSc library.",err,error,*999)
        CASE(SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER
        CASE(SOLVER_ITERATIVE_SOR_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_SOR_PRECONDITIONER
        CASE(SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER 
        CASE(SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
        CASE(SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER)
          iterativeSolver%iterativePreconditionerType=SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER
        CASE DEFAULT
          localError="The iterative preconditioner type of "// &
            & TRIM(NumberToVString(iterativePreconditionerType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(iterativeSolver%solverLibrary,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_LinearIterativePreconditionerTypeSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativePreconditionerTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativePreconditionerTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the relative tolerance for an iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeRelativeToleranceSet
  SUBROUTINE Solver_LinearIterativeRelativeToleranceSet(solver,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set 
    REAL(DP), INTENT(IN) :: relativeTolerance !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeRelativeToleranceSet",err,error,*999)

    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(relativeTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified relative tolerance of "//TRIM(NumberToVString(relativeTolerance,"*",err,error))// &
        & " is invalid. The relative tolerance must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    iterativeSolver%relativeTolerance=relativeTolerance
    
    EXITS("Solver_LinearIterativeRelativeToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeRelativeToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeRelativeToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution initialise type for an iterative linear solver
  SUBROUTINE Solver_LinearIterativeSolutionInitialiseTypeSet(solver,solutionInitialiseType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set 
    INTEGER(INTG), INTENT(IN) :: solutionInitialiseType !<The solution initialise type to set \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeSolutionInitialiseTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    SELECT CASE(solutionInitialiseType)
    CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
      iterativeSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_ZERO
    CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
      iterativeSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
      iterativeSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_NO_CHANGE
    CASE DEFAULT
      localError="The specified solution initialise type of "// &
        & TRIM(NumberToVString(solutionInitialiseType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_LinearIterativeSolutionInitialiseTypeSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeSolutionInitialiseTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeSolutionInitialiseTypeSet
        
  !
  !================================================================================================================================
  !

  !>Solves a linear iterative linear solver
  SUBROUTINE SolverLinearIterative_Solve(iterativeSolver,err,error,*)

    !Argument variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver !<A pointer the linear iterative solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: convergedReason,globalRow,localRow,localRowIdx,numberOfIterations,numberOfMatrices,numberOfRows,storageType
    REAL(DP) :: residualNorm,solverValue,matrixValue
    REAL(DP), POINTER :: rhsData(:)
    TYPE(DistributedMatrixType), POINTER :: distributedMatrix
    TYPE(DistributedVectorType), POINTER :: rhsVector,solverVector
    TYPE(DistributedVectorPETScType), POINTER :: rhsPETScVector,solverPETScVector
    TYPE(DomainMappingType), POINTER :: rowDOFsMapping
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverLinearIterative_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(iterativeSolver)) CALL FlagError("Linear iterative solver is not associated.",err,error,*999)
    NULLIFY(linearSolver)
    CALL SolverLinearIterative_LinearSolverGet(iterativeSolver,linearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
    IF(numberOfMatrices/=1) THEN
      localError="The number of solver matrices of "//TRIM(NumberToVString(numberOfMatrices,"*",err,error))// &
        & " is invalid. There should only be one solver matrix for a linear direct solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(solverMatrix)
    CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
    NULLIFY(distributedMatrix)
    CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,distributedMatrix,err,error,*999)
    NULLIFY(rhsVector)
    CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,rhsVector,err,error,*999)
    NULLIFY(solverVector)
    CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
    CALL DistributedMatrix_StorageTypeGet(distributedMatrix,storageType,err,error,*999)
    IF(storageType==DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE) THEN
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfRowsGet(solverMapping,numberOfRows,err,error,*999)
      NULLIFY(rowDOFsMapping)
      CALL SolverMapping_RowDOFSMappingGet(solverMapping,rowDOFSMapping,err,error,*999)
      CALL DistributedVector_DataGet(rhsVector,rhsData,err,error,*999)
      DO localRowIdx=1,numberOfRows
        CALL DomainMapping_LocalToGlobalGet(rowDOFsMapping,localRowIdx,globalRow,err,error,*999)
        CALL DistributedMatrix_ValuesGet(distributedMatrix,localRowIdx,globalRow,matrixValue,err,error,*999)
        IF(ABS(matrixValue)<=ZERO_TOLERANCE) THEN
          localError="The linear solver matrix has a zero pivot at local row "// &
            & TRIM(NumberToVString(localRowIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        solverValue=rhsData(localRowIdx)/matrixValue
        CALL DistributedVector_ValuesSet(solverVector,localRowIdx,solverValue,err,error,*999)
      ENDDO !localRowIdx
      CALL DistributedVector_DataRestore(rhsVector,rhsData,err,error,*999)
    ELSE
      SELECT CASE(iterativeSolver%solverLibrary)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        NULLIFY(rhsPETScVector)
        CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
        NULLIFY(solverPETScVector)
        CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
        SELECT CASE(iterativeSolver%solutionInitialiseType)
        CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
          !Zero the solution vector
          CALL DistributedVector_AllValuesSet(solverVector,0.0_DP,err,error,*999)
          !Tell PETSc that the solution vector is zero
          CALL PETSc_KSPSetInitialGuessNonZero(iterativeSolver%ksp,.FALSE.,err,error,*999)
        CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
          !Make sure the solver vector contains the current dependent field values
          CALL Solver_SolutionUpdate(solver,err,error,*999)
          !Tell PETSc that the solution vector is nonzero
          CALL PETSc_KSPSetInitialGuessNonZero(iterativeSolver%ksp,.TRUE.,err,error,*999)
        CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
          !Do nothing
        CASE DEFAULT
          localError="The linear iterative solver solution initialise type of "// &
            & TRIM(NumberToVString(iterativeSolver%solutionInitialiseType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT                              
       !Set the tolerances for the KSP solver in case user has updated them
        CALL PETSc_KSPSetTolerances(iterativeSolver%ksp,iterativeSolver%relativeTolerance, &
          & iterativeSolver%absoluteTolerance,iterativeSolver%divergenceTolerance, &
          & iterativeSolver%maximumNumberOfIterations,err,error,*999)
        !Solver the linear system                             
        CALL PETSc_KSPSolve(iterativeSolver%ksp,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999)
        !Check for convergence
        CALL PETSc_KSPGetConvergedReason(iterativeSolver%ksp,convergedReason,err,error,*999)
        SELECT CASE(convergedReason)
        CASE(PETSC_KSP_DIVERGED_NULL)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged null.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_ITS)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged its.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_DTOL)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged dtol.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_BREAKDOWN)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged breakdown.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_BREAKDOWN_BICG)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged breakdown BiCG.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_NONSYMMETRIC)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged nonsymmetric.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_INDEFINITE_PC)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged indefinite PC.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_NANORINF)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged NaN or Inf.",err,error,*999)
        CASE(PETSC_KSP_DIVERGED_INDEFINITE_MAT)
          CALL FlagWarning("Linear iterative solver did not converge. PETSc diverged indefinite mat.",err,error,*999)
        END SELECT
        IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
          !Output solution characteristics
          CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear iterative solver parameters:",err,error,*999)
          CALL PETSc_KSPGetIterationNumber(iterativeSolver%ksp,numberOfIterations,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final number of iterations = ",numberOfIterations,err,error,*999)
          CALL PETSc_KSPGetResidualNorm(iterativeSolver%ksp,residualNorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final residual norm = ",residualNorm,err,error,*999)
          SELECT CASE(convergedReason)
          CASE(PETSC_KSP_CONVERGED_RTOL)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged RTol",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_ATOL)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged ATol",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_ITS)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged its",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_CG_NEG_CURVE)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged CG neg curve",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_CG_CONSTRAINED)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged CG constrained",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_STEP_LENGTH)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged step length",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_HAPPY_BREAKDOWN)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged happy breakdown",err,error,*999)
          CASE(PETSC_KSP_CONVERGED_ITERATING)
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged iterating",err,error,*999)
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The solver library type of "// &
          & TRIM(NumberToVString(iterativeSolver%solverLibrary,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("SolverLinearIterative_Solve")
    RETURN
999 ERRORSEXITS("SolverLinearIterative_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinearIterative_Solve
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of iterative linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearIterativeTypeSet
  SUBROUTINE Solver_LinearIterativeTypeSet(solver,iterativeSolverType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the iterative linear solver type
    INTEGER(INTG), INTENT(IN) :: iterativeSolverType !<The type of iterative linear solver to set \see SolverRoutines_IterativeLinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_LinearIterativeTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsLinear(solver,err,error,*999)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*999)
    CALL SolverLinear_AssertIsIterative(linearSolver,err,error,*999)
    NULLIFY(iterativeSolver)
    CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
    IF(iterativeSolverType/=iterativeSolver%iterativeSolverType) THEN
      !Intialise the new solver type
      SELECT CASE(iterativeSolver%solverLibrary)
      CASE(SOLVER_PETSC_LIBRARY)
        SELECT CASE(iterativeSolverType)
        CASE(SOLVER_ITERATIVE_RICHARDSON)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_RICHARDSON
        CASE(SOLVER_ITERATIVE_CHEBYSHEV)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_CHEBYSHEV
        CASE(SOLVER_ITERATIVE_CONJUGATE_GRADIENT)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_CONJUGATE_GRADIENT
        CASE(SOLVER_ITERATIVE_BICONJUGATE_GRADIENT)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_BICONJUGATE_GRADIENT
        CASE(SOLVER_ITERATIVE_GMRES)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_GMRES
        CASE(SOLVER_ITERATIVE_BiCGSTAB)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_BiCGSTAB
        CASE(SOLVER_ITERATIVE_CONJGRAD_SQUARED)
          iterativeSolver%iterativeSolverType=SOLVER_ITERATIVE_CONJGRAD_SQUARED
        CASE DEFAULT
          localError="The iterative solver type of "//TRIM(NumberToVString(iterativeSolverType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The solver library type of "// &
          & TRIM(NumberToVString(solver%linearSolver%iterativeSolver%solverLibrary,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_LinearIterativeTypeSet")
    RETURN
999 ERRORSEXITS("Solver_LinearIterativeTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearIterativeTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a linear solver.
  SUBROUTINE SolverLinear_LibraryTypeSet(linearSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer the linear solver to get the library type for.
     INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the linear solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinear_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)
    
    SELECT CASE(linearSolver%linearSolveType)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      NULLIFY(directSolver)
      CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
      CALL SolverLinearDirect_LibraryTypeSet(directSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      NULLIFY(iterativeSolver)
      CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
      CALL SolverLinearIterative_LibraryTypeSet(iterativeSolver,solverLibraryType,err,error,*999)
    CASE DEFAULT
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("SolverLinear_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverLinear_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Solve a linear solver 
  SUBROUTINE SolverLinear_Solve(linearSolver,err,error,*)

    !Argument variables
    TYPE(LinearSolverType), POINTER :: linearSolver !<A pointer to the linear solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfMatrices,solverMatrixIdx
    TYPE(DistributedVectorType), POINTER :: distributedSolverVector
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverLinear_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(linearSolver)) CALL FlagError("Linear solver is not associated.",err,error,*999)

    NULLIFY(solver)
    CALL SolverLinear_SolverGet(linearSolver,solver,err,error,*999)
    IF(.NOT.ASSOCIATED(solver%linkingSolver)) THEN
      !Assemble the solver matrices
!!TODO: Work out what to assemble
      CALL Solver_StaticAssemble(solver,SOLVER_MATRICES_LINEAR_ONLY,err,error,*999)
    ENDIF
    SELECT CASE(linearSolver%linearSolveType)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      CALL SolverLinearDirect_Solve(linearSolver%directSolver,err,error,*999)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      CALL SolverLinearIterative_Solve(linearSolver%iterativeSolver,err,error,*999)
    CASE DEFAULT
      localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver solution vectors:",err,error,*999)
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfMatrices,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of solution vectors = ",numberOfMatrices,err,error,*999)
      DO solverMatrixIdx=1,numberOfMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solution vector for solver matrix : ",solverMatrixIdx,err,error,*999)
        NULLIFY(distributedSolverVector)
        CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,distributedSolverVector,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,distributedSolverVector,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
        
    IF(.NOT.ASSOCIATED(solver%linkingSolver)) THEN
      !Update depenent field with solution
      CALL Solver_VariablesFieldUpdate(solver,err,error,*999)
    ENDIF
        
    EXITS("SolverLinear_Solve")
    RETURN
999 ERRORSEXITS("SolverLinear_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverLinear_Solve
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of linear solver. \see OpenCMISS::Iron::cmfe_Solver_LinearTypeSet
  SUBROUTINE Solver_LinearTypeSet(solver,linearSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the linear solver type
    INTEGER(INTG), INTENT(IN) :: linearSolveType !<The type of linear solver to set \see SolverRoutines_LinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solver_LinearTypeSet",err,error,*998)

    CALL Solver_AssertNotFinished(solver,err,error,*998)
    CALL Solver_AssertIsLinear(solver,err,error,*998)
    NULLIFY(linearSolver)
    CALL Solver_LinearSolverGet(solver,linearSolver,err,error,*998)
    IF(linearSolveType/=linearSolver%linearSolveType) THEN
      !Intialise the new solver type
      SELECT CASE(linearSolveType)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SolverLinear_DirectInitialise(solver%linearSolver,err,error,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SolverLinear_IterativeInitialise(solver%linearSolver,err,error,*999)
      CASE DEFAULT
        localError="The linear solver type of "//TRIM(NumberToVString(linearSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old solver type
      SELECT CASE(solver%linearSolver%linearSolveType)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        CALL SolverLinear_DirectFinalise(solver%linearSolver%directSolver,err,error,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        CALL SolverLinear_IterativeFinalise(solver%linearSolver%iterativeSolver,err,error,*999)
      CASE DEFAULT
        localError="The linear solver type of "// &
          & TRIM(NumberToVString(solver%linearSolver%linearSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      solver%linearSolver%linearSolveType=linearSolveType
    ENDIF
    
    EXITS("Solver_LinearTypeSet")
    RETURN
999 SELECT CASE(linearSolveType)
    CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
      CALL SolverLinear_DirectFinalise(linearSolver%directSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
      CALL SolverLinear_IterativeFinalise(linearSolver%iterativeSolver,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("Solver_LinearTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_LinearTypeSet
        
  !
  !================================================================================================================================
  !

  !>Assembles the solver matrices and rhs from the dynamic equations.
  SUBROUTINE Solver_DynamicAssemble(solver,selectionType,err,error,*)
    
    !Argument variable
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: selectionType !<The type of matrix selection \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnNumber,currentIteration,dampingMatrixNumber,dependentVariableType,dirichletIdx,dynamicVariableType, &
      & equationsColumnNumber,equationsMatrixIdx,equationsMatrixIdx2,equationsMatrixNumber,equationsRowNumber, &
      & equationsSetIdx,inputIteration,interfaceColumnNumber,interfaceConditionIdx,interfaceConditionMethod,interfaceMatrixIdx, &
      & interfaceVariableType,jacobianMatrixIdx,lhsBoundaryCondition,lhsBoundaryFinish,lhsGlobalDOF,lhsVariableDOF, &
      & linearMatrixIdx,linearVariableIdx, &
      & linearVariableType,massMatrixNumber,numberOfDirichletConditions,numberOfDynamicMatrices,numberOfEquationsMatrices, &
      & numberOfEquationsSets, &
      & numberOfInterfaceConditions,numberOfInterfaceMatrices,numberOfJacobianMatrices,numberOfLinearMatrices, &
      & numberOfLinearVariables,numberOfResiduals,numberOfResidualVariables,numberOfRows, &
      & numberOfSolverMatrices,numberOfSources,outputIteration,residualIdx,residualVariableDOF,residualVariableIdx, &
      & rhsBoundaryCondition,rhsGlobalDOF, &
      & rhsVariableDOF,rhsVariableType,rowCondition,solverRowIdx,solverRowNumber,solverMatrixIdx,sourceIdx,stiffnessMatrixNumber, &
      & timeDependenceType,totalNumberOfRows,transposeTimeDependenceType,variableBoundaryCondition,variableDOF, &
      & variableGlobalDOF,variableIdx,variablePositionIdx,variableType
    INTEGER(INTG), POINTER :: equationsRowToLHSDOFMap(:),equationsRowToResidualDOFMap(:),equationsRowToRHSDOFMap(:), &
      & variableDOFToRowMap(:)
    REAL(SP) :: systemElapsed,systemTime1(1),systemTime2(1),userElapsed,userTime1(1),userTime2(1)
    REAL(DP) :: alphaValue,currentFunctionFactor,currentRHSValue,currentTime,dampingMatrixCoefficient,deltaT,dofValue, &
      & dynamicAccelerationFactor,dynamicDisplacementFactor,dynamicValue,dynamicVelocityFactor,equationsDampingCoefficient, &
      & equationsMassCoefficient,equationsStiffnessCoefficient,firstUpdateFactor,jacobianMatrixCoefficient,linearCoefficient, &
      & linearValue,linearValueSum,massMatrixCoefficient,matrixCoefficient,matrixCoefficients(2)=[0.0_DP,0.0_DP],nonlinearValue, &
      & previousFunctionFactor,previous2FunctionFactor,previous3FunctionFactor,previousRHSValue,previous2RHSValue, &
      & previous3RHSValue,residualCoefficient,residualValue,rhsValue,rowCouplingCoefficient,secondUpdateFactor, &
      & solverRHSValue,sourceValue,stiffnessMatrixCoefficient,sourceCoefficient,startTime,stopTime,timeIncrement, &
      & transposeMatrixCoefficient,vectorCoefficient
    REAL(DP), POINTER :: matrixCheckData(:),currentValuesVector(:),previousValuesVector(:),previousVelocityVector(:), &
      & previousAccelerationVector(:),previousResidualParameters(:),previous2ResidualParameters(:), &
      & previous3ResidualParameters(:),rhsIntegratedParameters(:),rhsParameters(:),solverRHSCheckData(:),solverResidualCheckData(:)
    LOGICAL :: hasIntegratedValues,hasTranspose,includeResidual,interfaceMatrixDynamic,rhsLinearMatrix,rhsResidual, &
      & updateResidual,updateRHS,updateSolverMatrix
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: dependentBoundaryConditions,lhsBoundaryConditionsVariable, &
      & rhsBoundaryConditionsVariable
    TYPE(BoundaryConditionsRowVariableType), POINTER :: lhsBoundaryConditionsRowVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DistributedMatrixType), POINTER :: dampingDistributedMatrix,equationsDistributedMatrix,interfaceDistributedMatrix, &
      & jacobianDistributedMatrix,linearDistributedMatrix,massDistributedMatrix,previousSolverDistributedMatrix, &
      & solverDistributedMatrix,stiffnessDistributedMatrix,transposeDistributedMatrix,transposeInterfaceDistributedMatrix
    TYPE(DistributedVectorType), POINTER :: currentDistributedVector,currentResidualVector,currentRHSVector,currentSourceVector, &
      & dependentDistributedVector,distributedResidualVector,distributedSourceVector,dynamicTempVector,equationsRHSVector, &
      & incrementalVector,interfaceRHSDistributedVector,interfaceTempVector,lagrangeDistributedVector,linearTempVector, &
      & nonlinearTempVector,previousDistributedVector,previous2DistributedVector,previous3DistributedVector, &
      & predictedMeanAccelerationVector, &
      & predictedMeanDisplacementVector,predictedMeanVelocityVector,previousResidualVector,previous2ResidualVector, &
      & previous3ResidualVector,previousRHSVector,previous2RHSVector,previous3RHSVector, &
      & previousSourceVector,previous2SourceVector,previous3SourceVector,residualDistributedVector,solverRHSVector, &
      & solverResidualVector,sourcesTempVector,transposeInterfaceTempVector
    TYPE(DomainMappingType), POINTER :: lhsDomainMapping,residualDomainMapping,rhsDomainMapping,variableDomainMapping
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,linearMatrix,massMatrix,stiffnessMatrix,equationsMatrix
    TYPE(EquationsMatrixToSolverMatricesMapType), POINTER :: equationsMatrixToSolverMatricesMap
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap,linearMatrixToSolverMatrixMap
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,dynamicVariable,interfaceVariable,lagrangeVariable,lhsVariable, &
      & linearVariable,residualVariable,rhsVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: interfaceRHSMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(InterfaceMatrixToVarMapType), POINTER :: interfaceMatrixToVarMap
    TYPE(InterfaceRHSType), POINTER :: interfaceRHSVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:),equationsRowToSolverRowsMap(:), &
      & interfaceColToSolverColsMap(:),interfaceColToSolverRowsMap(:),interfaceRowToSolverColsMap(:), &
      & interfaceRowToSolverRowsMap(:),jacobianColToSolverColsMap(:)
    TYPE(SolverType), POINTER :: linkingSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VarToEquationsMatricesMapType), POINTER :: dynamicVarToEquationsMatricesMap,linearVarToEquationsMatricesMap
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Solver_DynamicAssemble",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
       
    !Determine which dynamic solver needs to be used
    NULLIFY(dynamicSolver)
    IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
      CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    ELSE IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
      NULLIFY(linkingSolver)
      CALL Solver_LinkingSolverGet(solver,linkingSolver,err,error,*999)
      CALL Solver_DynamicSolverGet(linkingSolver,dynamicSolver,err,error,*999)
    ELSE
      localError="The solver solve type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
        & " is invalid for for dynamic assembly."
      CALL FlagError(localError,err,error,*999)
    ENDIF

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime, &
      & currentIteration,outputIteration,inputIteration,err,error,*999)
    
    deltaT=dynamicSolver%timeIncrement
    SELECT CASE(dynamicSolver%degree)
    CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
      stiffnessMatrixCoefficient=1.0_DP*dynamicSolver%theta(1)*deltaT
      dampingMatrixCoefficient=1.0_DP
      massMatrixCoefficient=0.0_DP
      jacobianMatrixCoefficient=stiffnessMatrixCoefficient
      dynamicDisplacementFactor=deltaT
      currentFunctionFactor=dynamicSolver%theta(1)
      previousFunctionFactor=1.0_DP-dynamicSolver%theta(1)
    CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
      stiffnessMatrixCoefficient=dynamicSolver%theta(2)*deltaT*deltaT/2.0_DP
      dampingMatrixCoefficient=dynamicSolver%theta(1)*deltaT
      massMatrixCoefficient=1.0_DP
      jacobianMatrixCoefficient=(dynamicSolver%theta(2)+dynamicSolver%theta(1))*deltaT*deltaT/4.0_DP
      firstUpdateFactor=deltaT
      dynamicDisplacementFactor=deltaT*deltaT/2.0_DP
      dynamicVelocityFactor=deltaT
      currentFunctionFactor=(dynamicSolver%theta(2)+dynamicSolver%theta(1))/2.0_DP
      previousFunctionFactor=1.0_DP-dynamicSolver%theta(2)
      previous2FunctionFactor=(dynamicSolver%theta(2)-dynamicSolver%theta(1))/2.0_DP
    CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
      stiffnessMatrixCoefficient=dynamicSolver%theta(3)*deltaT*deltaT*deltaT/6.0_DP
      dampingMatrixCoefficient=dynamicSolver%theta(2)*deltaT*deltaT/2.0_DP
      massMatrixCoefficient=dynamicSolver%theta(1)*deltaT
      jacobianMatrixCoefficient=(dynamicSolver%theta(3)+3.0_DP*dynamicSolver%theta(2)+2.0_DP*dynamicSolver%theta(1))* &
        & deltaT*deltaT*deltaT/36.0_DP
      firstUpdateFactor=deltaT
      secondUpdateFactor=deltaT*deltaT/2.0_DP
      dynamicDisplacementFactor=deltaT*deltaT*deltaT/6.0_DP
      dynamicVelocityFactor=deltaT
      dynamicAccelerationFactor=deltaT*deltaT/2.0_DP
      currentFunctionFactor=(dynamicSolver%theta(3)+3.0_DP*dynamicSolver%theta(2)+2.0_DP*dynamicSolver%theta(1))/6.0_DP
      previousFunctionFactor=1.0_DP-(dynamicSolver%theta(3)+2.0_DP*dynamicSolver%theta(2)-dynamicSolver%theta(1))/2.0_DP
      previous2FunctionFactor=(dynamicSolver%theta(3)+dynamicSolver%theta(2)-2.0_DP*dynamicSolver%theta(1))/2.0_DP
      previous3FunctionFactor=(dynamicSolver%theta(1)-dynamicSolver%theta(3))/6.0_DP            
    CASE DEFAULT
      localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)

    NULLIFY(previousSolverDistributedMatrix)
    
    !Assemble the solver matrices    
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
      IF(dynamicSolver%solverInitialised.OR.(.NOT.dynamicSolver%solverInitialised.AND. &
        & ((dynamicSolver%order==SOLVER_DYNAMIC_FIRST_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE).OR. &
        & (dynamicSolver%order==SOLVER_DYNAMIC_SECOND_ORDER.AND.dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE)))) THEN
        !Assemble solver matrices
        IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
          CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
          CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
        ENDIF
        !Just deal with one solver matrix for now. 
        solverMatrixIdx=1
        CALL SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*999)
        IF(numberOfSolverMatrices/=solverMatrixIdx) CALL FlagError("Invalid number of solver matrices.",err,error,*999)

        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL SolverMatrix_UpdateMatrixGet(solverMatrix,updateSolverMatrix,err,error,*999)
        IF(updateSolverMatrix) THEN
          NULLIFY(solverDistributedMatrix)
          CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)    
          !Initialise matrix to zero
          CALL DistributedMatrix_AllValuesSet(solverDistributedMatrix,0.0_DP,err,error,*999)
          !Get the check data
          NULLIFY(matrixCheckData)
          CALL DistributedMatrix_DataGet(solverDistributedMatrix,matrixCheckData,err,error,*999)
          !Loop over the equations sets
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            NULLIFY(vectorMapping)
            CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
            NULLIFY(vectorMatrices)
            CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
            NULLIFY(equationsSetToSolverMatricesMap)
            CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
              & err,error,*999)
            NULLIFY(equationsRowToSolverRowsMap)
            CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
              & equationsRowToSolverRowsMap,err,error,*999)
            NULLIFY(equationsMatricesToSolverMatrixMap)
            CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
              & equationsMatricesToSolverMatrixMap,err,error,*999)
            NULLIFY(dynamicMapping)
            CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
            IF(ASSOCIATED(dynamicMapping)) THEN
              NULLIFY(dynamicMatrices)
              CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
              IF(dynamicSolver%solverInitialised) THEN
                CALL EquationsMappingDynamic_StiffnessMatrixNumberGet(dynamicMapping,stiffnessMatrixNumber,err,error,*999)
                IF(stiffnessMatrixNumber/=0) THEN
                  NULLIFY(equationsMatrixToSolverMatricesMap)
                  CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                    & stiffnessMatrixNumber,equationsMatrixToSolverMatricesMap,err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatrixMap)
                  CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(equationsColToSolverColsMap)
                  CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                    & equationsColToSolverColsMap,err,error,*999)
                  NULLIFY(stiffnessMatrix)
                  CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,stiffnessMatrixNumber,stiffnessMatrix, &
                    & err,error,*999)
                  NULLIFY(stiffnessDistributedMatrix)
                  CALL EquationsMatrix_DistributedMatrixGet(stiffnessMatrix,stiffnessDistributedMatrix,err,error,*999)
                  CALL EquationsMatrix_MatrixCoefficientGet(stiffnessMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=matrixCoefficient*stiffnessMatrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,stiffnessDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ENDIF
                CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
                IF(dampingMatrixNumber/=0) THEN
                  NULLIFY(equationsMatrixToSolverMatricesMap)
                  CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                    & dampingMatrixNumber,equationsMatrixToSolverMatricesMap,err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatrixMap)
                  CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(equationsColToSolverColsMap)
                  CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                    & equationsColToSolverColsMap,err,error,*999)
                  NULLIFY(dampingMatrix)
                  CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix, &
                    & err,error,*999)
                  NULLIFY(dampingDistributedMatrix)
                  CALL EquationsMatrix_DistributedMatrixGet(dampingMatrix,dampingDistributedMatrix,err,error,*999)
                  CALL EquationsMatrix_MatrixCoefficientGet(dampingMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=matrixCoefficient*dampingMatrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,dampingDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ENDIF
                CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
                IF(massMatrixNumber/=0) THEN
                  NULLIFY(equationsMatrixToSolverMatricesMap)
                  CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                    & massMatrixNumber,equationsMatrixToSolverMatricesMap,err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatrixMap)
                  CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(equationsColToSolverColsMap)
                  CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                    & equationsColToSolverColsMap,err,error,*999)
                  NULLIFY(massMatrix)
                  CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
                  NULLIFY(massDistributedMatrix)
                  CALL EquationsMatrix_DistributedMatrixGet(massMatrix,massDistributedMatrix,err,error,*999)
                  CALL EquationsMatrix_MatrixCoefficientGet(massMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=matrixCoefficient*massMatrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,massDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ENDIF
              ELSE
                IF(dynamicSolver%order==SOLVER_DYNAMIC_SECOND_ORDER.AND. &
                  & dynamicSolver%degree==SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
                  IF(massMatrixNumber==0) CALL FlagError("Can not perform initial solve with no mass matrix.",err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatricesMap)
                  CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                    & massMatrixNumber,equationsMatrixToSolverMatricesMap,err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatrixMap)
                  CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(equationsColToSolverColsMap)
                  CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                    & equationsColToSolverColsMap,err,error,*999)
                  NULLIFY(massMatrix)
                  CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
                  NULLIFY(massDistributedMatrix)
                  CALL EquationsMatrix_DistributedMatrixGet(massMatrix,massDistributedMatrix,err,error,*999)
                  CALL EquationsMatrix_MatrixCoefficientGet(massMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=-1.0_DP*matrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,massDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ELSE
                  CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
                  IF(dampingMatrixNumber==0) CALL FlagError("Can not perform initial solve with no damping matrix.",err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatricesMap)
                  CALL SolverMappingESToSMSMap_EquationsMatrixToSolverMatricesMapGet(equationsSetToSolverMatricesMap, &
                    & dampingMatrixNumber,equationsMatrixToSolverMatricesMap,err,error,*999)
                  NULLIFY(equationsMatrixToSolverMatrixMap)
                  CALL SolverMappingEMToSMSMap_EquationsMatrixToSolverMatrixMapGet(equationsMatrixToSolverMatricesMap, &
                    & solverMatrixIdx,equationsMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(equationsColToSolverColsMap)
                  CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(equationsMatrixToSolverMatrixMap, &
                    & equationsColToSolverColsMap,err,error,*999)
                  NULLIFY(dampingMatrix)
                  CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix,err,error,*999)
                  NULLIFY(dampingDistributedMatrix)
                  CALL EquationsMatrix_DistributedMatrixGet(dampingMatrix,dampingDistributedMatrix,err,error,*999)
                  CALL EquationsMatrix_MatrixCoefficientGet(dampingMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=-1.0_DP*matrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,dampingDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ENDIF
              ENDIF !dynamic solver initialised
            ENDIF
            NULLIFY(nonlinearMapping)
            CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
            IF(ASSOCIATED(nonlinearMapping)) THEN
              IF(selectionType==SOLVER_MATRICES_ALL.OR. &
                & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
                & selectionType==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
                !Now set the values from the equations Jacobian
                CALL SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap, &
                  & numberOfJacobianMatrices,err,error,*999)
                DO jacobianMatrixIdx=1,numberOfJacobianMatrices
                  NULLIFY(jacobianMatrixToSolverMatrixMap)
                  CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                    & jacobianMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
                  NULLIFY(jacobianColToSolverColsMap)
                  CALL SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap, &
                    & jacobianColToSolverColsMap,err,error,*999)
                  NULLIFY(jacobianMatrix)
                  CALL SolverMappingJMToSMMap_JacobianMatrixGet(jacobianMatrixToSolverMatrixMap,jacobianMatrix,err,error,*999)
                  NULLIFY(jacobianDistributedMatrix)
                  CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,jacobianDistributedMatrix,err,error,*999)
                  CALL JacobianMatrix_MatrixCoefficientGet(jacobianMatrix,matrixCoefficient,err,error,*999)
                  matrixCoefficient=matrixCoefficient*jacobianMatrixCoefficient
                  CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & equationsRowToSolverRowsMap,jacobianColToSolverColsMap,matrixCoefficient,jacobianDistributedMatrix, &
                    & .FALSE.,err,error,*999)
                ENDDO !jacobianMatrixIdx
              ENDIF
            ENDIF
            NULLIFY(linearMapping)
            CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
            IF(ASSOCIATED(linearMapping)) THEN
              CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap, &
                & numberOfLinearMatrices,err,error,*999)
              DO linearMatrixIdx=1,numberOfLinearMatrices
                NULLIFY(linearMatrixToSolverMatrixMap)
                CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                  & linearMatrixIdx,linearMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(equationsColToSolverColsMap)
                CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(linearMatrixToSolverMatrixMap, &
                  & equationsColToSolverColsMap,err,error,*999)
                NULLIFY(linearMatrix)
                CALL SolverMappingEMToSMMap_EquationsMatrixGet(linearMatrixToSolverMatrixMap,linearMatrix,err,error,*999)
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)
                matrixCoefficient=matrixCoefficient*stiffnessMatrixCoefficient
                CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,linearDistributedMatrix, &
                  & .FALSE.,err,error,*999)
              ENDDO !equationsMatrixIdx
            ENDIF
          ENDDO !equationsSetIdx
          !Loop over any interface conditions
          DO interfaceConditionIdx=1,numberOfInterfaceConditions
            NULLIFY(interfaceConditionToSolverMatricesMap)
            CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
              & interfaceConditionToSolverMatricesMap,err,error,*999)
            NULLIFY(interfaceColToSolverRowsMap)
            CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceColToSolverRowsMap,err,error,*999)
            NULLIFY(interfaceMatricesToSolverMatrixMap)
            CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
              & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
            NULLIFY(interfaceColToSolverColsMap)
            CALL SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet(interfaceMatricesToSolverMatrixMap, &
              & interfaceColToSolverColsMap,err,error,*999)
            !Loop over the interface matrices
            CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
              & numberOfInterfaceMatrices,err,error,*999)
            DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
              NULLIFY(interfaceMatrixToSolverMatricesMap)
              CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
                & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
              NULLIFY(interfaceRowToSolverRowsMap)
              CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
                & interfaceRowToSolverRowsMap,err,error,*999)
              NULLIFY(interfaceMatrixToSolverMatrixMap)
              CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(interfaceMatrix)
              CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*999)
              NULLIFY(interfaceDistributedMatrix)
              CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
              CALL InterfaceMatrix_MatrixCoefficientGet(interfaceMatrix,matrixCoefficients(1),err,error,*999)
              CALL InterfaceMatrix_TimeDependenceTypeGet(interfaceMatrix,timeDependenceType,err,error,*999)
              SELECT CASE(timeDependenceType)
              CASE(INTERFACE_MATRIX_STATIC)
                matrixCoefficients(1)=matrixCoefficients(1)*stiffnessMatrixCoefficient
              CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                matrixCoefficients(1)=matrixCoefficients(1)*dampingMatrixCoefficient
              CASE(INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
                matrixCoefficients(1)=matrixCoefficients(1)*massMatrixCoefficient
              CASE DEFAULT
                localError="The interface matrix time dependence type of "// &
                  & TRIM(NumberToVString(timeDependenceType,"*",err,error))//" is invalid for interface matrix index "// &
                  & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" of interface condition index "// &
                  & TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END SELECT              
              CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                & interfaceRowToSolverRowsMap,interfaceColToSolverColsMap,matrixCoefficients(1),interfaceDistributedMatrix, &
                & .FALSE.,err,error,*999)
              CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
              IF(hasTranspose) THEN
                NULLIFY(interfaceRowToSolverColsMap)
                CALL SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap, &
                  & interfaceRowToSolverColsMap,err,error,*999)
                NULLIFY(transposeDistributedMatrix)
                CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*999)
                CALL InterfaceMatrix_TransposeMatrixCoefficientGet(interfaceMatrix,matrixCoefficients(2),err,error,*999)
                CALL InterfaceMatrix_TransposeTimeDependenceTypeGet(interfaceMatrix,transposeTimeDependenceType,err,error,*999)
                SELECT CASE(transposeTimeDependenceType)
                CASE(INTERFACE_MATRIX_STATIC)
                  matrixCoefficients(2)=matrixCoefficients(2)*stiffnessMatrixCoefficient
                CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                  matrixCoefficients(2)=matrixCoefficients(2)*dampingMatrixCoefficient
                CASE(INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
                  matrixCoefficients(2)=matrixCoefficients(2)*massMatrixCoefficient
                CASE DEFAULT
                  localError="The interface matrix transpose time dependence type of "// &
                    & TRIM(NumberToVString(transposeTimeDependenceType,"*",err,error))// &
                    & " is invalid for interface matrix index "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
                    & " of interface condition index "//TRIM(NumberToVString(interfaceConditionIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & interfaceColToSolverRowsMap,interfaceRowToSolverColsMap,matrixCoefficients(2),transposeDistributedMatrix, &
                  & .FALSE.,err,error,*999)
              ENDIF
            ENDDO !interfaceMatrixIdx
          ENDDO !interfaceConditionIdx
          !Update the solver matrix values
          CALL DistributedMatrix_UpdateStart(solverDistributedMatrix,err,error,*999)          
          IF(ASSOCIATED(previousSolverDistributedMatrix)) &
            & CALL DistributedMatrix_UpdateFinish(previousSolverDistributedMatrix,err,error,*999)
          previousSolverDistributedMatrix=>solverDistributedMatrix
          IF(solver%solveType==SOLVER_DYNAMIC_TYPE) THEN
            IF(dynamicSolver%solverInitialised) solverMatrix%updateMatrix=.FALSE.
          ELSE IF(solver%solveType==SOLVER_NONLINEAR_TYPE) THEN 
            IF(dynamicSolver%solverInitialised) solverMatrix%updateMatrix=.TRUE.
          ELSE
            CALL FlagError("Dynamic solver solve type is not associated.",err,error,*999)
          END IF
          !Restore check data
          CALL DistributedMatrix_DataRestore(solverDistributedMatrix,matrixCheckData,err,error,*999)
        ENDIF !Update matrix
        
        IF(ASSOCIATED(previousSolverDistributedMatrix)) &
          & CALL DistributedMatrix_UpdateFinish(previousSolverDistributedMatrix,err,error,*999)
        
        IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
          CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
          CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
          userElapsed=userTime2(1)-userTime1(1)
          systemElapsed=systemTime2(1)-systemTime1(1)
          IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
            & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
          CALL Profiling_TimingsOutput(1,"Solver matrices assembly",userElapsed,systemElapsed,err,error,*999)
        ENDIF
      ENDIF
    ENDIF !Calculate solver matrix

    NULLIFY(solverRHSVector)
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_ONLY) THEN
      !Assemble rhs vector
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF

      CALL SolverMatrices_UpdateRHSGet(solverMatrices,updateRHS,err,error,*999)
      IF(updateRHS) THEN
        CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,solverRHSVector,err,error,*999)
        !Initialise the RHS to zero
        CALL DistributedVector_AllValuesSet(solverRHSVector,0.0_DP,err,error,*999)          
        !Get the solver RHS check data                  
        NULLIFY(solverRHSCheckData)
        CALL DistributedVector_DataGet(solverRHSVector,solverRHSCheckData,err,error,*999)
        !Get the boundary conditions
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)

        !Loop over the equations sets
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(equationsSetToSolverMatricesMap)
          CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
            & err,error,*999)
          NULLIFY(equationsRowToSolverRowsMap)
          CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
            & equationsRowToSolverRowsMap,err,error,*999)

          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(lhsMapping)
          CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
          NULLIFY(lhsVariable)
          CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
          NULLIFY(equationsRowToLHSDOFMap)
          CALL EquationsMappingLHS_EquationsRowTOLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*999)
          NULLIFY(lhsDomainMapping)
          CALL FieldVariable_DomainMappingGet(lhsVariable,lhsDomainMapping,err,error,*999)
          CALL DomainMapping_BoundaryFinishGet(lhsDomainMapping,lhsBoundaryFinish,err,error,*999)
          CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*999)
          CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*999)
          NULLIFY(currentValuesVector)
          CALL FieldVariable_ParameterSetDataGet(lhsVariable,FIELD_VALUES_SET_TYPE,currentValuesVector,err,error,*999)
          NULLIFY(previousValuesVector)
          CALL FieldVariable_ParameterSetDataGet(lhsVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,previousValuesVector,err,error,*999)
          NULLIFY(previousVelocityVector)
          NULLIFY(previousAccelerationVector)
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
            CALL FieldVariable_ParameterSetDataGet(lhsVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,previousVelocityVector, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetDataGet(lhsVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,previousAccelerationVector, &
                & err,error,*999)
            ENDIF
          ENDIF

          interfaceMatrixDynamic=.FALSE.
          NULLIFY(dynamicTempVector)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
          IF(ASSOCIATED(dynamicMapping)) THEN
            !Calculate the dynamic contributions
            NULLIFY(dynamicVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
            CALL FieldVariable_VariableTypeGet(dynamicVariable,dynamicVariableType,err,error,*999)
            NULLIFY(dynamicVarToEquationsMatricesMap)
            CALL EquationsMappingDynamic_VariableToEquationsMatricesMapGet(dynamicMapping,dynamicVarToEquationsMatricesMap, &
              & err,error,*999)
            CALL EquationsMappingDynamic_NumberOfDynamicMatricesGet(dynamicMapping,numberOfDynamicMatrices,err,error,*999)
            NULLIFY(dynamicMatrices)
            CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
            CALL EquationsMatricesDynamic_TempDistributedVectorGet(dynamicMatrices,dynamicTempVector,err,error,*999)
            !Initialise the dynamic temporary vector to zero
            CALL DistributedVector_AllValuesSet(dynamicTempVector,0.0_DP,err,error,*999)
            CALL EquationsMappingDynamic_StiffnessMatrixNumberGet(dynamicMapping,stiffnessMatrixNumber,err,error,*999)
            IF(stiffnessMatrixNumber/=0) THEN
              NULLIFY(stiffnessMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,stiffnessMatrixNumber,stiffnessMatrix, &
                & err,error,*999)
              NULLIFY(stiffnessDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(stiffnessMatrix,stiffnessDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(stiffnessMatrix,equationsStiffnessCoefficient,err,error,*999)
              NULLIFY(predictedMeanDisplacementVector)
              CALL FieldVariable_ParameterSetVectorGet(dynamicVariable,FIELD_MEAN_PREDICTED_DISPLACEMENT_SET_TYPE, &
                & predictedMeanDisplacementVector,err,error,*999)              
              CALL DistributedMatrix_MatrixByVectorAdd(stiffnessDistributedMatrix,.FALSE.,predictedMeanDisplacementVector, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,equationsStiffnessCoefficient,dynamicTempVector,err,error,*999)         
            ENDIF
            CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
            IF(dampingMatrixNumber/=0) THEN
              NULLIFY(dampingMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix,err,error,*999)
              NULLIFY(dampingDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(dampingMatrix,dampingDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(dampingMatrix,equationsDampingCoefficient,err,error,*999)
              IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
                NULLIFY(predictedMeanVelocityVector)
                CALL FieldVariable_ParameterSetVectorGet(dynamicVariable,FIELD_MEAN_PREDICTED_VELOCITY_SET_TYPE, &
                  & predictedMeanVelocityVector,err,error,*999)
                CALL DistributedMatrix_MatrixByVectorAdd(dampingDistributedMatrix,.FALSE.,predictedMeanVelocityVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,equationsDampingCoefficient,dynamicTempVector,err,error,*999)
              ENDIF
            ENDIF
            CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
            IF(massMatrixNumber/=0.AND.dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE) THEN
              NULLIFY(massMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
              NULLIFY(massDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(massMatrix,massDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(massMatrix,equationsMassCoefficient,err,error,*999)
              IF(dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                NULLIFY(predictedMeanAccelerationVector)
                CALL FieldVariable_ParameterSetVectorGet(dynamicVariable,FIELD_MEAN_PREDICTED_ACCELERATION_SET_TYPE, &
                  & predictedMeanAccelerationVector,err,error,*999)
                CALL DistributedMatrix_MatrixByVectorAdd(massDistributedMatrix,.FALSE.,predictedMeanAccelerationVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,equationsMassCoefficient,dynamicTempVector,err,error,*999)
              ENDIF
            ENDIF
            !Work out if there are any interface matrices mapped to the dynamic variable.
            CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
            interfConditionLoop: DO interfaceConditionIdx=1,numberOfInterfaceConditions
              NULLIFY(interfaceCondition)
              CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
              NULLIFY(interfaceConditionToSolverMatricesMap)
              CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
                & interfaceConditionToSolverMatricesMap,err,error,*999)
              NULLIFY(interfaceMatricesToSolverMatrixMap)
              CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
                & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
              NULLIFY(interfaceEquations)
              CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
              NULLIFY(interfaceMapping)
              CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
              NULLIFY(interfaceMatrices)
              CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
              CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
                & numberOfInterfaceMatrices,err,error,*999)
              interfaceMatrixLoop: DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
                NULLIFY(interfaceMatrixToSolverMatrixMap)
                CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                  & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(interfaceMatrix)
                CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*999)
                CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
                IF(hasTranspose) THEN
                  NULLIFY(interfaceMatrixToVarMap)
                  CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx, &
                    & interfaceMatrixToVarMap,err,error,*999)
                  NULLIFY(dependentVariable)
                  CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,dependentVariable,err,error,*999)
                  IF(ASSOCIATED(dependentVariable,dynamicVariable)) THEN
                    interfaceMatrixDynamic=.TRUE.
                    EXIT interfConditionLoop
                  ENDIF
                ENDIF
              ENDDO interfaceMatrixLoop !interfaceMatrixIdx
            ENDDO interfConditionLoop !interfaceConditionIdx
          ENDIF !dynamic mapping

          NULLIFY(nonlinearMapping)
          CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
          IF(ASSOCIATED(nonlinearMapping)) THEN
            NULLIFY(nonlinearMatrices)
            CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
            NULLIFY(nonlinearTempVector)
            CALL EquationsMatricesNonlinear_TempDistributedVectorGet(nonlinearMatrices,nonlinearTempVector,err,error,*999)
            CALL DistributedVector_AllValuesSet(nonlinearTempVector,0.0_DP,err,error,*999)
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              rhsResidual=.TRUE.
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(solverMatrixToEquationsMap)
                CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
                  & err,error,*999)
                NULLIFY(solverMappingVariables)
                CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
                DO residualVariableIdx=1,numberOfResidualVariables
                  CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
                  NULLIFY(solverMappingVariable)
                  CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,residualVariable,solverMappingVariable, &
                    & err,error,*999)
                  IF(ASSOCIATED(solverMappingVariable)) rhsResidual=.FALSE.
                ENDDO !residualVariableIdx
              ENDDO !solverMatrixIdx
              IF(rhsResidual.OR.solver%solveType==SOLVER_NONLINEAR_TYPE) THEN
                CALL EquationsMappingResidual_VectorCoefficientGet(residualMapping,residualCoefficient,err,error,*999)
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                IF(rhsResidual) THEN
                  !The residual does not have any variables involved in the solver matrices so take it to the RHS
                  NULLIFY(currentResidualVector)
                  CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                    & currentResidualVector,err,error,*999)
                  vectorCoefficient=residualCoefficient*currentFunctionFactor
                  CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & vectorCoefficient,currentResidualVector,err,error,*999)
                ENDIF
                NULLIFY(previousResidualVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
                  & previousResidualVector,err,error,*999)
                vectorCoefficient=residualCoefficient*previousFunctionFactor
                CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & vectorCoefficient,previousResidualVector,err,error,*999)
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                  NULLIFY(previous2ResidualVector)
                  CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                    & previous2ResidualVector,err,error,*999)
                  vectorCoefficient=residualCoefficient*previous2FunctionFactor
                  CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                    & vectorCoefficient,previous2ResidualVector,err,error,*999)
                  IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                    NULLIFY(previous3ResidualVector)
                    CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                      & previous3ResidualVector,err,error,*999)
                    vectorCoefficient=residualCoefficient*previous3FunctionFactor
                    CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                      & vectorCoefficient,previous3ResidualVector,err,error,*999)
                  ENDIF
                ENDIF
              ENDIF !rhs or nonlinear residual
            ENDDO !residualIdx
          ENDIF !nonlinear mapping

          !Calculate the contributions from any linear matrices
          NULLIFY(linearTempVector)
          NULLIFY(linearMapping)
          NULLIFY(linearMatrices)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
            CALL EquationsMatricesLinear_TempDistributedVectorGet(linearMatrices,linearTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(linearTempVector,0.0_DP,err,error,*999)                  
            CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,linearVariable,err,error,*999)
              rhsLinearMatrix=.TRUE.
              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(solverMatrixToEquationsMap)
                CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
                  & err,error,*999)
                NULLIFY(solverMappingVariables)
                CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,linearVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) rhsLinearMatrix=.FALSE.
              ENDDO !solverMatrixIdx
              IF(rhsLinearMatrix) THEN
                !Linear matrix variable is not on the LHS so take the matrix times vector over to the RHS.
                NULLIFY(linearMatrix)
                CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,linearMatrix,err,error,*999)
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,linearCoefficient,err,error,*999)
!!TODO: Could store these mat-vec products so that we do not have to do the matrix vector products every time step.
                NULLIFY(currentDistributedVector)
                CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_VALUES_SET_TYPE,currentDistributedVector, &
                  & err,error,*999)
                matrixCoefficient=linearCoefficient*currentFunctionFactor
                CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,currentDistributedVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
                NULLIFY(previousDistributedVector)
                CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                  & previousDistributedVector,err,error,*999)
                matrixCoefficient=linearCoefficient*previousFunctionFactor
                CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,previousDistributedVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                  NULLIFY(previous2DistributedVector)
                  CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_PREVIOUS2_VALUES_SET_TYPE, &
                    & previous2DistributedVector,err,error,*999)
                  matrixCoefficient=linearCoefficient*previous2FunctionFactor
                  CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,previous2DistributedVector, &
                    & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
                  IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                    NULLIFY(previous3DistributedVector)
                    CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_PREVIOUS3_VALUES_SET_TYPE, &
                      & previous3DistributedVector,err,error,*999)
                    matrixCoefficient=linearCoefficient*previous3FunctionFactor
                    CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,previous3DistributedVector, &
                      & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
                  ENDIF
                ENDIF
              ENDIF !rhs linear matrix
            ENDDO !equationsMatrixIdx
          ENDIF !linear mapping

          !Add in any source vectors
          NULLIFY(sourcesTempVector)
          NULLIFY(sourcesMapping)
          NULLIFY(sourceVectors)
          CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
          IF(ASSOCIATED(sourcesMapping)) THEN
            CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
            CALL EquationsMatricesSources_TempDistributedVectorGet(sourceVectors,sourcesTempVector,err,error,*999)
            CALL DistributedVector_AllValuesSet(sourcesTempVector,0.0_DP,err,error,*999)
            CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
            DO sourceIdx=1,numberOfSources
              NULLIFY(sourceMapping)
              CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
              CALL EquationsMappingSource_VectorCoefficientGet(sourceMapping,sourceCoefficient,err,error,*999)
              NULLIFY(sourceVector)
              CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
              NULLIFY(distributedSourceVector)
              CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                & distributedSourceVector,err,error,*999)
              sourceCoefficient=sourceCoefficient*currentFunctionFactor
              CALL DistributedVector_VectorAdd(sourcesTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,sourceCoefficient, &
                & currentSourceVector,err,error,*999)
              NULLIFY(previousSourceVector)
              CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
                & previousSourceVector,err,error,*999)
              sourceCoefficient=sourceCoefficient*previousFunctionFactor
              CALL DistributedVector_VectorAdd(sourcesTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,sourceCoefficient, &
                & previousSourceVector,err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                NULLIFY(previous2SourceVector)
                CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                  & previous2SourceVector,err,error,*999)
                sourceCoefficient=sourceCoefficient*previous2FunctionFactor
                CALL DistributedVector_VectorAdd(sourcesTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,sourceCoefficient, &
                  & previous2SourceVector,err,error,*999)     
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  NULLIFY(previous3SourceVector)
                  CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                    & previous3SourceVector,err,error,*999)
                  sourceCoefficient=sourceCoefficient*previous3FunctionFactor
                  CALL DistributedVector_VectorAdd(sourcesTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,sourceCoefficient, &
                    & previous3SourceVector,err,error,*999) 
                ENDIF
              ENDIF
            ENDDO !sourceIdx
          ENDIF !source mapping

          NULLIFY(rhsVariable)
          NULLIFY(rhsMapping)
          NULLIFY(rhsVector)
          NULLIFY(currentRHSVector)
          NULLIFY(previousRHSVector)
          NULLIFY(previous2RHSVector)
          NULLIFY(previous3RHSVector)
          NULLIFY(rhsBoundaryConditionsVariable)
          NULLIFY(rhsIntegratedParameters)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL FieldVariable_ParameterSetCreated(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,hasIntegratedValues, &
              & err,error,*999)
            CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
            CALL EquationsMatricesRHS_DistributedVectorGet(rhsVector,EQUATIONS_MATRICES_CURRENT_VECTOR,currentRHSVector, &
              & err,error,*999)
            CALL EquationsMatricesRHS_DistributedVectorGet(rhsVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR,previousRHSVector, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
              CALL EquationsMatricesRHS_DistributedVectorGet(rhsVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR,previous2RHSVector, &
                & err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                CALL EquationsMatricesRHS_DistributedVectorGet(rhsVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR,previous3RHSVector, &
                  & err,error,*999)
              ENDIF
            ENDIF
            CALL FieldVariable_ParameterSetDataGet(rhsVariable,FIELD_VALUES_SET_TYPE,rhsParameters,err,error,*999) 
            IF(hasIntegratedValues) THEN
              !Update RHS field by integrating any point Neumann conditions
              CALL BoundaryConditions_VariableGet(boundaryConditions,rhsVariable,rhsBoundaryConditionsVariable,err,error,*999)
!!TODO: NEED TO FIX INTEGRATED FLUX MATRIX CALCULATION. THE MATRIX CURRENTLY USES THE RHS FOR THE ROWS BUT IT SHOULD USE THE
!!      LHS. IT SHOULD BE POSSIBLE TO USE A DIFFERENT BASIS FOR THE FLUX INTERPOLATION THAN FOR THE DEPENDENT INTERPOLATION.
              CALL BoundaryConditionsVariable_NeumannIntegrate(rhsBoundaryConditionsVariable,err,error,*999)
              CALL FieldVariable_ParameterSetDataGet(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhsIntegratedParameters, &
                & err,error,*999) 
            ENDIF
            NULLIFY(equationsRowToRHSDOFMap)
            CALL EquationsMappingRHS_EquationsRowToRHSDOFMapGet(rhsMapping,equationsRowToRHSDOFMap,err,error,*999)
          ENDIF

          NULLIFY(lhsBoundaryConditionsRowVariable)
          CALL BoundaryConditions_RowVariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsRowVariable,err,error,*999)
          NULLIFY(lhsBoundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsVariable,err,error,*999)
          CALL BoundaryCOnditionsVariable_NumberOfDirichletConditionsGet(lhsBoundaryConditionsVariable, &
            & numberOfDirichletConditions,err,error,*999)
          NULLIFY(dirichletBoundaryConditions)
          IF(numberOfDirichletConditions>0) CALL BoundaryConditionsVariable_DirichletConditionsExists( &
            & lhsBoundaryConditionsVariable,dirichletBoundaryConditions,err,error,*999)
          NULLIFY(neumannBoundaryConditions)
          CALL BoundaryConditionsVariable_NeumannConditionsExists(lhsBoundaryConditionsVariable, &
            & neumannBoundaryConditions,err,error,*999)

          !Loop over the rows in the equations set
          DO equationsRowNumber=1,totalNumberOfRows
            CALL BoundaryConditionsRowVariable_RowConditionTypeGet(lhsBoundaryConditionsRowVariable,equationsRowNumber, &
              & rowCondition,err,error,*999)

            !Get the dynamic contribution to the RHS values               
            IF(ASSOCIATED(dynamicTempVector)) THEN
              CALL DistributedVector_ValuesGet(dynamicTempVector,equationsRowNumber,dynamicValue,err,error,*999)
            ELSE
              dynamicValue=0.0_DP
            ENDIF
            !Get the nonlinear contribution to the RHS values               
            IF(ASSOCIATED(nonlinearTempVector)) THEN
              CALL DistributedVector_ValuesGet(nonlinearTempVector,equationsRowNumber,nonlinearValue,err,error,*999)
            ELSE
              nonlinearValue=0.0_DP
            ENDIF
            !Get the linear contribution to the RHS values               
            IF(ASSOCIATED(linearTempVector)) THEN
              CALL DistributedVector_ValuesGet(linearTempVector,equationsRowNumber,linearValue,err,error,*999)
            ELSE
              linearValue=0.0_DP
            ENDIF
            !Get the source contribution to the RHS values               
            IF(ASSOCIATED(sourcesTempVector)) THEN
              CALL DistributedVector_ValuesGet(sourcesTempVector,equationsRowNumber,sourceValue,err,error,*999)                
            ELSE
              sourceValue=0.0_DP
            ENDIF

            solverRHSValue=dynamicValue+nonlinearValue+linearValue+sourceValue

            IF(ASSOCIATED(rhsMapping)) THEN
              rhsVariableDOF=equationsRowToRHSDOFMap(equationsRowNumber)
              IF(hasIntegratedValues) THEN
                !Add any Neumann integrated values, b = f + N q
                CALL DistributedVector_ValuesAdd(currentRHSVector,equationsRowNumber,rhsIntegratedParameters(rhsVariableDOF), &
                  & err,error,*999)
              ELSE
                CALL DistributedVector_ValuesSet(currentRHSVector,equationsRowNumber,rhsParameters(rhsVariableDOF), &
                  & err,error,*999)
              ENDIF
              CALL DistributedVector_ValuesGet(currentRHSVector,equationsRowNumber,currentRHSValue,err,error,*999)
              rhsValue=currentRHSValue*currentFunctionFactor
              CALL DistributedVector_ValuesGet(previousRHSVector,equationsRowNumber,previousRHSValue,err,error,*999)
              rhsValue=rhsValue+previousRHSValue*previousFunctionFactor
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                CALL DistributedVector_ValuesGet(previous2RHSVector,equationsRowNumber,previous2RHSValue,err,error,*999)
                rhsValue=rhsValue+previous2RHSValue*previous2FunctionFactor
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  CALL DistributedVector_ValuesGet(previous3RHSVector,equationsRowNumber,previous3RHSValue,err,error,*999)
                  rhsValue=rhsValue+previous3RHSValue*previous3FunctionFactor
                ENDIF
              ENDIF
            ELSE
              rhsValue=0.0_DP
            ENDIF

            !Get the dynamic contribution to the the RHS values
            lhsVariableDOF=equationsRowToLHSDOFMap(equationsRowNumber)
            CALL DomainMapping_LocalToGlobalGet(lhsDomainMapping,lhsVariableDOF,lhsGlobalDOF,err,error,*999)
            CALL BoundaryConditionsVariable_DOFTypeGet(lhsBoundaryConditionsVariable,lhsGlobalDOF,lhsBoundaryCondition, &
              & err,error,*999)

            !Apply boundary conditions
            SELECT CASE(lhsBoundaryCondition)
            CASE(BOUNDARY_CONDITION_FREE_ROW)
              solverRHSValue=solverRHSValue+rhsValue
              !Loop over the solver rows associated with this equations set row
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,solverRHSValue,err,error,*999)
            CASE(BOUNDARY_CONDITION_DIRICHLET_ROW)
!!WHY ARE WE ADDING TO THE RHS ROW IF WE ARE GOING TO ELIMINATE IT?
              !Get the equations RHS values
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,rhsValue,err,error,*999)

              !Note: Find the incremental alpha values for Dirichlet DOFs
              !      For cases in which the boundary condition does not change with time then alpha will be
              !      implicitly zero and the contribution to the RHS will have been included above with the dynamic value.
              !      However, for cases in which the boundary condition is changing we need to compute the corresponding
              !      alpha value and include it on the RHS. This alpha value also needs to be set for nolinear problems as
              !      it will not have come from the nonlinear solver as A.alpha is included in the residual below. 
              IF(lhsBoundaryCondition==BOUNDARY_CONDITION_DOF_FIXED) THEN
                SELECT CASE(dynamicSolver%degree)
                CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                  alphaValue=(currentValuesVector(lhsVariableDOF)-previousValuesVector(lhsVariableDOF))/ &
                    & dynamicDisplacementFactor
                CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                  alphaValue=(currentValuesVector(lhsVariableDOF)- &
                    & previousValuesVector(lhsVariableDOF)- &
                    & dynamicVelocityFactor*previousVelocityVector(lhsVariableDOF))/ &
                    & dynamicDisplacementFactor
                CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                  alphaValue=(currentValuesVector(lhsVariableDOF)- &
                    & previousValuesVector(lhsVariableDOF)- &
                    & dynamicVelocityFactor*previousVelocityVector(lhsVariableDOF) - &
                    & dynamicAccelerationFactor*previousAccelerationVector(lhsVariableDOF))/ &
                    & dynamicDisplacementFactor
                CASE DEFAULT
                  localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT

                !Update the incremented values for non-linear problems
                IF(ASSOCIATED(nonlinearMapping)) THEN
                  !Only update if the DOF is not a ghost DOF
                  IF(lhsVariableDOF<=lhsBoundaryFinish) CALL FieldVariable_ParameterSetUpdateLocalDOF(lhsVariable, &
                    & FIELD_INCREMENTAL_VALUES_SET_TYPE,lhsVariableDOF,alphaValue,err,error,*999)
                ENDIF

                IF(ABS(alphaValue)>=ZERO_TOLERANCE) THEN
                  !Loop over the column entries of the matrices corresponding to the fixed DOF column and take over to the RHS.
                  !First for any equations matrices, then interface matrices.
                  IF(ASSOCIATED(dynamicMapping)) THEN
                    IF(ASSOCIATED(dynamicVariable,lhsVariable)) THEN
                      CALL EquationsMappingVectorVToEMSMap_NumberOfEquationsMatricesGet(dynamicVarToEquationsMatricesMap, &
                        & numberOfEquationsMatrices,err,error,*999)
                      DO equationsMatrixIdx=1,numberOfEquationsMatrices
                        CALL EquationsMappingVectorVToEMSMap_EquationsMatrixNumberGet(dynamicVarToEquationsMatricesMap, &
                          & equationsMatrixIdx,equationsMatrixNumber,err,error,*999)
                        CALL EquationsMappingVectorVToEMSMap_EquationsMatrixColumnNumberGet(dynamicVarToEquationsMatricesMap, &
                          & equationsMatrixIdx,lhsVariableDOF,columnNumber,err,error,*999)
                        IF(equationsMatrixNumber==stiffnessMatrixNumber) &
                          & dofValue=alphaValue*stiffnessMatrixCoefficient*equationsStiffnessCoefficient
                        IF(equationsMatrixNumber==dynamicMapping%dampingMatrixNumber) &
                          & dofValue=alphaValue*dampingMatrixCoefficient*equationsDampingCoefficient                    
                        IF(equationsMatrixNumber==dynamicMapping%massMatrixNumber) &
                          & dofValue=alphaValue*massMatrixCoefficient*equationsMassCoefficient
                        NULLIFY(equationsMatrix)
                        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,equationsMatrixNumber,equationsMatrix, &
                          & err,error,*999)
                        NULLIFY(equationsDistributedMatrix)
                        CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,equationsDistributedMatrix,err,error,*999)
!!WHAT IS THIS TRYING TO DO?
                        DO dirichletIdx=1,numberOfDirichletConditions
                          IF(dirichletBoundaryConditions%dirichletDOFIndices(dirichletIdx)==equationsColumnNumber) EXIT
                        ENDDO !dirichletIdx
                        CALL DistributedMatrix_MatrixColumnAdd(equationsDistributedMatrix,.FALSE.,equationsColumnNumber, &
                          & equationsRowToSolverRowsMap,-1.0_DP*dofValue,solverRHSVector,err,error,*999)
                      ENDDO !matrixIdx
                    ENDIF
                  ENDIF !dynamic mapping
                  IF(ASSOCIATED(linearMapping)) THEN
                    CALL EquationsMappingLinear_NumberOfLinearVariablesGet(linearMapping,numberOfLinearVariables,err,error,*999)
                    DO linearVariableIdx=1,numberOfLinearVariables
                      NULLIFY(linearVariable)
                      CALL EquationsMappingLinear_linearVariableGet(linearMapping,linearVariableIdx,linearVariable,err,error,*999)
                      IF(ASSOCIATED(linearVariable,lhsVariable)) THEN
                        NULLIFY(linearVarToEquationsMatricesMap)
                        CALL EquationsMappingLinear_VariableToEquationsMatricesMapGet(linearMapping,linearVariableIdx, &
                          & linearVarToEquationsMatricesMap,err,error,*999)
                        CALL EquationsMappingVectorVToEMSMap_NumberOfEquationsMatricesGet(linearVarToEquationsMatricesMap, &
                          & numberOfEquationsMatrices,err,error,*999)
                        DO equationsMatrixIdx=1,numberOfEquationsMatrices
                          CALL EquationsMappingVectorVToEMSMap_EquationsMatrixNumberGet(linearVarToEquationsMatricesMap, &
                            & equationsMatrixIdx,equationsMatrixNumber,err,error,*999)
                          CALL EquationsMappingVectorVToEMSMap_EquationsMatrixColumnNumberGet(linearVarToEquationsMatricesMap, &
                            & equationsMatrixIdx,lhsVariableDOF,columnNumber,err,error,*999)
                          NULLIFY(equationsMatrix)
                          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixNumber,equationsMatrix, &
                            & err,error,*999)
                          CALL EquationsMatrix_MatrixCoefficientGet(equationsMatrix,matrixCoefficient,err,error,*999)
                          dofValue=alphaValue*stiffnessMatrixCoefficient*equationsStiffnessCoefficient
                          NULLIFY(equationsDistributedMatrix)
                          CALL EquationsMatrix_DistributedMatrixGet(equationsMatrix,equationsDistributedMatrix,err,error,*999)
!!WHAT IS THIS TRYING TO DO?
                          DO dirichletIdx=1,numberOfDirichletConditions
                            IF(dirichletBoundaryConditions%dirichletDOFIndices(dirichletIdx)==equationsColumnNumber) EXIT
                          ENDDO !dirichletIdx
                          CALL DistributedMatrix_MatrixColumnAdd(equationsDistributedMatrix,.FALSE.,equationsColumnNumber, &
                            & equationsRowToSolverRowsMap,-1.0_DP*dofValue,solverRHSVector,err,error,*999)
                        ENDDO !matrixidx
                      ENDIF
                    ENDDO !linearVariableIdx
                  ENDIF !linear mapping
                  !Loop over any interface matrices which are mapped to the dependent variable with the fixed BC.
                  IF(interfaceMatrixDynamic) THEN
                    DO interfaceConditionIdx=1,numberOfInterfaceConditions
                      NULLIFY(interfaceCondition)
                      CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition, &
                        & err,error,*999)
                      NULLIFY(interfaceConditionToSolverMatricesMap)
                      CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
                        & interfaceConditionToSolverMatricesMap,err,error,*999)
                      NULLIFY(interfaceMatricesToSolverMatrixMap)
                      CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
                        & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
                      NULLIFY(interfaceEquations)
                      CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
                      NULLIFY(interfaceMapping)
                      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
                      NULLIFY(interfaceMatrices)
                      CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
                      NULLIFY(interfaceColToSolverRowsMap)
                      CALL SolverMapping_InterfaceColToSolverRowsMapGet(solverMapping,interfaceConditionIdx, &
                        & interfaceColToSolverRowsMap,err,error,*999)                     
                      CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
                        & numberOfInterfaceMatrices,err,error,*999)
                      DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
                        NULLIFY(interfaceMatrixToSolverMatrixMap)
                        CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                          & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
                        NULLIFY(interfaceMatrix)
                        CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix, &
                          & err,error,*999)
                        CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
                        IF(hasTranspose) THEN
                          NULLIFY(interfaceMatrixToVarMap)
                          CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx, &
                            & interfaceMatrixToVarMap,err,error,*999)
                          NULLIFY(dependentVariable)
                          CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,dependentVariable,err,error,*999)
                          IF(ASSOCIATED(dependentVariable,lhsVariable)) THEN
                            CALL InterfaceMatrix_TransposeMatrixCoefficientGet(interfaceMatrix,transposeMatrixCoefficient, &
                              & err,error,*999)
                            CALL InterfaceMatrix_TransposeTimeDependenceTypeGet(interfaceMatrix,transposeTimeDependenceType, &
                              & err,error,*999)
                            SELECT CASE(transposeTimeDependenceType)
                            CASE(INTERFACE_MATRIX_STATIC)
                              dofValue=transposeMatrixCoefficient*stiffnessMatrixCoefficient*alphaValue
                            CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                              dofValue=transposeMatrixCoefficient*dampingMatrixCoefficient*alphaValue
                            CASE(INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
                              dofValue=transposeMatrixCoefficient*massMatrixCoefficient*alphaValue
                            CASE DEFAULT
                              CALL FlagError("Not implemented.",err,error,*999)
                            END SELECT
                            NULLIFY(variableDOFToRowMap)
                            CALL InterfaceMappingIMToVMap_VariableDOFToRowMapGet(interfaceMatrixToVarMap, &
                              & variableDOFToRowMap,err,error,*999)
                            interfaceColumnNumber=variableDOFToRowMap(variableDof)
                            CALL DistributedMatrix_MatrixColumnAdd(interfaceMatrix%matrix,.TRUE.,interfaceColumnNumber, &
                              & interfaceColToSolverRowsMap,-1.0_DP*dofValue,solverRHSVector,err,error,*999)
                          ENDIF
                        ENDIF
                      ENDDO !interfaceMatrixIdx
                    ENDDO !interfaceConditionIdx
                  ENDIF !dynamic interface matrix
                ENDIF !alpha > 0
              ENDIF !LHS DOF fixed

            CASE(BOUNDARY_CONDITION_NEUMANN_ROW)
              !Set Neumann boundary conditions
              solverRHSValue=solverRHSValue+rhsValue
              !Loop over the solver rows associated with this equations set row
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,solverRHSValue,err,error,*999)
            CASE(BOUNDARY_CONDITION_ROBIN_ROW)
              !Set Robin boundary conditions
              CALL FlagError("Robin boundary conditions are not implemented.",err,error,*999)
            CASE(BOUNDARY_CONDITION_CAUCHY_ROW)
              !Set Cauchy boundary conditions
              CALL FlagError("Cauchy boundary conditions are not implemented.",err,error,*999)
            CASE(BOUNDARY_CONDITION_CONSTRAINED_ROW)
              !Set constrained row boundary conditions
              CALL FlagError("Constrained row boundary conditions are not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The RHS boundary condition of "// &
                & TRIM(NumberToVString(lhsBoundaryCondition,"*",err,error))// &
                & " for RHS variable dof number "// &
                & TRIM(NumberToVString(lhsVariableDOF,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsRowNumber
          IF(.NOT.dynamicSolver%solverInitialised) THEN
            !Copy current RHS i.e., RHS at time zero, to previous RHSs to initialise
            !Only do this for non Ghost rows
            CALL DistributedVector_VectorCopy(previousRHSVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
              & currentRHSVector,err,error,*999)
            IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
              CALL DistributedVector_VectorCopy(previous2RHSVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                & currentRHSVector,err,error,*999)
              IF(dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                CALL DistributedVector_VectorCopy(previous3RHSVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                  & currentRHSVector,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) &
            & CALL FieldVariable_ParameterSetUpdateStart(lhsVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(lhsVariable,FIELD_VALUES_SET_TYPE,currentValuesVector,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(lhsVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,previousValuesVector, &
            & err,error,*999)
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
            CALL FieldVariable_ParameterSetDataRestore(lhsVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,previousVelocityVector, &
              & err,error,*999)
            IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetDataRestore(lhsVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE, &
                & previousAccelerationVector,err,error,*999)
            ENDIF
          ENDIF

          IF(ASSOCIATED(rhsMapping)) THEN
            CALL FieldVariable_ParameterSetDataRestore(rhsVariable,FIELD_VALUES_SET_TYPE,rhsParameters,err,error,*999)
            IF(hasIntegratedValues) &
              & CALL FieldVariable_ParameterSetDataRestore(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhsIntegratedParameters, &
              & err,error,*999)
          ENDIF
          IF(ASSOCIATED(nonlinearMapping)) &
            & CALL FieldVariable_ParameterSetUpdateFinish(lhsVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
        ENDDO !equationsSetIdx
        
        !Add in any rows from any interface conditions
        DO interfaceConditionIdx=1,numberOfInterfaceConditions
          NULLIFY(interfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
          SELECT CASE(interfaceConditionMethod)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            NULLIFY(interfaceEquations)
            CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
            NULLIFY(interfaceMapping)
            CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
            NULLIFY(interfaceMatrices)
            CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
            NULLIFY(interfaceRHSMapping)
            CALL InterfaceMapping_RHSMappingGet(interfaceMapping,interfaceRHSMapping,err,error,*999)
            NULLIFY(interfaceRHSVector)
            CALL InterfaceMatrices_RHSVectorGet(interfaceMatrices,interfaceRHSVector,err,error,*999)
            NULLIFY(interfaceColToSolverRowsMap)
            CALL SolverMapping_InterfaceColToSolverRowsMapGet(solverMapping,interfaceConditionIdx, &
              & interfaceColToSolverRowsMap,err,error,*999)
            NULLIFY(interfaceRHSDistributedVector)
            CALL InterfaceMatricesRHS_DistributedVectorGet(interfaceRHSVector,interfaceRHSDistributedVector,err,error,*999)
            !Worry about BCs on the Lagrange variables later.
            CALL DistributedVector_VectorCoupleAdd(solverRHSVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
              & interfaceColToSolverRowsMap,1.0_DP,interfaceRHSDistributedVector,err,error,*999)
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interface condition method of "// &
              & TRIM(NumberToVString(interfaceCondition%METHOD,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !interfaceConditionIdx
        !        
        !Start the update the solver RHS vector values
        CALL DistributedVector_UpdateStart(solverRHSVector,err,error,*999)
        CALL DistributedVector_UpdateFinish(solverRHSVector,err,error,*999)
        !Restore the solver RHS check data                 
        CALL DistributedVector_DataRestore(solverRHSVector,solverRHSCheckData,err,error,*999)            
      ENDIF !Update RHS

      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
          & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
        CALL Profiling_TimingsOutput(1,"Solver RHS assembly",userElapsed,systemElapsed,err,error,*999)
      ENDIF

    ENDIF !Calculate solver RHS

    NULLIFY(solverResidualVector)
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN

      !Assemble residual vector
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF

      CALL SolverMatrices_UpdateResidualGet(solverMatrices,updateResidual,err,error,*999)
      IF(updateResidual) THEN
        NULLIFY(solverResidualVector)
        CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,solverResidualVector,err,error,*999)
        !Get the list of solver variables
        NULLIFY(solverMatrixToEquationsMap)
        CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,1,solverMatrixToEquationsMap,err,error,*999)
        NULLIFY(solverMappingVariables)
        CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
        !Initialise the residual to zero
        CALL DistributedVector_AllValuesSet(solverResidualVector,0.0_DP,err,error,*999)
        !Get the solver residual check data
        NULLIFY(solverResidualCheckData)
        CALL DistributedVector_DataGet(solverResidualVector,solverResidualCheckData,err,error,*999)
        !Loop over the equations sets
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(lhsMapping)
          CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
          NULLIFY(lhsVariable)
          CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
          NULLIFY(lhsDomainMapping)
          CALL FieldVariable_DomainMappingGet(lhsVariable,lhsDomainMapping,err,error,*999)

          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
          IF(ASSOCIATED(dynamicMapping)) THEN
            !Calculate the dynamic contributions
            NULLIFY(dynamicVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
            NULLIFY(dynamicMatrices)
            CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
            NULLIFY(dynamicTempVector)
            CALL EquationsMatricesDynamic_TempDistributedVectorGet(dynamicMatrices,dynamicTempVector,err,error,*999)
            !Initialise the dynamic temporary vector to zero
            CALL DistributedVector_AllValuesSet(dynamicTempVector,0.0_DP,err,error,*999)
            NULLIFY(incrementalVector)
            !Define the pointer to the incrementalVector
            CALL FieldVariable_ParameterSetVectorGet(dynamicVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,incrementalVector, &
              & err,error,*999)
            CALL EquationsMappingDynamic_StiffnessMatrixNumberGet(dynamicMapping,stiffnessMatrixNumber,err,error,*999)
            IF(stiffnessMatrixNumber/=0) THEN
              NULLIFY(stiffnessMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,stiffnessMatrixNumber,stiffnessMatrix, &
                & err,error,*999)
              NULLIFY(stiffnessDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(stiffnessMatrix,stiffnessDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(stiffnessMatrix,matrixCoefficient,err,error,*999)
              matrixCoefficient=matrixCoefficient*stiffnessMatrixCoefficient
              CALL DistributedMatrix_MatrixByVectorAdd(stiffnessDistributedMatrix,.FALSE.,incrementalVector, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,dynamicTempVector,err,error,*999)
            ENDIF
            CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
            IF(dampingMatrixNumber/=0.AND.dynamicSolver%degree>=SOLVER_DYNAMIC_FIRST_DEGREE) THEN
              NULLIFY(dampingMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix,err,error,*999)
              NULLIFY(dampingDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(dampingMatrix,dampingDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(dampingMatrix,matrixCoefficient,err,error,*999)
              matrixCoefficient=matrixCoefficient*dampingMatrixCoefficient
              CALL DistributedMatrix_MatrixByVectorAdd(dampingDistributedMatrix,.FALSE.,incrementalVector, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,dynamicTempVector,err,error,*999)
            ENDIF
            CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
            IF(massMatrixNumber/=0.AND.dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
              NULLIFY(massMatrix)
              CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
              NULLIFY(massDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(massMatrix,massDistributedMatrix,err,error,*999)
              CALL EquationsMatrix_MatrixCoefficientGet(massMatrix,matrixCoefficient,err,error,*999)
              matrixCoefficient=matrixCoefficient*massMatrixCoefficient
              CALL DistributedMatrix_MatrixByVectorAdd(massDistributedMatrix,.FALSE.,incrementalVector, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,dynamicTempVector,err,error,*999)
            ENDIF
          ENDIF

          !Calculate the contributions from any linear matrices
          NULLIFY(linearMapping)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            NULLIFY(linearMatrices)
            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
            !Get the linear temporary vector
            NULLIFY(linearTempVector)
            CALL EquationsMatricesLinear_TempDistributedVectorGet(linearMatrices,linearTempVector,err,error,*999)
            !Initialise the temp vector to zero.
            CALL DistributedVector_AllValuesSet(linearTempVector,0.0_DP,err,error,*999)
            CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(linearMatrix)
              CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,linearMatrix,err,error,*999)
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,linearVariable,err,error,*999)
              !Check if the linear variable is on the left or right hand sides
              CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,linearVariable,solverMappingVariable, &
                & err,error,*999)
              IF(ASSOCIATED(solverMappingVariable)) THEN
                !The linear variable is a LHS variable to add it to the residual
                !Get the matrix coefficient
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)
                matrixCoefficient=matrixCoefficient*stiffnessMatrixCoefficient
                !Get the incremental parameter set.
                NULLIFY(incrementalVector)
                CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,incrementalVector, &
                  & err,error,*999)
                !Get the distributed matrix
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                !Add the linear matrix times the incremental vector
                CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,incrementalVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
              ENDIF
            ENDDO !equationsMatrixIdx
          ENDIF

          !Calculate the contribution from nonlinear residuals
          NULLIFY(nonlinearMapping)
          CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
          IF(ASSOCIATED(nonlinearMapping)) THEN
            NULLIFY(nonlinearMatrices)
            CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
            !Get the temp vector
            NULLIFY(nonlinearTempVector)
            CALL EquationsMatricesNonlinear_TempDistributedVectorGet(nonlinearMatrices,nonlinearTempVector,err,error,*999)
            !Initialise the temp vector
            CALL DistributedVector_AllValuesSet(nonlinearTempVector,0.0_DP,err,error,*999)
            !Loop over the residuals
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)            
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              !Loop over residual variables
              includeResidual=.FALSE.              
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
              DO variableIdx=1,numberOfResidualVariables
                NULLIFY(residualVariable)
                CALL EquationsMappingResidual_VariableGet(residualMapping,variableIdx,residualVariable,err,error,*999)
                !Check if the variable is on the LHS
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,residualVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) THEN
                  includeResidual=.TRUE.
                  EXIT
                ENDIF
              ENDDO !variableIdx
              IF(includeResidual) THEN
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                CALL EquationsMatricesResidual_VectorCoefficientGet(residualVector,vectorCoefficient,err,error,*999)
                vectorCoefficient=vectorCoefficient*currentFunctionFactor
                NULLIFY(residualDistributedVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & residualDistributedVector,err,error,*999)
                CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
                  & vectorCoefficient,residualDistributedVector,err,error,*999)
              ENDIF
            ENDDO !residualIdx
          ENDIF

          !Calculate the solver residual
          NULLIFY(equationsRowToSolverRowsMap)
          CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
            & equationsRowToSolverRowsMap,err,error,*999)

          !Couple the equations set vectors to the solver residual vector
          IF(ASSOCIATED(dynamicMapping)) &
            & CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
            & equationsRowToSolverRowsMap,1.0_DP,dynamicTempVector,err,error,*999)
          IF(ASSOCIATED(linearMapping)) &
            & CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
            & equationsRowToSolverRowsMap,1.0_DP,linearTempVector,err,error,*999)
          IF(ASSOCIATED(nonlinearMapping)) &
            & CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
            & equationsRowToSolverRowsMap,1.0_DP,nonlinearTempVector,err,error,*999)

          IF(.NOT.dynamicSolver%solverInitialised) THEN
            IF(ASSOCIATED(nonlinearMapping)) THEN              
              !Copy current residual i.e., residual at time zero, to previous residuals to initialise
              DO residualIdx=1,numberOfResiduals
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                NULLIFY(currentResidualVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & currentResidualVector,err,error,*999)
                NULLIFY(previousResidualVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
                  & previousResidualVector,err,error,*999)
                CALL DistributedVector_VectorCopy(previousResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                  & currentResidualVector,err,error,*999)
                IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
                  NULLIFY(previous2ResidualVector)
                  CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                    & previous2ResidualVector,err,error,*999)
                  CALL DistributedVector_VectorCopy(previous2ResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                    & currentResidualVector,err,error,*999)
                  IF(dynamicSolver%degree>SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                    NULLIFY(previous3ResidualVector)
                    CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS3_VECTOR, &
                      & previous3ResidualVector,err,error,*999)
                    CALL DistributedVector_VectorCopy(previous3ResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP, &
                      & currentResidualVector,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !residualIdx
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx

        !Loop over the interface conditions
        DO interfaceConditionIdx=1,numberOfInterfaceConditions
          NULLIFY(interfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          NULLIFY(interfaceEquations)
          CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          NULLIFY(interfaceMapping)
          CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
          NULLIFY(interfaceMatrices)
          CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
          CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
          CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
          IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD) numberOfInterfaceMatrices=numberOfInterfaceMatrices-1
          
          NULLIFY(interfaceConditionToSolverMatricesMap)
          CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
            & interfaceConditionToSolverMatricesMap,err,error,*999)
          !Get the Lagrange vector
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          NULLIFY(lagrangeDistributedVector)
          CALL FieldVariable_ParameterSetVectorGet(lagrangeVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,lagrangeDistributedVector, &
            & err,error,*999)

          !Calculate the contributions from any interface matrices
          DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
            !Calculate the interface matrix-Lagrange vector product residual contribution
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
            CALL InterfaceMatrix_MatrixCoefficientGet(interfaceMatrix,matrixCoefficients(1),err,error,*999)
            NULLIFY(interfaceTempVector)
            CALL InterfaceMatrix_TempDistributedVectorGet(interfaceMatrix,interfaceTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(interfaceTempVector,0.0_DP,err,error,*999)
            !Get the interface matrix mappings to the solver matrix
            NULLIFY(interfaceMatrixToSolverMatricesMap)
            CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
            NULLIFY(interfaceRowToSolverRowsMap)
            CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
              & interfaceRowToSolverRowsMap,err,error,*999)
            CALL InterfaceMatrix_TimeDependenceTypeGet(interfaceMatrix,timeDependenceType,err,error,*999)
            SELECT CASE(timeDependenceType)
            CASE(INTERFACE_MATRIX_STATIC)
              matrixCoefficients(1)=matrixCoefficients(1)*stiffnessMatrixCoefficient
            CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
              matrixCoefficients(1)=matrixCoefficients(1)*dampingMatrixCoefficient
            CASE(INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
              matrixCoefficients(1)=matrixCoefficients(1)*massMatrixCoefficient
            CASE DEFAULT
              localError="The time dependence type of "//TRIM(NumberToVString(timeDependenceType,"*",err,error))// &
                & " for interface matrix "//TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))// &
                & " of interface condition number "//TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))// &
                & " is invalid or not implemented."
              CALL FlagError(localError,err,error,*999)
            END SELECT

            NULLIFY(interfaceDistributedMatrix)
            CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
            
            CALL DistributedMatrix_MatrixByVectorAdd(interfaceDistributedMatrix,.FALSE.,lagrangeDistributedVector, &
              & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficients(1),interfaceTempVector,err,error,*999)
            
            !Add interface matrix residual contribution to the solver residual
            CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
              & interfaceRowToSolverRowsMap,1.0_DP,interfaceTempVector,err,error,*999)
            
            CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
            IF(hasTranspose) THEN
              CALL InterfaceMatrix_TransposeMatrixCoefficientGet(interfaceMatrix,matrixCoefficients(2),err,error,*999)
              CALL InterfaceMatrix_TransposeTimeDependenceTypeGet(interfaceMatrix,transposeTimeDependenceType,err,error,*999)
              SELECT CASE(TransposeTimeDependenceType)
              CASE(INTERFACE_MATRIX_STATIC)
                matrixCoefficients(2)=matrixCoefficients(2)*stiffnessMatrixCoefficient
              CASE(INTERFACE_MATRIX_FIRST_ORDER_DYNAMIC)
                matrixCoefficients(2)=matrixCoefficients(2)*dampingMatrixCoefficient
              CASE(INTERFACE_MATRIX_SECOND_ORDER_DYNAMIC)
                matrixCoefficients(2)=matrixCoefficients(2)*massMatrixCoefficient
              CASE DEFAULT
                localError="The transpose time dependence type of "// &
                  & TRIM(NumberToVString(transposeTimeDependenceType,"*",err,error))//" for interface matrix "// &
                  & TRIM(NumberToVString(interfaceMatrixIdx,"*",err,error))//" of interface condition number "// &
                  & TRIM(NumberToVString(interfaceCondition%userNumber,"*",err,error))//" is invalid or not implemented."
                CALL FlagError(localError,err,error,*999)
              END SELECT
           
              !Calculate the transposed interface matrix-dependent variable product residual contribution
              NULLIFY(interfaceMatrixToVarMap)
              CALL InterfaceMapping_InterfaceMatrixToVarMapGet(interfaceMapping,interfaceMatrixIdx,interfaceMatrixToVarMap, &
                & err,error,*999)          
              NULLIFY(dependentVariable)
              CALL InterfaceMappingIMToVMap_VariableGet(interfaceMatrixToVarMap,dependentVariable,err,error,*999)
              !Get the temp vector
              NULLIFY(transposeInterfaceTempVector)
              CALL InterfaceMatrix_TempTransposeDistributedVectorGet(interfaceMatrix,transposeInterfaceTempVector,err,error,*999)
              !Initialise the interface temporary vector to zero
              CALL DistributedVector_AllValuesSet(transposeInterfaceTempVector,0.0_DP,err,error,*999)
              NULLIFY(dependentDistributedVector)
              CALL FieldVariable_ParameterSetVectorGet(dependentVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE, &
              & dependentDistributedVector,err,error,*999)
              NULLIFY(transposeInterfaceDistributedMatrix)
              CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeInterfaceDistributedMatrix, &
                & err,error,*999)
            
              CALL DistributedMatrix_MatrixByVectorAdd(transposeInterfaceDistributedMatrix,.FALSE.,dependentDistributedVector, &
                & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficients(2),transposeInterfaceTempVector,err,error,*999)
              
              !Add interface matrix residual contribution to the solver residual.
              !The number of columns in the interface matrix is equivalent to the
              !number of rows of the transposed interface matrices
              NULLIFY(interfaceColToSolverRowsMap)
              CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
                & interfaceColToSolverRowsMap,err,error,*999)
              
              CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                & interfaceColToSolverRowsMap,1.0_DP,transposeInterfaceTempVector,err,error,*999)
            ENDIF !has transpose
          ENDDO !interfaceMatrixIdx
          IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD) THEN
            CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,interfaceMatrixIdx,err,error,*999)         
            !Calculate the Lagrange-Lagrange vector product residual contribution from the penalty term
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
            NULLIFY(interfaceDistributedMatrix)
            CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
            NULLIFY(interfaceTempVector)
            CALL InterfaceMatrix_TempDistributedVectorGet(interfaceMatrix,interfaceTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(interfaceTempVector,0.0_DP,err,error,*999)
!!TODO: should this be the field values or the incremental vector???
            NULLIFY(lagrangeDistributedVector)
            CALL FieldVariable_ParameterSetVectorGet(lagrangeVariable,FIELD_VALUES_SET_TYPE,lagrangeDistributedVector, &
              & err,error,*999)
            
            CALL DistributedMatrix_MatrixByVectorAdd(interfaceDistributedMatrix,.FALSE.,lagrangeDistributedVector, &
              & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP,interfaceTempVector,err,error,*999)

            !Add interface matrix residual contribution to the solver residual
            NULLIFY(interfaceMatrixToSolverMatricesMap)
            CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
            NULLIFY(interfaceRowToSolverRowsMap)
            CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
              & interfaceRowToSolverRowsMap,err,error,*999)
            
            CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
              & interfaceRowToSolverRowsMap,1.0_DP,interfaceTempVector,err,error,*999)
          ENDIF
        ENDDO !interfaceConditionIdx
        !
        !Start the update the solver residual vector values
        CALL DistributedVector_UpdateStart(solverResidualVector,err,error,*999)
        !Restore the solver residual check data
        CALL DistributedVector_DataRestore(solverResidualVector,solverResidualCheckData,err,error,*999)
      ENDIF !Update residual
      
      IF(ASSOCIATED(solverResidualVector)) CALL DistributedVector_UpdateFinish(solverResidualVector,err,error,*999)
      
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
          & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
        CALL Profiling_TimingsOutput(1,"Solver residual assembly",userElapsed,systemElapsed,err,error,*999)
      ENDIF
      
    ENDIF !Calculate solver residual
    
    IF(dynamicSolver%solverInitialised) THEN
      !Set the first part of the next time step. Note that we do not have to add in the previous time value as it is
      !already there from when we copied the values to the previous time step.
      !Loop over the equations sets
      IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
          IF(ASSOCIATED(dynamicMapping)) THEN
            NULLIFY(dynamicVariable)
            CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
            SELECT CASE(dynamicSolver%degree)
            CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
              !Do nothing. Increment will be added after the solve.
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
              CALL FieldVariable_ParameterSetsAdd(dynamicVariable,firstUpdateFactor,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
              CALL FieldVariable_ParameterSetsAdd(dynamicVariable,[firstUpdateFactor,secondUpdateFactor], &
                & [FIELD_PREVIOUS_VELOCITY_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE],FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetsAdd(dynamicVariable,firstUpdateFactor,FIELD_PREVIOUS_ACCELERATION_SET_TYPE, &
                & FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
            CASE DEFAULT
              localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDDO !equationsSetIdx
      ENDIF
    ENDIF
      
    !If required output the solver matrices
    IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) THEN
      CALL SolverMatrices_Output(GENERAL_OUTPUT_TYPE,selectionType,solverMatrices,err,error,*999)
    ENDIF
    
    EXITS("Solver_DynamicAssemble")
    RETURN
999 ERRORSEXITS("Solver_DynamicAssemble",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_DynamicAssemble

  !
  !================================================================================================================================
  !

  !>Assembles the solver matrices and rhs from the static equations.
  SUBROUTINE Solver_StaticAssemble(solver,selectionType,err,error,*)

    !Argument variable
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(IN) :: selectionType !<The type of matrix selection \see SolverMatricesRoutines_SelectMatricesTypes,SolverMatricesRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dependentVariableType,dirichletIdx,equationsColumnNumber,equationsMatrixIdx,equationsMatrixIdx2, &
      & equationsMatrixNumber,equationsRowNumber,equationsSetIdx,interfaceConditionIdx,interfaceConditionMethod, &
      & interfaceMatrixIdx,interfaceVariableType,jacobianMatrixIdx,lhsBoundaryCondition,lhsBoundaryFinish,lhsGlobalDOF, &
      & lhsVariableDOF,linearMatrixIdx,linearMatrixNumber,linearVariableIdx,linearVariableType,numberOfDirichletConditions, &
      & numberOfEquationsMatrices,numberOfEquationsSets,numberOfInterfaceConditions,numberOfInterfaceMatrices, &
      & numberOfJacobianMatrices,numberOfLinearMatrices,numberOfLinearVariables,numberOfResiduals,numberOfResidualVariables, &
      & numberOfRows,numberOfSolverMatrices,numberOfSources,penaltyMatrixIdx,residualIdx,residualVariableIdx, &
      & rhsBoundaryCondition,rhsGlobalDOF,rhsVariableDOF,rhsVariableType,rowCondition,solverMatrixIdx,solverRowNumber, &
      & sourceIdx,totalNumberOfRows,variableDOF,variableGlobalDOF,variableBoundaryCondition,variableIdx,variableType
    INTEGER(INTG), POINTER :: equationsRowToLHSDOFMap(:),equationsRowToRHSDOFMap(:)
    REAL(SP) :: systemElapsed,systemTime1(1),systemTime2(1),userElapsed,userTime1(1),userTime2(1)
    REAL(DP) :: alphaValue,currentRHSValue,dependentValue,dofValue,linearValue,linearValueSum,matrixCoefficient, &
      & matrixCoefficients(2),nonlinearValue,residualCoefficient,residualValue,rhsIntegratedValue,rhsValue,solverRHSValue, &
      & sourceCoefficient,sourceValue
    REAL(DP), POINTER :: checkData(:),checkData2(:),checkData3(:),checkData4(:),matrixCheckData(:),rhsIntegratedParameters(:), &
      & rhsParameters(:),solverResidualCheckData(:),solverRHSCheckData(:)
    TYPE(RealDPPtrType), ALLOCATABLE :: dependentParameters(:)
    LOGICAL :: hasIntegratedValues,hasTranspose,rhsLinearMatrix,rhsResidual,solverResidual,subtractFixedBCsFromResidual, &
      & updateMatrix,updateResidual,updateRHS
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: dependentBoundaryConditions,lhsBoundaryConditionsVariable, &
      & rhsBoundaryConditionsVariable
    TYPE(BoundaryConditionsRowVariableType), POINTER :: lhsBoundaryConditionsRowVariable
    TYPE(DistributedMatrixType), POINTER :: interfaceDistributedMatrix,jacobianDistributedMatrix,linearDistributedMatrix, &
      & previousSolverDistributedMatrix,transposeDistributedMatrix,solverDistributedMatrix
    TYPE(DistributedVectorType), POINTER :: currentRHSVector,currentSourceVector,dependentDistributedVector,dependentVector, &
      & distributedSourceVector,equationsRHSVector,interfaceRHSDistributedVector,interfaceTempVector,lagrangeVector, &
      & linearTempVector,nonlinearTempVector,residualDistributedVector,solverResidualVector,solverRHSVector,sourcesTempVector
    TYPE(DomainMappingType), POINTER :: lhsDomainMapping,rhsDomainMapping,variableDomainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix,linearMatrix
    TYPE(EquationsMatrixToSolverMatrixMapType), POINTER :: equationsMatrixToSolverMatrixMap,linearMatrixToSolverMatrixMap
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,interfaceVariable,lagrangeVariable,lhsVariable,linearVariable, &
      & residualVariable,rhsVariable
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceLagrangeType), POINTER :: interfaceLagrange
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMappingRHSType), POINTER :: interfaceRHSMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMatrixToSolverMatricesMapType), POINTER :: interfaceMatrixToSolverMatricesMap
    TYPE(InterfaceMatrixToSolverMatrixMapType), POINTER :: interfaceMatrixToSolverMatrixMap
    TYPE(InterfaceRHSType), POINTER :: interfaceRHSVector
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(JacobianMatrixToSolverMatrixMapType), POINTER :: jacobianMatrixToSolverMatrixMap
    TYPE(MatrixRowColCouplingType), POINTER :: equationsColToSolverColsMap(:),equationsRowToSolverRowsMap(:), &
      & interfaceColToSolverColsMap(:),interfaceColToSolverRowsMap(:),interfaceRowToSolverColsMap(:), &
      & interfaceRowToSolverRowsMap(:),jacobianColToSolverColsMap(:)
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VarToEquationsMatricesMapType), POINTER :: linearVarToEquationsMatricesMap
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("Solver_StaticAssemble",err,error,*999)
  
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
    CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
    
    !Assemble the solver matrices
    NULLIFY(previousSolverDistributedMatrix)
    
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
      !Assemble solver matrices
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF
      !Loop over the solver matrices
      CALL SolverMapping_NumberOfSolverMatricesGet(solverMapping,numberOfSolverMatrices,err,error,*999)
      DO solverMatrixIdx=1,numberOfSolverMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        CALL SolverMatrix_UpdateMatrixGet(solverMatrix,updateMatrix,err,error,*999)
        IF(updateMatrix) THEN
          NULLIFY(solverDistributedMatrix)
          CALL SolverMatrix_SolverDistributedMatrixGet(solverMatrix,solverDistributedMatrix,err,error,*999)
          !Initialise matrix to zero
          CALL DistributedMatrix_AllValuesSet(solverDistributedMatrix,0.0_DP,err,error,*999)
          !Get the check data
          NULLIFY(matrixCheckData)
          CALL DistributedMatrix_DataGet(solverDistributedMatrix,matrixCheckData,err,error,*999)
          
          !Loop over the equations sets
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSetToSolverMatricesMap)
            CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
              & err,error,*999)
            NULLIFY(equationsRowToSolverRowsMap)
            CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
              & equationsRowToSolverRowsMap,err,error,*999)
            NULLIFY(equationsMatricesToSolverMatrixMap)
            CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
              & equationsMatricesToSolverMatrixMap,err,error,*999)
            IF(selectionType==SOLVER_MATRICES_ALL.OR. &
              & selectionType==SOLVER_MATRICES_LINEAR_ONLY) THEN
              CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap, &
                & numberOfLinearMatrices,err,error,*999)
              DO linearMatrixIdx=1,numberOfLinearMatrices
                NULLIFY(linearMatrixToSolverMatrixMap)
                CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                  & linearMatrixIdx,linearMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(equationsColToSolverColsMap)
                CALL SolverMappingEMToSMMap_EquationsColToSolverColsMapGet(linearMatrixToSolverMatrixMap, &
                  & equationsColToSolverColsMap,err,error,*999)
                NULLIFY(linearMatrix)
                CALL SolverMappingEMToSMMap_EquationsMatrixGet(linearMatrixToSolverMatrixMap,linearMatrix,err,error,*999)
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)
                CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & equationsRowToSolverRowsMap,equationsColToSolverColsMap,matrixCoefficient,linearDistributedMatrix, &
                  & .FALSE.,err,error,*999)
              ENDDO !equationsMatrixIdx
            ENDIF
            IF(selectionType==SOLVER_MATRICES_ALL.OR. &
              & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
              & selectionType==SOLVER_MATRICES_JACOBIAN_ONLY) THEN
              !Now set the values from the equations Jacobian
              CALL SolverMappingEMSToSMMap_NumberOfJacobianMatricesGet(equationsMatricesToSolverMatrixMap, &
                & numberOfJacobianMatrices,err,error,*999)
              DO jacobianMatrixIdx=1,numberOfJacobianMatrices
                NULLIFY(jacobianMatrixToSolverMatrixMap)
                CALL SolverMappingEMSToSMMap_JacobianMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                  & jacobianMatrixIdx,jacobianMatrixToSolverMatrixMap,err,error,*999)
                NULLIFY(jacobianColToSolverColsMap)
                CALL SolverMappingJMToSMMap_JacobianColToSolverColsMapGet(jacobianMatrixToSolverMatrixMap, &
                  & jacobianColToSolverColsMap,err,error,*999)
                NULLIFY(jacobianMatrix)
                CALL SolverMappingJMToSMMap_JacobianMatrixGet(jacobianMatrixToSolverMatrixMap,jacobianMatrix,err,error,*999)
                NULLIFY(jacobianDistributedMatrix)
                CALL JacobianMatrix_DistributedMatrixGet(jacobianMatrix,jacobianDistributedMatrix,err,error,*999)
                CALL JacobianMatrix_MatrixCoefficientGet(jacobianMatrix,matrixCoefficient,err,error,*999)
                CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & equationsRowToSolverRowsMap,jacobianColToSolverColsMap,matrixCoefficient,jacobianDistributedMatrix, &
                  & .FALSE.,err,error,*999)
              ENDDO !jacobianMatrixIdx
            ENDIF
             
          ENDDO !equationsSetIdx
          
          !Loop over any interface conditions
          DO interfaceConditionIdx=1,numberOfInterfaceConditions
            NULLIFY(interfaceConditionToSolverMatricesMap)
            CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
              & interfaceConditionToSolverMatricesMap,err,error,*999)
            NULLIFY(interfaceColToSolverRowsMap)
            CALL SolverMappingICToSMSMap_InterfaceColToSolverRowsMapGet(interfaceConditionToSolverMatricesMap, &
              & interfaceColToSolverRowsMap,err,error,*999)
             NULLIFY(interfaceMatricesToSolverMatrixMap)
            CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap, &
              & solverMatrixIdx,interfaceMatricesToSolverMatrixMap,err,error,*999)
            NULLIFY(interfaceColToSolverColsMap)
            CALL SolverMappingIMSToSMMap_InterfaceColToSolverColsMapGet(interfaceMatricesToSolverMatrixMap, &
              & interfaceColToSolverColsMap,err,error,*999)
            !Loop over the interface matrices mapped to this solver matrix
            CALL SolverMappingIMSToSMMap_NumberOfInterfaceMatricesGet(interfaceMatricesToSolverMatrixMap, &
              & numberOfInterfaceMatrices,err,error,*999)
            DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
              NULLIFY(interfaceMatrixToSolverMatricesMap)
              CALL SolverMappingICToSMSMap_InterfaceMatrixToSolverMatricesMapGet(interfaceConditionToSolverMatricesMap, &
                & interfaceMatrixIdx,interfaceMatrixToSolverMatricesMap,err,error,*999)
              NULLIFY(interfaceRowToSolverRowsMap)
              CALL SolverMappingIMToSMSMap_InterfaceRowToSolverRowsMapGet(interfaceMatrixToSolverMatricesMap, &
                & interfaceRowToSolverRowsMap,err,error,*999)
              NULLIFY(interfaceMatrixToSolverMatrixMap)
              CALL SolverMappingIMSToSMMap_InterfaceMatrixToSolverMatrixMapGet(interfaceMatricesToSolverMatrixMap, &
                & interfaceMatrixIdx,interfaceMatrixToSolverMatrixMap,err,error,*999)
              NULLIFY(interfaceMatrix)
              CALL SolverMappingIMToSMMap_InterfaceMatrixGet(interfaceMatrixToSolverMatrixMap,interfaceMatrix,err,error,*999)
              NULLIFY(interfaceDistributedMatrix)
              CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
              CALL InterfaceMatrix_MatrixCoefficientGet(interfaceMatrix,matrixCoefficients(1),err,error,*999)
              CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                & interfaceRowToSolverRowsMap,interfaceColToSolverColsMap,matrixCoefficients(1),interfaceDistributedMatrix, &
                & .FALSE.,err,error,*999)
              CALL InterfaceMatrix_HasTransposeGet(interfaceMatrix,hasTranspose,err,error,*999)
              IF(hasTranspose) THEN
                NULLIFY(interfaceRowToSolverColsMap)
                CALL SolverMappingIMToSMMap_InterfaceRowToSolverColsMapGet(interfaceMatrixToSolverMatrixMap, &
                  & interfaceRowToSolverColsMap,err,error,*999)
                NULLIFY(transposeDistributedMatrix)
                CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,transposeDistributedMatrix,err,error,*999)
                CALL InterfaceMatrix_TransposeMatrixCoefficientGet(interfaceMatrix,matrixCoefficients(2),err,error,*999)
                CALL DistributedMatrix_MatrixCoupleAdd(solverDistributedMatrix,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & interfaceColToSolverRowsMap,interfaceRowToSolverColsMap,matrixCoefficients(2),transposeDistributedMatrix, &
                  & .FALSE.,err,error,*999)
              ENDIF
            ENDDO !interfaceMatrixIdx
          ENDDO !interfaceConditionIdx
          
          !Update the solver matrix values
          CALL DistributedMatrix_UpdateStart(solverDistributedMatrix,err,error,*999)
          IF(ASSOCIATED(previousSolverDistributedMatrix)) &
            & CALL DistributedMatrix_UpdateFinish(previousSolverDistributedMatrix,err,error,*999)                      
          previousSolverDistributedMatrix=>solverDistributedMatrix
          
        ENDIF !Update matrix
      ENDDO !solverMatrixIdx
      
      IF(ASSOCIATED(previousSolverDistributedMatrix)) &
        & CALL DistributedMatrix_UpdateFinish(previousSolverDistributedMatrix,err,error,*999)
      
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
          & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
        CALL Profiling_TimingsOutput(1,"Solver matrices assembly",userElapsed,systemElapsed,err,error,*999)
      ENDIF
      
    ENDIF
        
    NULLIFY(solverRHSVector)
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_LINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
      !Assemble rhs vector
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF
      
      CALL SolverMatrices_UpdateRHSGet(solverMatrices,updateRHS,err,error,*999)
      IF(updateRHS) THEN
        NULLIFY(solverRHSVector)
        CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,solverRHSVector,err,error,*999)
        !Initialise the RHS to zero
        CALL DistributedVector_AllValuesSet(solverRHSVector,0.0_DP,err,error,*999)
        !Get the solver RHS check data                  
        NULLIFY(solverRHSCheckData)
        CALL DistributedVector_DataGet(solverRHSVector,solverRHSCheckData,err,error,*999)
        !Get the boundary conditions
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        
        subtractFixedBCsFromResidual=.FALSE.
        IF(selectionType==SOLVER_MATRICES_ALL.OR. &
          & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
          & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
          IF(solverMatrices%updateResidual) subtractFixedBCsFromResidual=.TRUE.
        ENDIF
        
        !Loop over the equations sets
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(equationsSetToSolverMatricesMap)
          CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
            & err,error,*999)
          NULLIFY(equationsRowToSolverRowsMap)
          CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
            & equationsRowToSolverRowsMap,err,error,*999)
          
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(lhsMapping)
          CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
          NULLIFY(lhsVariable)
          CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
          NULLIFY(equationsRowToLHSDOFMap)
          CALL EquationsMappingLHS_EquationsRowTOLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*999)
          NULLIFY(lhsDomainMapping)
          CALL FieldVariable_DomainMappingGet(lhsVariable,lhsDomainMapping,err,error,*999)
          CALL DomainMapping_BoundaryFinishGet(lhsDomainMapping,lhsBoundaryFinish,err,error,*999)
          CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*999)
          CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*999)
          
          !Calculate the contributions from any nonlinear residuals
          NULLIFY(nonlinearTempVector)
          NULLIFY(nonlinearMapping)
          NULLIFY(nonlinearMatrices)
          CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
          IF(ASSOCIATED(nonlinearMapping)) THEN
            CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
            CALL EquationsMatricesNonlinear_TempDistributedVectorGet(nonlinearMatrices,nonlinearTempVector,err,error,*999)
            CALL DistributedVector_AllValuesSet(nonlinearTempVector,0.0_DP,err,error,*999)
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              rhsResidual=.TRUE.
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(solverMatrixToEquationsMap)
                CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
                  & err,error,*999)
                NULLIFY(solverMappingVariables)
                CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
                DO residualVariableIdx=1,numberOfResidualVariables
                  CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
                  NULLIFY(solverMappingVariable)
                  CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,residualVariable,solverMappingVariable, &
                    & err,error,*999)
                  IF(ASSOCIATED(solverMappingVariable)) rhsResidual=.FALSE.
                ENDDO !residualVariableIdx
              ENDDO !solverMatrixIdx
              IF(rhsResidual) THEN
                !The residual does not have any variables involved in the solver matrices so take it over to the RHS
                CALL EquationsMappingResidual_VectorCoefficientGet(residualMapping,residualCoefficient,err,error,*999)
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                NULLIFY(residualDistributedVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & residualDistributedVector,err,error,*999)
                CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & residualCoefficient,residualDistributedVector,err,error,*999)
              ENDIF
            ENDDO !residualIdx
          ENDIF !nonlinear mapping

          !Calculate the contributions from any linear matrices
          NULLIFY(linearTempVector)
          NULLIFY(linearMapping)
          NULLIFY(linearMatrices)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
            CALL EquationsMatricesLinear_TempDistributedVectorGet(linearMatrices,linearTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(linearTempVector,0.0_DP,err,error,*999)                  
            CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
            DO equationsMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,linearVariable,err,error,*999)
              rhsLinearMatrix=.TRUE.
              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(solverMatrixToEquationsMap)
                CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
                  & err,error,*999)
                NULLIFY(solverMappingVariables)
                CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,linearVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) rhsLinearMatrix=.FALSE.
              ENDDO !solverMatrixIdx
              IF(rhsLinearMatrix) THEN
                !Linear matrix variable is not on the LHS so take the matrix times vector over to the RHS.
                NULLIFY(linearMatrix)
                CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,linearMatrix,err,error,*999)
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)
                NULLIFY(dependentDistributedVector)
                CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_VALUES_SET_TYPE,dependentDistributedVector, &
                  & err,error,*999)
                CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,dependentDistributedVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
              ENDIF !rhs linear matrix
            ENDDO !equationsMatrixIdx
          ENDIF !linear mapping

          !Add in any source vectors
          NULLIFY(sourcesTempVector)
          NULLIFY(sourcesMapping)
          NULLIFY(sourceVectors)
          CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
          IF(ASSOCIATED(sourcesMapping)) THEN
            CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
            CALL EquationsMatricesSources_TempDistributedVectorGet(sourceVectors,sourcesTempVector,err,error,*999)
            CALL DistributedVector_AllValuesSet(sourcesTempVector,0.0_DP,err,error,*999)
            CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
            DO sourceIdx=1,numberOfSources
              NULLIFY(sourceMapping)
              CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,sourceIdx,sourceMapping,err,error,*999)
              CALL EquationsMappingSource_VectorCoefficientGet(sourceMapping,sourceCoefficient,err,error,*999)
              NULLIFY(sourceVector)
              CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
              NULLIFY(distributedSourceVector)
              CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                & distributedSourceVector,err,error,*999)
              CALL DistributedVector_VectorAdd(sourcesTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,sourceCoefficient, &
                & currentSourceVector,err,error,*999)
            ENDDO !sourceIdx
          ENDIF !source mapping
          
          NULLIFY(rhsVariable)
          NULLIFY(rhsMapping)
          NULLIFY(rhsVector)
          NULLIFY(currentRHSVector)
          NULLIFY(rhsBoundaryConditionsVariable)
          NULLIFY(rhsIntegratedParameters)
          CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
            CALL FieldVariable_ParameterSetCreated(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,hasIntegratedValues, &
              & err,error,*999)
            CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
            CALL EquationsMatricesRHS_DistributedVectorGet(rhsVector,EQUATIONS_MATRICES_CURRENT_VECTOR,currentRHSVector, &
              & err,error,*999)
            CALL FieldVariable_ParameterSetDataGet(rhsVariable,FIELD_VALUES_SET_TYPE,rhsParameters,err,error,*999) 
            IF(hasIntegratedValues) THEN
              !Update RHS field by integrating any point Neumann conditions
              CALL BoundaryConditions_VariableGet(boundaryConditions,rhsVariable,rhsBoundaryConditionsVariable,err,error,*999)
!!TODO: NEED TO FIX INTEGRATED FLUX MATRIX CALCULATION. THE MATRIX CURRENTLY USES THE RHS FOR THE ROWS BUT IT SHOULD USE THE
!!      LHS. IT SHOULD BE POSSIBLE TO USE A DIFFERENT BASIS FOR THE FLUX INTERPOLATION THAN FOR THE DEPENDENT INTERPOLATION.
              CALL BoundaryConditionsVariable_NeumannIntegrate(rhsBoundaryConditionsVariable,err,error,*999)
              CALL FieldVariable_ParameterSetDataGet(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhsIntegratedParameters, &
                & err,error,*999) 
            ENDIF
            NULLIFY(equationsRowToRHSDOFMap)
            CALL EquationsMappingRHS_EquationsRowToRHSDOFMapGet(rhsMapping,equationsRowToRHSDOFMap,err,error,*999)
          ENDIF
          
          NULLIFY(lhsBoundaryConditionsRowVariable)
          CALL BoundaryConditions_RowVariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsRowVariable,err,error,*999)
          NULLIFY(lhsBoundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsVariable,err,error,*999)
          CALL BoundaryCOnditionsVariable_NumberOfDirichletConditionsGet(lhsBoundaryConditionsVariable, &
            & numberOfDirichletConditions,err,error,*999)
          NULLIFY(dirichletBoundaryConditions)
          IF(numberOfDirichletConditions>0) CALL BoundaryConditionsVariable_DirichletConditionsExists( &
            & lhsBoundaryConditionsVariable,dirichletBoundaryConditions,err,error,*999)
          NULLIFY(neumannBoundaryConditions)
          CALL BoundaryConditionsVariable_NeumannConditionsExists(lhsBoundaryConditionsVariable, &
            & neumannBoundaryConditions,err,error,*999)

          !Loop over the rows in the equations set
          DO equationsRowNumber=1,totalNumberOfRows
            CALL BoundaryConditionsRowVariable_RowConditionTypeGet(lhsBoundaryConditionsRowVariable,equationsRowNumber, &
              & rowCondition,err,error,*999)

            !Get the nonlinear contribution to the RHS values               
            IF(ASSOCIATED(nonlinearTempVector)) THEN
              CALL DistributedVector_ValuesGet(nonlinearTempVector,equationsRowNumber,nonlinearValue,err,error,*999)
            ELSE
              nonlinearValue=0.0_DP
            ENDIF
            !Get the linear contribution to the RHS values               
            IF(ASSOCIATED(linearTempVector)) THEN
              CALL DistributedVector_ValuesGet(linearTempVector,equationsRowNumber,linearValue,err,error,*999)
            ELSE
              linearValue=0.0_DP
            ENDIF
            !Get the source contribution to the RHS values               
            IF(ASSOCIATED(sourcesTempVector)) THEN
              CALL DistributedVector_ValuesGet(sourcesTempVector,equationsRowNumber,sourceValue,err,error,*999)                
            ELSE
              sourceValue=0.0_DP
            ENDIF

            solverRHSValue=nonlinearValue+linearValue+sourceValue

            IF(ASSOCIATED(rhsMapping)) THEN
              rhsVariableDOF=equationsRowToRHSDOFMap(equationsRowNumber)
              IF(hasIntegratedValues) THEN
                !Add any Neumann integrated values, b = f + N q
                CALL DistributedVector_ValuesAdd(currentRHSVector,equationsRowNumber,rhsIntegratedParameters(rhsVariableDOF), &
                  & err,error,*999)
              ELSE
                CALL DistributedVector_ValuesSet(currentRHSVector,equationsRowNumber,rhsParameters(rhsVariableDOF), &
                  & err,error,*999)
              ENDIF
              CALL DistributedVector_ValuesGet(currentRHSVector,equationsRowNumber,rhsValue,err,error,*999)
            ELSE
              rhsValue=0.0_DP
            ENDIF

            !Get the dynamic contribution to the the RHS values
            lhsVariableDOF=equationsRowToLHSDOFMap(equationsRowNumber)
            CALL DomainMapping_LocalToGlobalGet(lhsDomainMapping,lhsVariableDOF,lhsGlobalDOF,err,error,*999)
            CALL BoundaryConditionsVariable_DOFTypeGet(lhsBoundaryConditionsVariable,lhsGlobalDOF,lhsBoundaryCondition, &
              & err,error,*999)

            !Apply boundary conditions
            SELECT CASE(lhsBoundaryCondition)
            CASE(BOUNDARY_CONDITION_FREE_ROW)
              solverRHSValue=solverRHSValue+rhsValue
              !Loop over the solver rows associated with this equations set row
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,solverRHSValue,err,error,*999)
            CASE(BOUNDARY_CONDITION_DIRICHLET_ROW)
!!WHY ARE WE ADDING TO THE RHS ROW IF WE ARE GOING TO ELIMINATE IT?
              !Get the equations RHS values
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,rhsValue,err,error,*999)

              DO solverMatrixIdx=1,numberOfSolverMatrices
                NULLIFY(solverMatrixToEquationsMap)
                CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap, &
                  & err,error,*999)
                NULLIFY(solverMappingVariables)
                CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,lhsVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) THEN
                  IF(ASSOCIATED(linearMapping)) THEN
                    NULLIFY(equationsMatricesToSolverMatrixMap)
                    CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap, &
                      & solverMatrixIdx,equationsMatricesToSolverMatrixMap,err,error,*999)
                    CALL SolverMappingEMSToSMMap_NumberOfLinearMatricesGet(equationsMatricesToSolverMatrixMap, &
                      & numberOfLinearMatrices,err,error,*999)
                    DO linearMatrixIdx=1,numberOfLinearMatrices
                      NULLIFY(linearMatrixToSolverMatrixMap)
                      CALL SolverMappingEMSToSMMap_LinearMatrixToSolverMatrixMapGet(equationsMatricesToSolverMatrixMap, &
                        & linearMatrixIdx,linearMatrixToSolverMatrixMap,err,error,*999)                    
                      NULLIFY(linearMatrix)
                      CALL SolverMappingEMToSMMap_EquationsMatrixGet(linearMatrixToSolverMatrixMap,linearMatrix,err,error,*999)
                      CALL EquationsMatrix_MatrixNumberGet(linearMatrix,linearMatrixNumber,err,error,*999)
                      NULLIFY(linearVariable)
                      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixNumber,linearVariable, &
                        & err,error,*999)
                      IF(ASSOCIATED(linearVariable,lhsVariable)) THEN
                        CALL FieldVariable_ParameterSetGetLocalDOF(linearVariable,FIELD_VALUES_SET_type,lhsVariableDOF, &
                          & dofValue,err,error,*999)
                        CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)
                        dofValue=dofValue*matrixCoefficient
                        NULLIFY(linearDistributedMatrix)
                        CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
!!WHAT IS THIS TRYING TO DO?
                        DO dirichletIdx=1,numberOfDirichletConditions
                          IF(dirichletBoundaryConditions%dirichletDOFIndices(dirichletIdx)==lhsGlobalDOF) EXIT
                        ENDDO !dirichletIdx
                        CALL DistributedMatrix_MatrixColumnAdd(linearDistributedMatrix,.FALSE.,lhsGlobalDOF, &
                          & equationsRowToSolverRowsMap,-1.0_DP*dofValue,solverRHSVector,err,error,*999)

                      ENDIF
                    ENDDO !linearMatrixIdx
                  ENDIF !linear mapping
                ENDIF !solver mapping variable
              ENDDO !solverMatrixIdx

            CASE(BOUNDARY_CONDITION_NEUMANN_ROW)
              !Set Neumann boundary conditions
              solverRHSValue=solverRHSValue+rhsValue
              !Loop over the solver rows associated with this equations set row
              CALL DistributedVector_VectorRowCoupleAdd(solverRHSVector,equationsRowToSolverRowsMap(equationsRowNumber), &
                & -1.0_DP,solverRHSValue,err,error,*999)
            CASE(BOUNDARY_CONDITION_ROBIN_ROW)
              !Set Robin boundary conditions
              CALL FlagError("Robin boundary conditions are not implemented.",err,error,*999)
            CASE(BOUNDARY_CONDITION_CAUCHY_ROW)
              !Set Cauchy boundary conditions
              CALL FlagError("Cauchy boundary conditions are not implemented.",err,error,*999)
            CASE(BOUNDARY_CONDITION_CONSTRAINED_ROW)
              !Set constrained row boundary conditions
              CALL FlagError("Constrained row boundary conditions are not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The RHS boundary condition of "// &
                & TRIM(NumberToVString(rhsBoundaryCondition,"*",err,error))// &
                & " for RHS variable dof number "// &
                & TRIM(NumberToVString(rhsVariableDOF,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsRowNumber
          
          IF(ASSOCIATED(rhsMapping)) THEN
            CALL FieldVariable_ParameterSetDataRestore(rhsVariable,FIELD_VALUES_SET_TYPE,rhsParameters,err,error,*999)
            IF(hasIntegratedValues) &
              & CALL FieldVariable_ParameterSetDataRestore(rhsVariable,FIELD_INTEGRATED_NEUMANN_SET_TYPE,rhsIntegratedParameters, &
              & err,error,*999)
          ENDIF
        ENDDO !equationsSetIdx
      
        !Add in any rows from any interface conditions
        DO interfaceConditionIdx=1,numberOfInterfaceConditions
          NULLIFY(interfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
          SELECT CASE(interfaceConditionMethod)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
            NULLIFY(interfaceEquations)
            CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
            NULLIFY(interfaceMapping)
            CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
            NULLIFY(interfaceMatrices)
            CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
            NULLIFY(interfaceRHSMapping)
            CALL InterfaceMapping_RHSMappingGet(interfaceMapping,interfaceRHSMapping,err,error,*999)
            NULLIFY(interfaceRHSVector)
            CALL InterfaceMatrices_RHSVectorGet(interfaceMatrices,interfaceRHSVector,err,error,*999)
            NULLIFY(interfaceColToSolverRowsMap)
            CALL SolverMapping_InterfaceColToSolverRowsMapGet(solverMapping,interfaceConditionIdx, &
              & interfaceColToSolverRowsMap,err,error,*999)
            NULLIFY(interfaceRHSDistributedVector)
            CALL InterfaceMatricesRHS_DistributedVectorGet(interfaceRHSVector,interfaceRHSDistributedVector,err,error,*999)
            !Worry about BCs on the Lagrange variables later.
            CALL DistributedVector_VectorCoupleAdd(solverRHSVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE, &
              & interfaceColToSolverRowsMap,1.0_DP,interfaceRHSDistributedVector,err,error,*999)
          CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The interface condition method of "// &
              & TRIM(NumberToVString(interfaceCondition%METHOD,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !interfaceConditionIdx
        !        
        !Start the update the solver RHS vector values
        CALL DistributedVector_UpdateStart(solverRHSVector,err,error,*999)
        CALL DistributedVector_UpdateFinish(solverRHSVector,err,error,*999)
        !Restore the solver RHS check data                 
        CALL DistributedVector_DataRestore(solverRHSVector,solverRHSCheckData,err,error,*999)            
      ENDIF !Update RHS

      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
          & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
        CALL Profiling_TimingsOutput(1,"Solver RHS assembly",userElapsed,systemElapsed,err,error,*999)
      ENDIF

    ENDIF !Calculate solver RHS
          
    !The solver matrices have only one residual vector
    NULLIFY(solverResidualVector)
    IF(selectionType==SOLVER_MATRICES_ALL.OR. &
      & selectionType==SOLVER_MATRICES_NONLINEAR_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RESIDUAL_ONLY.OR. &
      & selectionType==SOLVER_MATRICES_RHS_RESIDUAL_ONLY) THEN
      !Assemble residual vector
      !We assemble residual vector before RHS vector, then when assembling the RHS vector we subtract
      !the RHS terms for fixed BCs from the residual vector as this residual evaluation uses a matrix
      !vector product of the full equations matrix rather than the reduced solver matrix
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
      ENDIF
      
      CALL SolverMatrices_UpdateResidualGet(solverMatrices,updateResidual,err,error,*999)
      IF(updateResidual) THEN
        NULLIFY(solverResidualVector)
        CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,solverResidualVector,err,error,*999)
        !Initialise the residual to zero
        CALL DistributedVector_AllValuesSet(solverResidualVector,0.0_DP,err,error,*999)
        !Get the solver RHS check data                  
        NULLIFY(solverResidualCheckData)
        CALL DistributedVector_DataGet(solverResidualVector,solverResidualCheckData,err,error,*999)
        
        !Loop over the equations sets
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(equationsSetToSolverMatricesMap)
          CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
            & err,error,*999)
          NULLIFY(equationsRowToSolverRowsMap)
          CALL SolverMappingESToSMSMap_EquationsRowToSolverRowsMapGet(equationsSetToSolverMatricesMap, &
            & equationsRowToSolverRowsMap,err,error,*999)
          
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(lhsMapping)
          CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
          NULLIFY(equationsRowToLHSDOFMap)
          CALL EquationsMappingLHS_EquationsRowTOLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*999)
          NULLIFY(lhsDomainMapping)
          CALL FieldVariable_DomainMappingGet(lhsVariable,lhsDomainMapping,err,error,*999)
          CALL DomainMapping_BoundaryFinishGet(lhsDomainMapping,lhsBoundaryFinish,err,error,*999)
          CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*999)
          CALL EquationsMappingLHS_TotalNumberOfRowsGet(lhsMapping,totalNumberOfRows,err,error,*999)
          
          !Calculate the contributions from any linear matrices
          NULLIFY(linearMapping)
          NULLIFY(linearMatrices)
          NULLIFY(linearTempVector)
          CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
          IF(ASSOCIATED(linearMapping)) THEN
            NULLIFY(linearMatrices)
            CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
            CALL EquationsMatricesLinear_TempDistributedVectorGet(linearMatrices,linearTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(linearTempVector,0.0_DP,err,error,*999)                  
            NULLIFY(solverMatrixToEquationsMap)
            CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,1,solverMatrixToEquationsMap,err,error,*999)
            NULLIFY(solverMappingVariables)
            CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
            CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
            DO linearMatrixIdx=1,numberOfLinearMatrices
              NULLIFY(linearMatrix)
              CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,linearMatrixIdx,linearMatrix,err,error,*999)
              NULLIFY(linearDistributedMatrix)
              CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixIdx,linearVariable,err,error,*999)
              NULLIFY(solverMappingVariable)
              CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,linearVariable,solverMappingVariable, &
                & err,error,*999)
              IF(ASSOCIATED(solverMappingVariable)) THEN
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,matrixCoefficient,err,error,*999)                
                NULLIFY(dependentDistributedVector)
                CALL FieldVariable_ParameterSetVectorGet(linearVariable,FIELD_VALUES_SET_TYPE,dependentDistributedVector, &
                  & err,error,*999)
                CALL DistributedMatrix_MatrixByVectorAdd(linearDistributedMatrix,.FALSE.,dependentDistributedVector, &
                  & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,matrixCoefficient,linearTempVector,err,error,*999)
              ENDIF
            ENDDO !linearMatrixIdx
          ENDIF
          
          !Calculate the solver residual
          NULLIFY(nonlinearMapping)
          NULLIFY(nonlinearMatrices)
          NULLIFY(nonlinearTempVector)
          CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
          IF(ASSOCIATED(nonlinearMapping)) THEN
            NULLIFY(nonlinearMatrices)
            CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
            CALL EquationsMatricesNonlinear_TempDistributedVectorGet(nonlinearMatrices,nonlinearTempVector,err,error,*999)
            CALL DistributedVector_AllValuesSet(nonlinearTempVector,0.0_DP,err,error,*999)
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualMapping)
              CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
              !Check if any of the residual variables are part of the (1st) solver matrix (Jacobian)
              solverResidual=.FALSE.
              CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
              NULLIFY(solverMatrixToEquationsMap)
              CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,1,solverMatrixToEquationsMap,err,error,*999)
              NULLIFY(solverMappingVariables)
              CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
              DO residualVariableIdx=1,numberOfResidualVariables
                CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
                NULLIFY(solverMappingVariable)
                CALL SolverMappingVariables_VariableInListCheck(solverMappingVariables,residualVariable,solverMappingVariable, &
                  & err,error,*999)
                IF(ASSOCIATED(solverMappingVariable)) THEN
                  solverResidual=.TRUE.
                  EXIT
                ENDIF                
              ENDDO !residualVariableIdx
              IF(solverResidual) THEN
                !The residual has variables involved in the solver matrices so add it to the solver residual
                CALL EquationsMappingResidual_VectorCoefficientGet(residualMapping,residualCoefficient,err,error,*999)
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                NULLIFY(residualDistributedVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & residualDistributedVector,err,error,*999)
                CALL DistributedVector_VectorAdd(nonlinearTempVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
                  & residualCoefficient,residualDistributedVector,err,error,*999)
              ENDIF
            ENDDO !residualIdx
          ENDIF !nonlinear mapping

          !Couple the linear and nonlinear contributions to the solver residual
          IF(ASSOCIATED(linearTempVector)) CALL DistributedVector_VectorCoupleAdd(solverResidualVector, &
            & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,equationsRowToSolverRowsMap,1.0_DP,linearTempVector,err,error,*999)
          IF(ASSOCIATED(nonlinearTempVector)) CALL DistributedVector_VectorCoupleAdd(solverResidualVector, &
            & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,equationsRowToSolverRowsMap,1.0_DP,nonlinearTempVector,err,error,*999)
          
        ENDDO !equationsSetIdx
        
        !Loop over the interface conditions
        DO interfaceConditionIdx=1,numberOfInterfaceConditions
          NULLIFY(interfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
          NULLIFY(lagrangeField)
          CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
          NULLIFY(interfaceEquations)
          CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
          NULLIFY(interfaceMatrices)
          CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
          NULLIFY(interfaceMapping)
          CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
          CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
          CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)         
          SELECT CASE(interfaceConditionMethod)
          CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
            !Do nothing
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            penaltyMatrixIdx=numberOfInterfaceMatrices
            numberOfInterfaceMatrices=numberOfInterfaceMatrices-1
          CASE DEFAULT
            localError="The interface condition method of "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
              & " is invalid or not implemented."
            CALL FlagError(localError,err,error,*999)
          ENDSELECT
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          NULLIFY(lagrangeVector)
          CALL FieldVariable_ParameterSetVectorGet(lagrangeVariable,FIELD_VALUES_SET_TYPE,lagrangeVector,err,error,*999)
          !Calculate the contributions from any interface matrices
          DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
            !Calculate the interface matrix-Lagrange vector product residual contribution
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
            NULLIFY(lagrangeVariable)
            CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
            NULLIFY(interfaceDistributedMatrix)
            CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
            NULLIFY(interfaceTempVector)
            CALL InterfaceMatrix_TempDistributedVectorGet(interfaceMatrix,interfaceTempVector,err,error,*999)
           !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(interfaceTempVector,0.0_DP,err,error,*999)
            CALL DistributedMatrix_MatrixByVectorAdd(interfaceDistributedMatrix,.FALSE.,lagrangeVector, &
              & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP,interfaceTempVector,err,error,*999)
            NULLIFY(interfaceRowToSolverRowsMap)
            CALL SolverMapping_InterfaceRowToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceMatrixIdx, &
              & interfaceRowToSolverRowsMap,err,error,*999)
            !Add interface matrix residual contribution to the solver residual
            CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
              & interfaceRowToSolverRowsMap,1.0_DP,interfaceTempVector,err,error,*999)
            !Calculate the transposed interface matrix-dependent variable product residual contribution
            NULLIFY(dependentVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,dependentVariable,err,error,*999)
            NULLIFY(interfaceDistributedMatrix)
            CALL InterfaceMatrix_TransposeDistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
            NULLIFY(interfaceTempVector)
            CALL InterfaceMatrix_TempTransposeDistributedVectorGet(interfaceMatrix,interfaceTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(interfaceTempVector,0.0_DP,err,error,*999)
            NULLIFY(dependentVector)
            CALL FieldVariable_ParameterSetVectorGet(dependentVariable,FIELD_VALUES_SET_TYPE,dependentVector,err,error,*999)
            CALL DistributedMatrix_MatrixByVectorAdd(interfaceDistributedMatrix,.FALSE.,dependentVector, &
              & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP,interfaceTempVector,err,error,*999)
            NULLIFY(interfaceColToSolverRowsMap)
            CALL SolverMapping_InterfaceColToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceColToSolverRowsMap, &
              & err,error,*999)             
            !Add interface matrix residual contribution to the solver residual.
            CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
              & interfaceColToSolverRowsMap,1.0_DP,interfaceTempVector,err,error,*999)
          ENDDO !interfaceMatrixIdx
          SELECT CASE(interfaceConditionMethod)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            !Calculate the Lagrange-Lagrange vector product residual contribution from the penalty term
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,penaltyMatrixIdx,interfaceMatrix,err,error,*999)
            NULLIFY(lagrangeVariable)
            NULLIFY(interfaceDistributedMatrix)
            CALL InterfaceMatrix_DistributedMatrixGet(interfaceMatrix,interfaceDistributedMatrix,err,error,*999)
            NULLIFY(interfaceTempVector)
            CALL InterfaceMatrix_TempDistributedVectorGet(interfaceMatrix,interfaceTempVector,err,error,*999)
            !Initialise the linear temporary vector to zero
            CALL DistributedVector_AllValuesSet(interfaceTempVector,0.0_DP,err,error,*999)
            CALL DistributedMatrix_MatrixByVectorAdd(interfaceDistributedMatrix,.FALSE.,lagrangeVector, &
              & DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,1.0_DP,interfaceTempVector,err,error,*999)
            NULLIFY(interfaceRowToSolverRowsMap)
            CALL SolverMapping_InterfaceRowToSolverRowsMapGet(solverMapping,interfaceConditionIdx,interfaceMatrixIdx, &
              & interfaceRowToSolverRowsMap,err,error,*999)
            !Add interface matrix residual contribution to the solver residual
            CALL DistributedVector_VectorCoupleAdd(solverResidualVector,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE, &
              & interfaceRowToSolverRowsMap,1.0_DP,interfaceTempVector,err,error,*999)
          END SELECT
        ENDDO !interfaceConditionIdx
        !Start the update the solver residual vector values
        CALL DistributedVector_UpdateStart(solverResidualVector,err,error,*999)
      ENDIF !update residual
      IF(ASSOCIATED(solverResidualVector)) CALL DistributedVector_UpdateFinish(solverResidualVector,err,error,*999)
      
      IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
        CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
        CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
        userElapsed=userTime2(1)-userTime1(1)
        systemElapsed=systemTime2(1)-systemTime1(1)
        IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
          & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
        CALL Profiling_TimingsOutput(1,"Solver residual assembly",userElapsed,systemElapsed,err,error,*999)
      ENDIF
      
    ENDIF !Calculate residual

    !If required output the solver matrices
    IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
      & CALL SolverMatrices_Output(GENERAL_OUTPUT_TYPE,selectionType,solverMatrices,err,error,*999)

    EXITS("Solver_StaticAssemble")
    RETURN
999 IF(ALLOCATED(dependentParameters)) DEALLOCATE(dependentParameters)    
    ERRORSEXITS("Solver_StaticAssemble",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_StaticAssemble

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for a nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonAbsoluteToleranceSet
  SUBROUTINE Solver_QuasiNewtonAbsoluteToleranceSet(solver,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the absolute tolerance for
    REAL(DP), INTENT(IN) :: absoluteTolerance !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonAbsoluteToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(absoluteTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified absolute tolerance of "//TRIM(NumberToVString(absoluteTolerance,"*",err,error))// &
        & " is invalid. The absolute tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    quasiNewtonSolver%absoluteTolerance=absoluteTolerance
   
    EXITS("Solver_QuasiNewtonAbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonAbsoluteToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonAbsoluteToleranceSet

  !
  !================================================================================================================================
  !

  !>Enables/disables output monitoring for a nonlinear Quasi-Newton line search solver.
  SUBROUTINE Solver_QuasiNewtonLineSearchMonitorOutputSet(solver,linesearchMonitorOutputFlag,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the absolute tolerance for
    LOGICAL, INTENT(IN) :: linesearchMonitorOutputFlag !<Flag to determine whether to enable/disable linsearch monitor output.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    
    ENTERS("Solver_QuasiNewtonLineSearchMonitorOutputSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    CALL SolverNonlinearQuasiNewton_AssertIsLinesearch(quasiNewtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,lineSearchSolver,err,error,*999)
    linesearchSolver%linesearchMonitorOutput=linesearchMonitorOutputFlag
    
    EXITS("Solver_QuasiNewtonLineSearchMonitorOutputSet")
    RETURN
999 ERRORS("Solver_QuasiNewtonLineSearchMonitorOutputSet",err,error)    
    EXITS("Solver_QuasiNewtonLineSearchMonitorOutputSet")    
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonLineSearchMonitorOutputSet

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a Quasi-Newton solver 
  SUBROUTINE SolverNonlinearQuasiNewton_CreateFinish(quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the Quasi-Newton solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearQuasiNewton_CreateFinish",err,error,*999)

    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi-Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
      CALL SolverQuasiNewtonLinesearch_CreateFinish(quasiNewtonSolver%linesearchSolver,err,error,*999)
    CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
      CALL SolverNonlinearQuasiNewtonTrustregion_CreateFinish(quasiNewtonSolver%trustregionSolver,err,error,*999)
    CASE DEFAULT
      localError="The Quasi-Newton solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearQuasiNewton_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewton_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a Quasi-Newton solver and deallocate all memory
  RECURSIVE SUBROUTINE SolverNonlinear_QuasiNewtonFinalise(quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer the Quasi-Newton solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_QuasiNewtonFinalise",err,error,*999)

    IF(ASSOCIATED(quasiNewtonSolver)) THEN
      CALL SolverNonlinearQuasiNewton_LinesearchFinalise(quasiNewtonSolver%linesearchSolver,err,error,*999)
      CALL SolverNonlinearQuasiNewton_TrustregionFinalise(quasiNewtonSolver%trustregionSolver,err,error,*999)
      CALL Solver_Finalise(quasiNewtonSolver%linearSolver,err,error,*999)
      DEALLOCATE(quasiNewtonSolver)
    ENDIF
         
    EXITS("SolverNonlinear_QuasiNewtonFinalise")
    RETURN
999 ERRORSEXITS("SolverNonlinear_QuasiNewtonFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_QuasiNewtonFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Quasi-Newton solver for a nonlinear solver
  SUBROUTINE SolverNonlinear_QuasiNewtonInitialise(nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer the solver to initialise the Quasi-Newton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("SolverNonlinear_QuasiNewtonInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*998)
    IF(ASSOCIATED(nonlinearSolver%quasiNewtonSolver)) &
      & CALL FlagError("Quasi-Newton solver is already associated for this nonlinear solver.",err,error,*998)

    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    !Allocate and initialise a Quasi-Newton solver
    ALLOCATE(nonlinearSolver%quasiNewtonSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nonlinear solver Quasi-Newton solver.",err,error,*999)
    nonlinearSolver%quasiNewtonSolver%nonlinearSolver=>nonlinearSolver
    nonlinearSolver%quasiNewtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    nonlinearSolver%quasiNewtonSolver%totalNumberOfFunctionEvaluations=0
    nonlinearSolver%quasiNewtonSolver%totalNumberOfJacobianEvaluations=0
    nonlinearSolver%quasiNewtonSolver%maximumNumberOfIterations=50
    nonlinearSolver%quasiNewtonSolver%maximumNumberOfFunctionEvaluations=1000
    nonlinearSolver%quasiNewtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
    nonlinearSolver%quasiNewtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT
    nonlinearSolver%quasiNewtonSolver%absoluteTolerance=1.0E-10_DP
    nonlinearSolver%quasiNewtonSolver%relativeTolerance=1.0E-05_DP
    nonlinearSolver%quasiNewtonSolver%solutionTolerance=1.0E-05_DP
    NULLIFY(nonlinearSolver%quasiNewtonSolver%linesearchSolver)
    NULLIFY(nonlinearSolver%quasiNewtonSolver%trustregionSolver)
    NULLIFY(nonlinearSolver%quasiNewtonSolver%cellMLEvaluatorSolver)
    NULLIFY(nonlinearSolver%quasiNewtonSolver%convergenceTest)
    ALLOCATE(nonlinearSolver%quasiNewtonSolver%convergenceTest,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate convergence test object.",err,error,*999)
    nonlinearSolver%quasiNewtonSolver%convergenceTest%energyFirstIter = 0.0_DP
    nonlinearSolver%quasiNewtonSolver%convergenceTest%normalisedEnergy = 0.0_DP
    !Default to a Quasi-Newton linesearch solver
    nonlinearSolver%quasiNewtonSolver%quasiNewtonSolveType=SOLVER_QUASI_NEWTON_LINESEARCH
    CALL SolverNonlinearQuasiNewton_LinesearchInitialise(nonlinearSolver%quasiNewtonSolver,err,error,*999)
    !Default to a Quasi-Newton Good Broyden variant
    nonlinearSolver%quasiNewtonSolver%quasiNewtonType=SOLVER_QUASI_NEWTON_GOODBROYDEN
    nonlinearSolver%quasiNewtonSolver%restartType=SOLVER_QUASI_NEWTON_RESTART_PERIODIC
    nonlinearSolver%quasiNewtonSolver%restart=10
    nonlinearSolver%quasiNewtonSolver%scaleType=SOLVER_QUASI_NEWTON_SCALE_JACOBIAN
    !Create the linked linear solver
    ALLOCATE(nonlinearSolver%quasiNewtonSolver%linearSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Quasi-Newton solver linear solver.",err,error,*999)
    NULLIFY(nonlinearSolver%quasiNewtonSolver%linearSolver%solvers)
    CALL Solver_Initialise(nonlinearSolver%quasiNewtonSolver%linearSolver,err,error,*999)
    CALL Solver_LinearInitialise(nonlinearSolver%quasiNewtonSolver%linearSolver,err,error,*999)
    CALL Solver_LinkedSolverAdd(solver,nonlinearSolver%quasiNewtonSolver%linearSolver,SOLVER_LINEAR_TYPE,err,error,*999)
        
    EXITS("SolverNonlinear_QuasiNewtonInitialise")
    RETURN
999 CALL SolverNonlinear_QuasiNewtonFinalise(nonlinearSolver%quasiNewtonSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverNonlinear_QuasiNewtonInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_QuasiNewtonInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of Jacobian calculation type for a Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonJacobianCalculationSet
  SUBROUTINE Solver_QuasiNewtonJacobianCalculationTypeSet(solver,jacobianCalculationType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the Jacobian calculation type
    INTEGER(INTG), INTENT(IN) :: jacobianCalculationType !<The type of Jacobian calculation type to set \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonJacobianCalculationTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(jacobianCalculationType/=quasiNewtonSolver%jacobianCalculationType) THEN
      SELECT CASE(jacobianCalculationType)
      CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
        quasiNewtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED
      CASE(SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED)
        quasiNewtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED
      CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
        quasiNewtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
      CASE DEFAULT
        localError="The Jacobian calculation type of "// &
          & TRIM(NumberToVString(jacobianCalculationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_QuasiNewtonJacobianCalculationTypeSet")
    RETURN
999 ERRORS("Solver_QuasiNewtonJacobianCalculationTypeSet",err,error)
    EXITS("Solver_QuasiNewtonJacobianCalculationTypeSet")
    RETURN 1
    
  END SUBROUTINE Solver_QuasiNewtonJacobianCalculationTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a Quasi-Newton solver.
  SUBROUTINE SolverNonlinearQuasiNewton_LibraryTypeSet(quasiNewtonSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer the Quasi-Newton solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the Quasi-Newton solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearQuasiNewton_LibraryTypeSet",err,error,*999)

    
    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi-Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
      NULLIFY(linesearchSolver)
      CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        linesearchSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Quasi-Newton linesearch solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
      NULLIFY(trustregionSolver)
      CALL SolverNonlinearQuasiNewton_TrustregionSolverGet(quasiNewtonSolver,trustregionSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        trustregionSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Quasi-Newton trustregion solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The Quasi-Newton solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearQuasiNewton_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewton_LibraryTypeSet

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the convergence test for a Quasi-Newton nonlinear solver \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonConvergenceTestSet
  SUBROUTINE Solver_QuasiNewtonConvergenceTestTypeSet(solver,convergenceTestType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the convergence test for
    INTEGER(INTG), INTENT(IN) :: convergenceTestType !<The convergence test type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver 
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonConvergenceTestTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    SELECT CASE(convergenceTestType)
    CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
      quasiNewtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT
    CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
      quasiNewtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM
    CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
      quasiNewtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO
    CASE DEFAULT
      localError="The specified convergence test type of "//TRIM(NumberToVString(convergenceTestType, &
        & "*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_QuasiNewtonConvergenceTestTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonConvergenceTestTypeSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonConvergenceTestTypeSet
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Quasi-Newton line search solver
  SUBROUTINE SolverQuasiNewtonLinesearch_CreateFinish(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer the nonlinear Quasi-Newton line search solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx,groupCommunicator, &
      & numberOfEquationsSets,numberOfInterfaceConditions,numberOfInterfaceMatrices,numberOfLinearMatrices,sparsityType, &
      & symmetryType
    EXTERNAL :: Problem_SolverJacobianEvaluatePetsc
    EXTERNAL :: Problem_SolverJacobianFDCalculatePetsc
    EXTERNAL :: Problem_SolverResidualEvaluatePetsc
    EXTERNAL :: Problem_SolverConvergenceTestPetsc
    EXTERNAL :: Problem_SolverNonlinearMonitorPETSC
    TYPE(DistributedMatrixType), POINTER :: jacobianMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: jacobianPETScMatrix
    TYPE(DistributedVectorType), POINTER :: matrixTempVector,residualVector
    TYPE(DistributedVectorPETScType), POINTER :: residualPETScVector
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: linearVariable,interfaceVariable,lagrangeVariable,lhsVariable
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(SolverType), POINTER :: linkedLinearSolver,solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverJacobian
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("SolverQuasiNewtonLinesearch_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)

    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet(linesearchSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearQuasiNewton_NonlinearSolverGet(quasiNewtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    SELECT CASE(linesearchSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Loop over the equations set in the solver equations
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        NULLIFY(lhsVariable)
        CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lhsVariable,domainMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          !If there are any linear matrices create temporary vector for matrix-vector products
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
          DO equationsMatrixIdx=1,numberOfLinearMatrices
            NULLIFY(equationsMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
            NULLIFY(matrixTempVector)
            CALL EquationsMatrix_TempDistributedVectorExists(equationsMatrix,matrixTempVector,err,error,*999)
            IF(.NOT.ASSOCIATED(matrixTempVector)) THEN
              CALL DistributedVector_CreateStart(domainMapping,equationsMatrix%tempVector,err,error,*999)
              CALL DistributedVector_DataTypeSet(equationsMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
              CALL DistributedVector_CreateFinish(equationsMatrix%tempVector,err,error,*999)
            ENDIF
          ENDDO !equationsMatrixIdx
        ENDIF
      ENDDO !equationsSetIdx
      !Loop over the interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMatrices)
        CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
        !Create temporary vector for matrix-vector products
        CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
        DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
          IF(.NOT.ASSOCIATED(interfaceMatrix%tempVector)) THEN
            NULLIFY(interfaceVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,interfaceVariable,err,error,*999)
            NULLIFY(domainMapping)
            CALL FieldVariable_DomainMappingGet(interfaceVariable,domainMapping,err,error,*999)
            !Set up the temporary interface distributed vector to be used with interface matrices
            CALL DistributedVector_CreateStart(domainMapping,interfaceMatrix%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(interfaceMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(interfaceMatrix%tempVector,err,error,*999)
            !Set up the temporary interface distributed vector to be used with transposed interface matrices
            NULLIFY(domainMapping)
            CALL FieldVariable_DomainMappingGet(lagrangeVariable,domainMapping,err,error,*999)
            CALL DistributedVector_CreateStart(domainMapping,interfaceMatrix%tempTransposeVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(interfaceMatrix%tempTransposeVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
              & err,error,*999)
            CALL DistributedVector_CreateFinish(interfaceMatrix%tempTransposeVector,err,error,*999)
          ENDIF
        ENDDO !interfaceMatrixIdx
      ENDDO !interfaceConiditionIdx
      !Create the PETSc SNES solver
      CALL PETSc_SNESCreate(groupCommunicator,linesearchSolver%snes,err,error,*999)
      !Set the nonlinear solver type to be a Quasi-Newton line search solver
      CALL PETSc_SNESSetType(linesearchSolver%snes,PETSC_SNESQN,err,error,*999)
      !Following routines don't work for petsc version < 3.5.
      !Set the nonlinear Quasi-Newton type
      SELECT CASE(quasiNewtonSolver%quasiNewtonType)
      CASE(SOLVER_QUASI_NEWTON_LBFGS)
        CALL PETSc_SNESQNSetType(linesearchSolver%snes,PETSC_SNES_QN_LBFGS,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_GOODBROYDEN)
        CALL PETSc_SNESQNSetType(linesearchSolver%snes,PETSC_SNES_QN_BROYDEN,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_BADBROYDEN)
        CALL PETSc_SNESQNSetType(linesearchSolver%snes,PETSC_SNES_QN_BADBROYDEN,err,error,*999)
      CASE DEFAULT
        localError="The specified nonlinear Quasi-Newton type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the nonlinear Quasi-Newton restart type
      SELECT CASE(quasiNewtonSolver%restartType)
      CASE(SOLVER_QUASI_NEWTON_RESTART_NONE)
        CALL PETSc_SNESQNSetRestartType(linesearchSolver%snes,PETSC_SNES_QN_RESTART_NONE,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_RESTART_POWELL)
        CALL PETSc_SNESQNSetRestartType(linesearchSolver%snes,PETSC_SNES_QN_RESTART_POWELL,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_RESTART_PERIODIC)
        CALL PETSc_SNESQNSetRestartType(linesearchSolver%snes,PETSC_SNES_QN_RESTART_PERIODIC,err,error,*999)
      CASE DEFAULT
        localError="The specified nonlinear Quasi-Newton restart type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%restartType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the nonlinear Quasi-Newton scale type
      SELECT CASE(quasiNewtonSolver%scaleType)
      CASE(SOLVER_QUASI_NEWTON_SCALE_NONE)
        CALL PETSc_SNESQNSetScaleType(linesearchSolver%snes,PETSC_SNES_QN_SCALE_NONE,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_SCALE_SHANNO)
        CALL PETSc_SNESQNSetScaleType(linesearchSolver%snes,PETSC_SNES_QN_SCALE_SHANNO,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_SCALE_LINESEARCH)
        CALL PETSc_SNESQNSetScaleType(linesearchSolver%snes,PETSC_SNES_QN_SCALE_LINESEARCH,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_SCALE_JACOBIAN)
        CALL PETSc_SNESQNSetScaleType(linesearchSolver%snes,PETSC_SNES_QN_SCALE_JACOBIAN,err,error,*999)
      CASE DEFAULT
        localError="The specified nonlinear Quasi-Newton scale type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%scaleType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the Quasi-Newton restart
      !Not implemented yet, as there is currently no routine in PETSc for this. If need be, this can be set in your petscrc file.
      !Create the solver matrices and vectors
      NULLIFY(linkedLinearSolver)
      CALL SolverNonlinearQuasiNewton_LinkedLinearSolverGet(quasiNewtonSolver,linkedLinearSolver,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
      CALL SolverMatrices_LibraryTypeSet(solverMatrices,SOLVER_PETSC_LIBRARY,err,error,*999)
      CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
      SELECT CASE(sparsityType)
      CASE(SOLVER_SPARSE_MATRICES)
        CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
      CASE(SOLVER_FULL_MATRICES)
        CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
      SELECT CASE(symmetryType)
      CASE(SOLVER_SYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
      CASE(SOLVER_UNSYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      !Link linear solver
      linkedLinearSolver%solverEquations=>solver%solverEquations
      !Finish the creation of the linear solver
      CALL SolverLinear_CreateFinish(linkedLinearSolver%linearSolver,err,error,*999)
      !Associate linear solver's KSP to nonlinear solver's SNES
      NULLIFY(linearSolver)
      CALL Solver_LinearSolverGet(linkedLinearSolver,linearSolver,err,error,*999)
      SELECT CASE(linearSolver%linearSolveType)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        NULLIFY(directSolver)
        CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
        CALL PETSc_SNESSetKsp(linesearchSolver%snes,directSolver%ksp,err,error,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        NULLIFY(iterativeSolver)
        CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
        CALL PETSc_SNESSetKsp(linesearchSolver%snes,iterativeSolver%ksp,err,error,*999)
      CASE DEFAULT
        localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the nonlinear function
      NULLIFY(residualVector)
      CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualVector,err,error,*999)
      NULLIFY(residualPETScVector)
      CALL DistributedVector_PETScVectorGet(residualVector,residualPETScVector,err,error,*999)
      !Pass the linesearch solver object rather than the temporary solver
      CALL PETSc_SNESSetFunction(linesearchSolver%snes,residualPETScVector%vector,Problem_SolverResidualEvaluatePetsc, &
        & linesearchSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
      SELECT CASE(linesearchSolver%quasiNewtonSolver%convergenceTestType)
      CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
        !Default convergence test, do nothing
      CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM,SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
        CALL PETSc_SNESSetConvergenceTest(linesearchSolver%snes,Problem_SolverConvergenceTestPetsc, &
          & linesearchSolver%quasiNewtonSolver%nonlinearSolver%SOLVER,err,error,*999)
      CASE DEFAULT
        localError="The specified convergence test type of "//TRIM(NumberToVString(linesearchSolver% &
          & quasiNewtonSolver%convergenceTestType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT                 
      !Set the Jacobian
      IF(solverMatrices%numberOfMatrices/=1) THEN
        localError="Invalid number of solver matrices. The number of solver matrices is "// &
          & TRIM(NumberToVString(solverMatrices%numberOfMatrices,"*",err,error))//" and it should be 1."     
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(solverJacobian)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverJacobian,err,error,*999)
      NULLIFY(jacobianMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverJacobian,jacobianMatrix,err,error,*999)
      NULLIFY(jacobianPETScMatrix)
      CALL DistributedMatrix_PETScMatrixGet(jacobianMatrix,jacobianPETScMatrix,err,error,*999)
      SELECT CASE(quasiNewtonSolver%jacobianCalculationType)
      CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
        CALL FlagError("Cannot have no Jacobian calculation for a PETSc nonlinear linesearch solver.",err,error,*999)
      CASE(SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED)
        solverJacobian%updateMatrix=.TRUE. !CMISS will fill in the Jacobian values
        !Pass the linesearch solver object rather than the temporary solver
        CALL PETSc_SNESSetJacobian(linesearchSolver%snes,jacobianPETScMatrix%matrix,jacobianPETScMatrix%matrix, &
          & Problem_SolverJacobianEvaluatePetsc,linesearchSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
      CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
        solverJacobian%updateMatrix=.FALSE. !Petsc will fill in the Jacobian values
        CALL DistributedMatrix_Form(jacobianMatrix,err,error,*999)
        CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
        SELECT CASE(sparsityType)
        CASE(SOLVER_SPARSE_MATRICES)
          CALL PETSc_MatColoringCreate(jacobianPETScMatrix%matrix,linesearchSolver%jacobianMatColoring,err,error,*999)
          CALL PETSc_MatColoringSetType(linesearchSolver%jacobianMatColoring,PETSC_MATCOLORING_SL,err,error,*999)
          CALL PETSc_MatColoringSetFromOptions(linesearchSolver%jacobianMatColoring,err,error,*999)
          CALL PETSc_MatColoringApply(linesearchSolver%jacobianMatColoring,linesearchSolver%jacobianISColoring,err,error,*999)
          CALL PETSc_MatColoringDestroy(linesearchSolver%jacobianMatColoring,err,error,*999)
          CALL PETSc_MatFDColoringCreate(jacobianPETScMatrix%matrix,linesearchSolver%jacobianISColoring, &
            & linesearchSolver%jacobianMatFDColoring,err,error,*999)
          CALL PETSc_MatFDColoringSetFunction(linesearchSolver%jacobianMatFDColoring,Problem_SolverResidualEvaluatePetsc, &
            & linesearchSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
          CALL PETSc_MatFDColoringSetFromOptions(linesearchSolver%jacobianMatFDColoring,err,error,*999)
          CALL PETSc_MatFDColoringSetup(jacobianPETScMatrix%matrix,linesearchSolver%jacobianISColoring, &
            & linesearchSolver%jacobianMatFDColoring,err,error,*999)
          CALL PETSc_ISColoringDestroy(linesearchSolver%jacobianISColoring,err,error,*999)
        CASE(SOLVER_FULL_MATRICES)
          !Do nothing
        CASE DEFAULT
          localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL PETSc_SNESSetJacobian(linesearchSolver%snes,jacobianPETScMatrix%matrix,jacobianPETScMatrix%matrix, &
          & Problem_SolverJacobianFDCalculatePetsc,linesearchSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
      CASE DEFAULT
        localError="The Jacobian calculation type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%jacobianCalculationType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(solver%outputType>=SOLVER_MONITOR_OUTPUT) THEN
        !Set the monitor
        !Pass the linesearch solver object rather than the temporary solver
        CALL PETSc_SNESMonitorSet(linesearchSolver%snes,Problem_SolverNonlinearMonitorPETSC, &
          & linesearchSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
      ENDIF
      CALL PETSc_SNESGetLineSearch(linesearchSolver%snes,linesearchSolver%snesLineSearch,err,error,*999)
      !Set the line search type and order where applicable
      SELECT CASE(linesearchSolver%linesearchType)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH_BASIC)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_BASIC,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH_L2)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_L2,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH_CP)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_CP,err,error,*999)
      CASE DEFAULT
        localError="The nonlinear Quasi-Newton line search type of "// &
          & TRIM(NumberToVString(linesearchSolver%linesearchType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set step tolerances, leave iterative line search options as defaults
      CALL PETSc_SNESLineSearchSetTolerances(linesearchSolver%snesLineSearch,linesearchSolver%linesearchStepTolerance, &
        & linesearchSolver%linesearchMaxstep,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER, &
        & err,error,*999)
      IF(linesearchSolver%linesearchMonitorOutput) THEN
        CALL PETSc_SNESLineSearchSetMonitor(linesearchSolver%snesLineSearch,PETSC_TRUE,err,error,*999)
      ELSE
        CALL PETSc_SNESLineSearchSetMonitor(linesearchSolver%snesLineSearch,PETSC_FALSE,err,error,*999)
      ENDIF
      !Set the tolerances for the SNES solver
      CALL PETSc_SNESSetTolerances(linesearchSolver%snes,quasiNewtonSolver%absoluteTolerance, &
        & quasiNewtonSolver%relativeTolerance,quasiNewtonSolver%solutionTolerance, &
        & quasiNewtonSolver%maximumNumberOfIterations,quasiNewtonSolver%maximumNumberOfFunctionEvaluations,err,error,*999)
      !Set any further SNES options from the command line options
      CALL PETSc_SNESSetFromOptions(linesearchSolver%snes,err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "// &
        & TRIM(NumberToVString(linesearchSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverQuasiNewtonLinesearch_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverQuasiNewtonLinesearch_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverQuasiNewtonLinesearch_CreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise a nonlinear Quasi-Newton line search solver and deallocate all memory
  SUBROUTINE SolverNonlinearQuasiNewton_LinesearchFinalise(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer the nonlinear Quasi-Newton line search solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("SolverNonlinearQuasiNewton_LinesearchFinalise",err,error,*999)

    IF(ASSOCIATED(linesearchSolver)) THEN
      CALL PETSc_ISColoringFinalise(linesearchSolver%jacobianISColoring,err,error,*999)
      CALL PETSc_MatColoringFinalise(linesearchSolver%jacobianMatColoring,err,error,*999)
      CALL PETSc_MatFDColoringFinalise(linesearchSolver%jacobianMatFDColoring,err,error,*999)
      CALL PETSc_SNESLineSearchFinalise(linesearchSolver%snesLineSearch,err,error,*999)
      CALL PETSc_SNESFinalise(linesearchSolver%snes,err,error,*999)
      DEALLOCATE(linesearchSolver)
    ENDIF
        
    EXITS("SolverNonlinearQuasiNewton_LinesearchFinalise")
    RETURN
999 ERRORS("SolverNonlinearQuasiNewton_LinesearchFinalise",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinesearchFinalise")
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_LinesearchFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a nonlinear Quasi-Newton line search solver for a Quasi-Newton solver
  SUBROUTINE SolverNonlinearQuasiNewton_LinesearchInitialise(quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer the nonlinear Quasi-Newton solver to initialise the Quasi-Newton line search solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
  
    ENTERS("SolverNonlinearQuasiNewton_LinesearchInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi-Newton solver is not associated.",err,error,*998)
    IF(ASSOCIATED(quasiNewtonSolver%linesearchSolver)) &
      & CALL FlagError("Quasi-Newton line search solver is already associated for this Quasi-Newton solver.",err,error,*998)
      
    !Allocate and initialise the Quasi-Newton linesearch solver
    ALLOCATE(quasiNewtonSolver%linesearchSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nonlinear solver Quasi-Newton line search solver.",err,error,*999)
    quasiNewtonSolver%linesearchSolver%quasiNewtonSolver=>quasiNewtonSolver
    quasiNewtonSolver%linesearchSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    quasiNewtonSolver%linesearchSolver%linesearchType=SOLVER_QUASI_NEWTON_LINESEARCH_CP
    quasiNewtonSolver%linesearchSolver%linesearchMaxstep=1.0E8_DP
    quasiNewtonSolver%linesearchSolver%linesearchStepTolerance=CONVERGENCE_TOLERANCE
    CALL PETSc_MatColoringInitialise(quasiNewtonSolver%linesearchSolver%jacobianMatColoring,err,error,*999)
    CALL PETSc_ISColoringInitialise(quasiNewtonSolver%linesearchSolver%jacobianISColoring,err,error,*999)
    CALL PETSc_MatFDColoringInitialise(quasiNewtonSolver%linesearchSolver%jacobianMatFDColoring,err,error,*999)
    CALL PETSc_SNESInitialise(quasiNewtonSolver%linesearchSolver%snes,err,error,*999)
    CALL PETSc_SNESLinesearchInitialise(quasiNewtonSolver%linesearchSolver%snesLineSearch,err,error,*999)
    quasiNewtonSolver%linesearchSolver%linesearchMonitorOutput=.FALSE.
        
    EXITS("SolverNonlinearQuasiNewton_LinesearchInitialise")
    RETURN
999 CALL SolverNonlinearQuasiNewton_LinesearchFinalise(quasiNewtonSolver%linesearchSolver,dummyErr,dummyError,*998)
998 ERRORS("SolverNonlinearQuasiNewton_LinesearchInitialise",err,error)
    EXITS("SolverNonlinearQuasiNewton_LinesearchInitialise")
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewton_LinesearchInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the line search maximum step for a nonlinear Quasi-Newton linesearch solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonLineSearchMaxStepSet
  SUBROUTINE Solver_QuasiNewtonLinesearchMaxStepSet(solver,linesearchMaximumStep,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search maximum step for
    REAL(DP), INTENT(IN) :: linesearchMaximumStep !<The line search maximum step to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonLinesearchMaxStepSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    CALL SolverNonlinearQuasiNewton_AssertIsLinesearch(quasiNewtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*999)
    IF(linesearchMaximumStep<=ZERO_TOLERANCE) THEN
      localError="The specified line search maximum step of "// &
        & TRIM(NumberToVString(linesearchMaximumStep,"*",err,error))// &
        & " is invalid. The line search maximum step must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    linesearchSolver%linesearchMaxstep=linesearchMaximumStep
    
    EXITS("Solver_QuasiNewtonLinesearchMaxStepSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonLinesearchMaxStepSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonLinesearchMaxStepSet
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Quasi-Newton line search solver 
  SUBROUTINE SolverNonlinearQuasiNewtonLinesearch_Solve(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the nonlinear Quasi-Newton line search solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    !EXTERNAL :: Problem_SolverResidualEvaluatePetsc
    INTEGER(INTG) :: convergedReason,numberOfIterations,numberOfSolverMatrices
    REAL(DP) :: functionNorm
    TYPE(DistributedVectorType), POINTER :: rhsVector,solverVector
    TYPE(DistributedVectorPETScType), POINTER :: rhsPETScVector,solverPETScVector
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(PetscVecType) :: functionPETScVector
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearQuasiNewtonLinesearch_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)

    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinearQuasiNewtonLinesearch_QuasiNewtonSolverGet(linesearchSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearQuasiNewton_NonlinearSolverGet(quasiNewtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    IF(numberOfSolverMatrices/=1) THEN
      localError="The number of solver matrices of "//TRIM(NumberToVString(numberOfSolverMatrices,"*",err,error))// &
        & " is invalid. There should only be one solver matrix for a Quasi-Newton linesearch solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(solverMatrix)
    CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
    NULLIFY(rhsVector)
    CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,rhsVector,err,error,*999)
    NULLIFY(solverVector)
    CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
    SELECT CASE(linesearchSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      SELECT CASE(quasiNewtonSolver%solutionInitialiseType)
      CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
        !Zero the solution vector
        CALL DistributedVector_AllValuesSet(solverVector,0.0_DP,err,error,*999)
      CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
        !Make sure the solver vector contains the current dependent field values
        CALL Solver_SolutionUpdate(solver,err,error,*999)
      CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
        !Do nothing
      CASE DEFAULT
        localError="The Quasi-Newton solver solution initialise type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%solutionInitialiseType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      NULLIFY(rhsPETScVector)
      CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
      NULLIFY(solverPETScVector)
      CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
      !Set step tolerances, leave iterative line search options as defaults in case the user has changed them
      CALL PETSc_SNESLineSearchSetTolerances(linesearchSolver%snesLineSearch,linesearchSolver%linesearchStepTolerance, &
        & linesearchSolver%linesearchMaxstep,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER, &
        & err,error,*999)
      !Set the tolerances for the SNES solver in case the user has changed them
      CALL PETSc_SNESSetTolerances(linesearchSolver%snes,quasiNewtonSolver%absoluteTolerance, &
        & quasiNewtonSolver%relativeTolerance,quasiNewtonSolver%solutionTolerance, &
        & quasiNewtonSolver%maximumNumberOfIterations,quasiNewtonSolver%maximumNumberOfFunctionEvaluations,err,error,*999)
      !Solve the nonlinear equations     
      CALL PETSc_SNESSolve(linesearchSolver%snes,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999)
      !Check for convergence
      CALL PETSc_SNESGetConvergedReason(linesearchSolver%snes,convergedReason,err,error,*999)
      SELECT CASE(convergedReason)
      CASE(PETSC_SNES_DIVERGED_FUNCTION_COUNT)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged function count.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LINEAR_SOLVE)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged linear solve.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_FNORM_NAN)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged F Norm NaN.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_MAX_IT)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged maximum iterations.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LINE_SEARCH)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged line search.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LOCAL_MIN)
        CALL FlagWarning("Nonlinear line search solver did not converge. PETSc diverged local minimum.",err,error,*999)
      END SELECT
      IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
        !Output solution characteristics
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Quasi-Newton linesearch solver parameters:",err,error,*999)
        CALL PETSc_SNESGetIterationNumber(linesearchSolver%snes,numberOfIterations,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final number of iterations = ",numberOfIterations,err,error,*999)
        CALL PETSc_SNESGetFunction(linesearchSolver%snes,functionPETScVector,err,error,*999)
        CALL PETSc_VecNorm(functionPETScVector,PETSC_NORM_2,functionNorm,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final function norm = ",functionNorm,err,error,*999)
        SELECT CASE(convergedReason)
        CASE(PETSC_SNES_CONVERGED_FNORM_ABS)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm absolute.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_FNORM_RELATIVE)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm relative.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_ITS)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged its.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_ITERATING)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged iterating.",err,error,*999)
        END SELECT
      ENDIF
    CASE DEFAULT
      localError="The Quasi-Newton line search solver library type of "// &
        & TRIM(NumberToVString(linesearchSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinearQuasiNewtonLinesearch_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewtonLinesearch_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewtonLinesearch_Solve
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search step tolerance for a nonlinear Quasi-Newton line search solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonLineSearchStepTolSet
  SUBROUTINE Solver_QuasiNewtonLinesearchStepToleranceSet(solver,linesearchStepTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search step tolerance for
    REAL(DP), INTENT(IN) :: linesearchStepTolerance !<The line search step tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonLinesearchStepToleranceSet",err,error,*999)
    
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*999)
    IF(linesearchStepTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified line search step tolerance of "// &
        & TRIM(NumberToVString(linesearchStepTolerance,"*",err,error))// &
        & " is invalid. The line search step tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    linesearchSolver%linesearchStepTolerance=linesearchStepTolerance
    
    EXITS("Solver_QuasiNewtonLinesearchStepToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonLinesearchStepToleranceSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_QuasiNewtonLinesearchStepToleranceSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search type for a nonlinear Quasi-Newton linesearch solver \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonLineSearchTypeSet
  SUBROUTINE Solver_QuasiNewtonLinesearchTypeSet(solver,linesearchType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search type for
    INTEGER(INTG), INTENT(IN) :: linesearchType !<The line search type to set \see SolverRoutines_QuasiNewtonLineSearchTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonLinesearchTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,linesearchSolver,err,error,*999)
    SELECT CASE(linesearchType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH_BASIC)
      linesearchSolver%linesearchType=SOLVER_QUASI_NEWTON_LINESEARCH_BASIC
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH_L2)
      linesearchSolver%linesearchType=SOLVER_QUASI_NEWTON_LINESEARCH_L2
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH_CP)
      linesearchSolver%linesearchType=SOLVER_QUASI_NEWTON_LINESEARCH_CP
    CASE DEFAULT
      localError="The specified line search type of "//TRIM(NumberToVString(linesearchType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_QuasiNewtonLinesearchTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonLinesearchTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_QuasiNewtonLinesearchTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of function evaluations for a nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonMaximumFunctionEvaluationsSet
  SUBROUTINE Solver_QuasiNewtonMaximumFunctionEvaluationsSet(solver,maximumFunctionEvaluations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the maximum function evaluations for
    INTEGER(INTG), INTENT(IN) :: maximumFunctionEvaluations !<The maximum function evaluations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonMaximumFunctionEvaluationsSet",err,error,*999)
    
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(maximumFunctionEvaluations<1) THEN
      localError="The specified maximum number of function evaluations of "// &
        & TRIM(NumberToVString(maximumFunctionEvaluations,"*",err,error))// &
        & " is invalid. The maximum number of function evaluations must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    quasiNewtonSolver%maximumNumberOfFunctionEvaluations=maximumFunctionEvaluations
    
    EXITS("Solver_QuasiNewtonMaximumFunctionEvaluationsSet")
    RETURN
999 ERRORS("Solver_QuasiNewtonMaximumFunctionEvaluationsSet",err,error)
    EXITS("Solver_QuasiNewtonMaximumFunctionEvaluationsSet")
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonMaximumFunctionEvaluationsSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of iterations for a nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonMaximumIterationsSet
  SUBROUTINE Solver_QuasiNewtonMaxNumberOfIterationsSet(solver,maxNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: maxNumberOfIterations !<The maximum iterations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonMaxNumberOfIterationsSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(maxNumberOfIterations<1) THEN
      localError="The specified maximum iterations of "//TRIM(NumberToVString(maxNumberOfIterations,"*",err,error))// &
        & " is invalid. The maximum number of iterations must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    quasiNewtonSolver%maximumNumberOfIterations=maxNumberOfIterations
    
    EXITS("Solver_QuasiNewtonMaxNumberOfIterationsSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonMaxNumberOfIterationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_QuasiNewtonMaxNumberOfIterationsSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the relative tolerance for a nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonRelativeToleranceSet
  SUBROUTINE Solver_QuasiNewtonRelativeToleranceSet(solver,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the relative tolerance for
    REAL(DP), INTENT(IN) :: relativeTolerance !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonRelativeToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(relativeTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified relative tolerance of "//TRIM(NumberToVString(relativeTolerance,"*",err,error))// &
        & " is invalid. The relative tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    quasiNewtonSolver%relativeTolerance=relativeTolerance
    
    EXITS("Solver_QuasiNewtonRelativeToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonRelativeToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonRelativeToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution initialisation for a nonlinear Quasi-Newton solver
  SUBROUTINE Solver_QuasiNewtonSolutionInitialiseTypeSet(solver,solutionInitialiseType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the solution tolerance for
    INTEGER(INTG), INTENT(IN) :: solutionInitialiseType !<The solution initialise type to set \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonSolutionInitialiseTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    SELECT CASE(solutionInitialiseType)
    CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
      quasiNewtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_ZERO
    CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
      quasiNewtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
      quasiNewtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_NO_CHANGE
    CASE DEFAULT
      localError="The specified solution initialise type  of "// &
        & TRIM(NumberToVString(solutionInitialiseType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_QuasiNewtonSolutionInitialiseTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonSolutionInitialiseTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonSolutionInitialiseTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution tolerance for a nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonSolutionToleranceSet
  SUBROUTINE Solver_QuasiNewtonSolutionToleranceSet(solver,solutionTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the solution tolerance for
    REAL(DP), INTENT(IN) :: solutionTolerance !<The solution tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonSolutionToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(solutionTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified solution tolerance of "//TRIM(NumberToVString(solutionTolerance,"*",err,error))// &
        & " is invalid. The relative tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    quasiNewtonSolver%solutionTolerance=solutionTolerance
    
    EXITS("Solver_QuasiNewtonSolutionToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonSolutionToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonSolutionToleranceSet
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Quasi-Newton solver 
  SUBROUTINE SolverNonlinearQuasiNewton_Solve(quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer to the nonlinear Quasi-Newton solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearQuasiNewton_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi-Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
      CALL SolverNonlinearQuasiNewtonLinesearch_Solve(quasiNewtonSolver%linesearchSolver,err,error,*999)
    CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
      CALL SolverNonlinearQuasiNewtonTrustregion_Solve(quasiNewtonSolver%trustregionSolver,err,error,*999)
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearQuasiNewton_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_Solve
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Quasi-Newton trust region solver
  SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_CreateFinish(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer the nonlinear Quasi-Newton trust region solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx,groupCommunicator,numberOfEquationsSets,numberOfLinearMatrices,symmetryType
    EXTERNAL :: Problem_SolverResidualEvaluatePetsc
    TYPE(DistributedVectorType), POINTER :: residualVector
    TYPE(DistributedVectorPETScType), POINTER :: residualPETScVector
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: lhsVariable,linearVariable
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("SolverNonlinearQuasiNewtonTrustregion_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trust region solver is not associated.",err,error,*999)

    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet(trustregionSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearQuasiNewton_NonlinearSolverGet(quasiNewtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    SELECT CASE(trustregionSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Loop over the equations set in the solver equations
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        NULLIFY(lhsVariable)
        CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lhsVariable,domainMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          !If there are any linear matrices create temporary vector for matrix-vector products
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
          DO equationsMatrixIdx=1,numberOfLinearMatrices
            NULLIFY(equationsMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
            IF(.NOT.ASSOCIATED(equationsMatrix%tempVector)) THEN
              CALL DistributedVector_CreateStart(domainMapping,equationsMatrix%tempVector,err,error,*999)
              CALL DistributedVector_DataTypeSet(equationsMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
              CALL DistributedVector_CreateFinish(equationsMatrix%tempVector,err,error,*999)
            ENDIF
          ENDDO !equationsMatrixIdx
        ENDIF
      ENDDO !equationsSetIdx                  
      !Create the solver matrices and vectors
      NULLIFY(solverMatrices)
      CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
      CALL SolverMatrices_LibraryTypeSet(solverMatrices,SOLVER_PETSC_LIBRARY,err,error,*999)
!!TODO: set up the matrix structure if using an analytic Jacobian
      CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
      SELECT CASE(symmetryType)
      CASE(SOLVER_SYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
      CASE(SOLVER_UNSYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      !Create the PETSc SNES solver
      CALL PETSc_SNESCreate(groupCommunicator,trustregionSolver%snes,err,error,*999)
      !Set the nonlinear solver type to be a Quasi-Newton trust region solver
      CALL PETSc_SNESSetType(trustregionSolver%snes,PETSC_SNESNEWTONTR,err,error,*999)
      !Set the nonlinear function
      NULLIFY(residualVector)
      CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualVector,err,error,*999)
      NULLIFY(residualPETScVector)
      CALL DistributedVector_PETScVectorGet(residualVector,residualPETScVector,err,error,*999)
      CALL PETSc_SNESSetFunction(trustregionSolver%snes,residualPETScVector%vector, &
        & Problem_SolverResidualEvaluatePetsc,trustregionSolver%quasiNewtonSolver%nonlinearSolver%solver,err,error,*999)
      !Set the Jacobian if necessary
      !Set the trust region delta ???
                  
      !Set the trust region tolerance
      CALL PETSc_SNESSetTrustRegionTolerance(trustregionSolver%snes,trustregionSolver%trustregionTolerance, &
        & err,error,*999)
      !Set the tolerances for the SNES solver
      CALL PETSc_SNESSetTolerances(trustregionSolver%snes,quasiNewtonSolver%absoluteTolerance, &
        & quasiNewtonSolver%relativeTolerance,quasiNewtonSolver%solutionTolerance, &
        & quasiNewtonSolver%maximumNumberOfIterations,quasiNewtonSolver%maximumNumberOfFunctionEvaluations, &
        & err,error,*999)
      !Set any further SNES options from the command line options
      CALL PETSc_SNESSetFromOptions(trustregionSolver%snes,err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "// &
        & TRIM(NumberToVString(trustregionSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearQuasiNewtonTrustregion_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewtonTrustregion_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_CreateFinish

  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region delta0 for a nonlinear Quasi-Newton trust region solver solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonTrustRegionDelta0Set
  SUBROUTINE Solver_QuasiNewtonTrustregionDelta0Set(solver,trustregionDelta0,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the trust region delta0 for
    REAL(DP), INTENT(IN) :: trustregionDelta0 !<The trust region delta0 to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonTrustregionDelta0Set",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    CALL SolverNonlinearQuasiNewton_AssertIsTrustregion(quasiNewtonSolver,err,error,*999)
    NULLIFY(trustregionSolver)
    CALL SolverNonlinearQuasiNewton_TrustregionSolverGet(quasiNewtonSolver,trustregionSolver,err,error,*999)
    IF(trustregionDelta0<=ZERO_TOLERANCE) THEN
      localError="The specified trust region delta0 of "//TRIM(NumberToVString(trustregionDelta0,"*",err,error))// &
        & " is invalid. The trust region delta0 must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    trustregionSolver%trustregionDelta0=trustregionDelta0
    
    EXITS("Solver_QuasiNewtonTrustregionDelta0Set")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonTrustregionDelta0Set",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_QuasiNewtonTrustregionDelta0Set
        
  !
  !================================================================================================================================
  !
  
  !>Finalise a nonlinear Quasi-Newton trust region solver and deallocate all memory
  SUBROUTINE SolverNonlinearQuasiNewton_TrustregionFinalise(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer the non linear trust region solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("SolverNonlinearQuasiNewton_TrustregionFinalise",err,error,*999)

    IF(ASSOCIATED(trustregionSolver)) THEN      
      CALL PETSc_SNESFinalise(trustregionSolver%snes,err,error,*999)
      DEALLOCATE(trustregionSolver)
    ENDIF
    
    EXITS("SolverNonlinearQuasiNewton_TrustregionFinalise")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewton_TrustregionFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewton_TrustregionFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Quaso-Newton trust region solver for a nonlinear solver
  SUBROUTINE SolverNonlinearQuasiNewton_TrustregionInitialise(quasiNewtonSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver !<A pointer the Quasi-Newton solver to initialise the trust region solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
  
    ENTERS("SolverNonlinearQuasiNewton_TrustregionInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(quasiNewtonSolver)) CALL FlagError("Quasi-Newton solver is not associated.",err,error,*998)
    IF(ASSOCIATED(quasiNewtonSolver%trustregionSolver)) &
      & CALL FlagError("Trust region solver is already associated for this nonlinear solver.",err,error,*998)
     
    ALLOCATE(quasiNewtonSolver%trustregionSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Quasi-Newton solver trust region solver.",err,error,*999)
    quasiNewtonSolver%trustregionSolver%quasiNewtonSolver=>quasiNewtonSolver
    quasiNewtonSolver%trustregionSolver%solverLibrary=SOLVER_PETSC_LIBRARY
!!TODO: set this properly
    quasiNewtonSolver%trustregionSolver%trustregionDelta0=0.01_DP
    CALL PETSc_SNESInitialise(quasiNewtonSolver%trustregionSolver%snes,err,error,*999)
       
    EXITS("SolverNonlinearQuasiNewton_TrustregionInitialise")
    RETURN
999 CALL SolverNonlinearQuasiNewton_TrustregionFinalise(quasiNewtonSolver%trustregionSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverNonlinearQuasiNewton_TrustregionInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearQuasiNewton_TrustregionInitialise

  !
  !================================================================================================================================
  !

  !Solves a nonlinear Quasi-Newton trust region solver 
  SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_Solve(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the nonlinear Quasi-Newton trust region solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearQuasiNewtonTrustregion_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is not associated.",err,error,*999)
    
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinearQuasiNewtonTrustregion_QuasiNewtonSolverGet(trustregionSolver,quasiNewtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearQuasiNewton_NonlinearSolverGet(quasiNewtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    SELECT CASE(trustregionSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The nonlinear Quasi-Newton trust region solver library type of "// &
        & TRIM(NumberToVString(trustregionSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearQuasiNewtonTrustregion_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearQuasiNewtonTrustregion_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearQuasiNewtonTrustregion_Solve
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region tolerance for a nonlinear Quasi-Newton trust region solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonTrustRegionToleranceSet
  SUBROUTINE Solver_QuasiNewtonTrustRegionToleranceSet(solver,trustregionTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the trust region tolerance for
    REAL(DP), INTENT(IN) :: trustregionTolerance !<The trust region tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonTrustRegionToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    CALL SolverNonlinearQuasiNewton_AssertIsTrustregion(quasiNewtonSolver,err,error,*999)
    NULLIFY(trustregionSolver)
    CALL SolverNonlinearQuasiNewton_TrustregionSolverGet(quasiNewtonSolver,trustregionSolver,err,error,*999)
    IF(trustregionTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified trust region tolerance of "//TRIM(NumberToVString(trustregionTolerance,"*",err,error))// &
        & " is invalid. The trust region tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    trustregionSolver%trustregionTolerance=trustregionTolerance
    
    EXITS("Solver_QuasiNewtonTrustRegionToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonTrustRegionToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonTrustRegionToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the restart of nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonRestartSet
  SUBROUTINE Solver_QuasiNewtonRestartNumberSet(solver,restartNumber,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Quasi-Newton solver type
    INTEGER(INTG), INTENT(IN) :: restartNumber !<Sets the number of stored updates and the restart period for periodic restart type 
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    
    ENTERS("Solver_QuasiNewtonRestartNumberSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    quasiNewtonSolver%restart=restartNumber
    
    EXITS("Solver_QuasiNewtonRestartNumberSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonRestartNumberSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonRestartNumberSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the restart type of nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonRestartTypeSet
  SUBROUTINE Solver_QuasiNewtonRestartTypeSet(solver,quasiNewtonRestartType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Quasi-Newton solver type
    INTEGER(INTG), INTENT(IN) :: quasiNewtonRestartType !<The restart type of nonlinear Quasi-Newton to set \see SolverRoutines_QuasiNewtonRestartTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonRestartTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(quasiNewtonRestartType/=quasiNewtonSolver%restartType) THEN
      !Intialise the new type
      SELECT CASE(quasiNewtonRestartType)
      CASE(SOLVER_QUASI_NEWTON_RESTART_NONE)
        quasiNewtonSolver%restartType=SOLVER_QUASI_NEWTON_RESTART_NONE
      CASE(SOLVER_QUASI_NEWTON_RESTART_POWELL)
        quasiNewtonSolver%restartType=SOLVER_QUASI_NEWTON_RESTART_POWELL
      CASE(SOLVER_QUASI_NEWTON_RESTART_PERIODIC)
        quasiNewtonSolver%restartType=SOLVER_QUASI_NEWTON_RESTART_PERIODIC
      CASE DEFAULT
        localError="The Quasi-Newton restart type of "//TRIM(NumberToVString( &
          & quasiNewtonRestartType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_QuasiNewtonRestartTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonRestartTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonRestartTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the scale type of nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonScaleTypeSet
  SUBROUTINE Solver_QuasiNewtonScaleTypeSet(solver,quasiNewtonScaleType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Quasi-Newton solver type
    INTEGER(INTG), INTENT(IN) :: quasiNewtonScaleType !<The scale type of nonlinear Quasi-Newton to set \see SolverRoutines_QuasiNewtonScaleTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
   TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonScaleTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(quasiNewtonScaleType/=quasiNewtonSolver%scaleType) THEN
      !Intialise the new type
      SELECT CASE(quasiNewtonScaleType)
      CASE(SOLVER_QUASI_NEWTON_SCALE_NONE)
        quasiNewtonSolver%scaleType=SOLVER_QUASI_NEWTON_SCALE_NONE
      CASE(SOLVER_QUASI_NEWTON_SCALE_SHANNO)
        quasiNewtonSolver%scaleType=SOLVER_QUASI_NEWTON_SCALE_SHANNO
      CASE(SOLVER_QUASI_NEWTON_SCALE_LINESEARCH)
        quasiNewtonSolver%scaleType=SOLVER_QUASI_NEWTON_SCALE_LINESEARCH
      CASE(SOLVER_QUASI_NEWTON_SCALE_JACOBIAN)
        quasiNewtonSolver%scaleType=SOLVER_QUASI_NEWTON_SCALE_JACOBIAN
      CASE DEFAULT
        localError="The Quasi-Newton scale type of "//TRIM(NumberToVString(quasiNewtonScaleType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_QuasiNewtonScaleTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonScaleTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonScaleTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonTypeSet
  SUBROUTINE Solver_QuasiNewtonTypeSet(solver,quasiNewtonType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Quasi-Newton solver type
    INTEGER(INTG), INTENT(IN) :: quasiNewtonType !<The type of nonlinear Quasi-Newton to set \see SolverRoutines_QuasiNewtonTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_QuasiNewtonTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(quasiNewtonType/=quasiNewtonSolver%quasiNewtonType) THEN
      !Intialise the new type
      SELECT CASE(quasiNewtonType)
      CASE(SOLVER_QUASI_NEWTON_LBFGS)
        quasiNewtonSolver%quasiNewtonType=SOLVER_QUASI_NEWTON_LBFGS
      CASE(SOLVER_QUASI_NEWTON_GOODBROYDEN)
        quasiNewtonSolver%quasiNewtonType=SOLVER_QUASI_NEWTON_GOODBROYDEN
      CASE(SOLVER_QUASI_NEWTON_BADBROYDEN)
        quasiNewtonSolver%quasiNewtonType=SOLVER_QUASI_NEWTON_BADBROYDEN
      CASE DEFAULT
        localError="The Quasi-Newton type of "//TRIM(NumberToVString(quasiNewtonType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_QuasiNewtonTypeSet")
    RETURN
999 ERRORSEXITS("Solver_QuasiNewtonTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solve type of nonlinear Quasi-Newton solver. \see OpenCMISS::Iron::cmfe_Solver_QuasiNewtonSolveTypeSet
  SUBROUTINE Solver_QuasiNewtonSolveTypeSet(solver,quasiNewtonSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Quasi-Newton solver type
    INTEGER(INTG), INTENT(IN) :: quasiNewtonSolveType !<The type of nonlinear solver to set \see SolverRoutines_QuasiNewtonSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solver_QuasiNewtonSolveTypeSet",err,error,*998)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsQuasiNewton(nonlinearSolver,err,error,*999)
    NULLIFY(quasiNewtonSolver)
    CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
    IF(quasiNewtonSolveType/=quasiNewtonSolver%quasiNewtonSolveType) THEN
      !Intialise the new solver type
      SELECT CASE(quasiNewtonSolveType)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
        CALL SolverNonlinearQuasiNewton_LinesearchInitialise(quasiNewtonSolver,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
        CALL SolverNonlinearQuasiNewton_TrustregionInitialise(quasiNewtonSolver,err,error,*999)
      CASE DEFAULT
        localError="The Quasi-Newton solver type of "//TRIM(NumberToVString(quasiNewtonSolveType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old solver type
      SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
        CALL SolverNonlinearQuasiNewton_LinesearchFinalise(quasiNewtonSolver%linesearchSolver,err,error,*999)
      CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
        CALL SolverNonlinearQuasiNewton_TrustregionFinalise(quasiNewtonSolver%trustregionSolver,err,error,*999)
      CASE DEFAULT
        localError="The Quasi-Newton solver type of "// &
          & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      quasiNewtonSolver%quasiNewtonSolveType=quasiNewtonSolveType
    ENDIF
    
    EXITS("Solver_QuasiNewtonSolveTypeSet")
    RETURN
999 SELECT CASE(quasiNewtonSolveType)
    CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
      CALL SolverNonlinearQuasiNewton_LinesearchFinalise(quasiNewtonSolver%linesearchSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_QUASI_NEWTON_TRUSTREGION)
      CALL SolverNonlinearQuasiNewton_TrustregionFinalise(quasiNewtonSolver%trustregionSolver,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("Solver_QuasiNewtonSolveTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_QuasiNewtonSolveTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum absolute tolerance for a nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonAbsoluteToleranceSet
  SUBROUTINE Solver_NewtonAbsoluteToleranceSet(solver,absoluteTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the absolute tolerance for
    REAL(DP), INTENT(IN) :: absoluteTolerance !<The absolute tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonAbsoluteToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(absoluteTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified absolute tolerance of "//TRIM(NumberToVString(absoluteTolerance,"*",err,error))// &
        & " is invalid. The absolute tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    newtonSolver%absoluteTolerance=absoluteTolerance
    
    EXITS("Solver_NewtonAbsoluteToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonAbsoluteToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonAbsoluteToleranceSet

  !
  !================================================================================================================================
  !

  !>Enables/disables output monitoring for a nonlinear Newton line search solver.
  SUBROUTINE Solver_NewtonLineSearchMonitorOutputSet(solver,linesearchMonitorOutputFlag,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the absolute tolerance for
    LOGICAL, INTENT(IN) :: linesearchMonitorOutputFlag !<Flag to determine whether to enable/disable linsearch monitor output.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    
    ENTERS("Solver_NewtonLineSearchMonitorOutputSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
    linesearchSolver%linesearchMonitorOutput=linesearchMonitorOutputFlag
    
    EXITS("Solver_NewtonLineSearchMonitorOutputSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonLineSearchMonitorOutputSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Solver_NewtonLineSearchMonitorOutputSet

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a Newton solver 
  SUBROUTINE SolverNonlinearNewton_CreateFinish(newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the Newton solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewton_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(newtonSolver%newtonSolveType)
    CASE(SOLVER_NEWTON_LINESEARCH)
      CALL SolverNonlinearNewtonLinesearch_CreateFinish(newtonSolver%linesearchSolver,err,error,*999)
    CASE(SOLVER_NEWTON_TRUSTREGION)
      CALL SolverNonlinearNewtonTrustregion_CreateFinish(newtonSolver%trustregionSolver,err,error,*999)
    CASE DEFAULT
      localError="The Newton solver type of "// &
        & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearNewton_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewton_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a Newton solver and deallocate all memory
  RECURSIVE SUBROUTINE SolverNonlinear_NewtonFinalise(newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer the Newton solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverNonlinear_NewtonFinalise",err,error,*999)

    IF(ASSOCIATED(newtonSolver)) THEN
      CALL SolverNonlinearNewton_LinesearchFinalise(newtonSolver%linesearchSolver,err,error,*999)
      CALL SolverNonlinearNewton_TrustregionFinalise(newtonSolver%trustregionSolver,err,error,*999)
      CALL Solver_Finalise(newtonSolver%linearSolver,err,error,*999)
      DEALLOCATE(newtonSolver)
    ENDIF
         
    EXITS("SolverNonlinear_NewtonFinalise")
    RETURN
999 ERRORSEXITS("SolverNonlinear_NewtonFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_NewtonFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Newton solver for a nonlinear solver
  SUBROUTINE SolverNonlinear_NewtonInitialise(nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer the solver to initialise the Newton solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(SolverType), POINTER :: solver
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("SolverNonlinear_NewtonInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*998)
    IF(ASSOCIATED(nonlinearSolver%newtonSolver)) &
      & CALL FlagError("Newton solver is already associated for this nonlinear solver.",err,error,*998)

    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    !Allocate and initialise a Newton solver
    ALLOCATE(nonlinearSolver%newtonSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nonlinear solver Newton solver.",err,error,*999)
    nonlinearSolver%newtonSolver%nonlinearSolver=>nonlinearSolver
    nonlinearSolver%newtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    nonlinearSolver%newtonSolver%totalNumberOfFunctionEvaluations=0
    nonlinearSolver%newtonSolver%totalNumberOfJacobianEvaluations=0
    nonlinearSolver%newtonSolver%maximumNumberOfIterations=50
    nonlinearSolver%newtonSolver%maximumNumberOfFunctionEvaluations=1000
    nonlinearSolver%newtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
    nonlinearSolver%newtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT
    nonlinearSolver%newtonSolver%absoluteTolerance=1.0E-10_DP
    nonlinearSolver%newtonSolver%relativeTolerance=1.0E-05_DP
    nonlinearSolver%newtonSolver%solutionTolerance=1.0E-05_DP
    NULLIFY(nonlinearSolver%newtonSolver%linesearchSolver)
    NULLIFY(nonlinearSolver%newtonSolver%trustregionSolver)
    NULLIFY(nonlinearSolver%newtonSolver%cellMLEvaluatorSolver)
    NULLIFY(nonlinearSolver%newtonSolver%convergenceTest)
    ALLOCATE(nonlinearSolver%newtonSolver%convergenceTest,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate convergence test object.",err,error,*999)
    nonlinearSolver%newtonSolver%convergenceTest%energyFirstIter = 0.0_DP
    nonlinearSolver%newtonSolver%convergenceTest%normalisedEnergy = 0.0_DP
    !Default to a Newton linesearch solver
    nonlinearSolver%newtonSolver%newtonSolveType=SOLVER_NEWTON_LINESEARCH
    CALL SolverNonlinearNewton_LinesearchInitialise(nonlinearSolver%newtonSolver,err,error,*999)
    !Create the linked linear solver
    ALLOCATE(nonlinearSolver%newtonSolver%linearSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Newton solver linear solver.",err,error,*999)
    NULLIFY(nonlinearSolver%newtonSolver%linearSolver%solvers)
    CALL Solver_Initialise(nonlinearSolver%newtonSolver%linearSolver,err,error,*999)
    CALL Solver_LinearInitialise(nonlinearSolver%newtonSolver%linearSolver,err,error,*999)
    CALL Solver_LinkedSolverAdd(solver,nonlinearSolver%newtonSolver%linearSolver,SOLVER_LINEAR_TYPE,err,error,*999)
        
    EXITS("SolverNonlinear_NewtonInitialise")
    RETURN
999 CALL SolverNonlinear_NewtonFinalise(nonlinearSolver%newtonSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverNonlinear_NewtonInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_NewtonInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of Jacobian calculation type for a Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonJacobianCalculationSet
  SUBROUTINE Solver_NewtonJacobianCalculationTypeSet(solver,jacobianCalculationType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the Jacobian calculation type
    INTEGER(INTG), INTENT(IN) :: jacobianCalculationType !<The type of Jacobian calculation type to set \see SolverRoutines_JacobianCalculationTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonJacobianCalculationTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(jacobianCalculationType/=newtonSolver%jacobianCalculationType) THEN
      SELECT CASE(jacobianCalculationType)
      CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
        newtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED
      CASE(SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED)
        newtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED
      CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
        newtonSolver%jacobianCalculationType=SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
      CASE DEFAULT
        localError="The Jacobian calculation type of "// &
          & TRIM(NumberToVString(jacobianCalculationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_NewtonJacobianCalculationTypeSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonJacobianCalculationTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonJacobianCalculationTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a Newton solver.
  SUBROUTINE SolverNonlinearNewton_LibraryTypeSet(newtonSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer the Newton solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the Newton solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewton_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(newtonSolver%newtonSolveType)
    CASE(SOLVER_NEWTON_LINESEARCH)
      NULLIFY(linesearchSolver)
      CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        linesearchSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Newton linesearch solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(SOLVER_NEWTON_TRUSTREGION)
      NULLIFY(trustregionSolver)
      CALL SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,trustregionSolver,err,error,*999)
      SELECT CASE(solverLibraryType)
      CASE(SOLVER_CMISS_LIBRARY)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_PETSC_LIBRARY)
        trustregionSolver%solverLibrary=SOLVER_PETSC_LIBRARY
      CASE DEFAULT
        localError="The solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
          & " is invalid for a Newton trustregion solver."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The Newton solver type of "// &
        & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearNewton_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewton_LibraryTypeSet

  !
  !================================================================================================================================
  !
  
  !>Sets/changes the convergence test for a Newton nonlinear solver \see OpenCMISS::Iron::cmfe_Solver_NewtonConvergenceTestSet
  SUBROUTINE Solver_NewtonConvergenceTestTypeSet(solver,convergenceTestType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the convergence test for
    INTEGER(INTG), INTENT(IN) :: convergenceTestType !<The convergence test type to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver 
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonConvergenceTestTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    SELECT CASE(convergenceTestType)
    CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
      newtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT
    CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
      newtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM
    CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
      newtonSolver%convergenceTestType=SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO
    CASE DEFAULT
      localError="The specified convergence test type of "//TRIM(NumberToVString(convergenceTestType, &
        & "*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_NewtonConvergenceTestTypeSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonConvergenceTestTypeSet",err,error)    
    RETURN 1
   
  END SUBROUTINE Solver_NewtonConvergenceTestTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the line search alpha for a Newton linesearch solver \see OpenCMISS::Iron::cmfe_Solver_NewtonLineSearchAlphaSet
  SUBROUTINE Solver_NewtonLinesearchAlphaSet(solver,linesearchAlpha,err,error,*)
    
    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search alpha for
    REAL(DP), INTENT(IN) :: linesearchAlpha !<The line search alpha to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonLinesearchAlphaSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
    IF(linesearchAlpha<=ZERO_TOLERANCE) THEN
      localError="The specified line search alpha of "//TRIM(NumberToVString(linesearchAlpha,"*",err,error))// &
        & " is invalid. The line search alpha must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    linesearchSolver%linesearchAlpha=linesearchAlpha
   
    EXITS("Solver_NewtonLinesearchAlphaSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonLinesearchAlphaSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonLinesearchAlphaSet
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Newton line search solver
  SUBROUTINE SolverNonlinearNewtonLinesearch_CreateFinish(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer the nonlinear Newton line search solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx,interfaceConditionIdx,interfaceMatrixIdx,groupCommunicator, &
      & numberOfEquationsSets,numberOfInterfaceConditions,numberOfInterfaceMatrices,numberOfLinearMatrices, &
      & numberOfSolverMatrices,sparsityType,symmetryType
    EXTERNAL :: Problem_SolverJacobianEvaluatePetsc
    EXTERNAL :: Problem_SolverJacobianFDCalculatePetsc
    EXTERNAL :: Problem_SolverResidualEvaluatePetsc
    EXTERNAL :: Problem_SolverConvergenceTestPetsc
    EXTERNAL :: Problem_SolverNonlinearMonitorPetsc
    TYPE(DistributedMatrixType), POINTER :: jacobianMatrix
    TYPE(DistributedMatrixPETScType), POINTER :: jacobianPETScMatrix
    TYPE(DistributedVectorType), POINTER :: residualVector
    TYPE(DistributedVectorPETScType), POINTER :: residualPETScVector
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: linearVariable,interfaceVariable,lagrangeVariable,lhsVariable
    TYPE(LinearDirectSolverType), POINTER :: directSolver
    TYPE(LinearIterativeSolverType), POINTER :: iterativeSolver
    TYPE(LinearSolverType), POINTER :: linearSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: linkedLinearSolver,solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverJacobian
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("SolverNonlinearNewtonLinesearch_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)
    
    NULLIFY(newtonSolver)
    CALL SolverNonlinearNewtonLinesearch_NewtonSolverGet(linesearchSolver,newtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearNewton_NonlinearSolverGet(newtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    SELECT CASE(linesearchSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Loop over the equations set in the solver equations
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(lhsMapping)
        CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
        NULLIFY(lhsVariable)
        CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lhsVariable,domainMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          !If there are any linear matrices create temporary vector for matrix-vector products
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
          DO equationsMatrixIdx=1,numberOfLinearMatrices
            NULLIFY(equationsMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
            IF(.NOT.ASSOCIATED(equationsMatrix%tempVector)) THEN
              CALL DistributedVector_CreateStart(domainMapping,equationsMatrix%tempVector,err,error,*999)
              CALL DistributedVector_DataTypeSet(equationsMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
              CALL DistributedVector_CreateFinish(equationsMatrix%tempVector,err,error,*999)
            ENDIF
          ENDDO !equationsMatrixIdx
        ENDIF
      ENDDO !equationsSetIdx
      !Loop over the interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceCondition)
        CALL SolverMapping_InterfaceConditionGet(solverMapping,interfaceConditionIdx,interfaceCondition,err,error,*999)
        NULLIFY(lagrangeField)
        CALL InterfaceCondition_LagrangeFieldGet(interfaceCondition,lagrangeField,err,error,*999)
        NULLIFY(interfaceEquations)
        CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
        NULLIFY(interfaceMatrices)
        CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
        NULLIFY(interfaceMapping)
        CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
        !Create temporary vector for matrix-vector products
        CALL InterfaceMapping_NumberOfInterfaceMatricesGet(interfaceMapping,numberOfInterfaceMatrices,err,error,*999)
        DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
          IF(.NOT.ASSOCIATED(interfaceMatrix%tempVector)) THEN
            NULLIFY(interfaceVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,interfaceMatrixIdx,interfaceVariable,err,error,*999)
            NULLIFY(domainMapping)
            CALL FieldVariable_DomainMappingGet(interfaceVariable,domainMapping,err,error,*999)
            !Set up the temporary interface distributed vector to be used with interface matrices
            CALL DistributedVector_CreateStart(domainMapping,interfaceMatrix%tempVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(interfaceMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
            CALL DistributedVector_CreateFinish(interfaceMatrix%tempVector,err,error,*999)
            !Set up the temporary interface distributed vector to be used with transposed interface matrices
            NULLIFY(domainMapping)
            CALL FieldVariable_DomainMappingGet(lagrangeVariable,domainMapping,err,error,*999)
            CALL DistributedVector_CreateStart(domainMapping,interfaceMatrix%tempTransposeVector,err,error,*999)
            CALL DistributedVector_DataTypeSet(interfaceMatrix%tempTransposeVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE, &
              & err,error,*999)
            CALL DistributedVector_CreateFinish(interfaceMatrix%tempTransposeVector,err,error,*999)
          ENDIF
        ENDDO !interfaceMatrixIdx
      ENDDO !interfaceConiditionIdx
      !Create the PETSc SNES solver
      CALL PETSc_SNESCreate(groupCommunicator,linesearchSolver%snes,err,error,*999)
      !Set the nonlinear solver type to be a Newton line search solver
      CALL PETSc_SNESSetType(linesearchSolver%snes,PETSC_SNESNEWTONLS,err,error,*999)
      !Create the solver matrices and vectors
      NULLIFY(linkedLinearSolver)
      CALL SolverNonlinearNewton_LinkedLinearSolverGet(newtonSolver,linkedLinearSolver,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
      CALL SolverMatrices_LibraryTypeSet(solverMatrices,SOLVER_PETSC_LIBRARY,err,error,*999)
      CALL SolverEquations_SparsityTypeGet(solverEquations,sparsityType,err,error,*999)
      SELECT CASE(sparsityType)
      CASE(SOLVER_SPARSE_MATRICES)
        CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
      CASE(SOLVER_FULL_MATRICES)
        CALL SolverMatrices_StorageTypesSet(solverMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverEquations_SymmetryTypeGet(solverEquations,symmetryType,err,error,*999)
      SELECT CASE(symmetryType)
      CASE(SOLVER_SYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
      CASE(SOLVER_UNSYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations symmetry type of "//TRIM(NumberToVString(symmetryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      !Link linear solver
      linkedLinearSolver%solverEquations=>solver%solverEquations
      !Finish the creation of the linear solver
      NULLIFY(linearSolver)
      CALL Solver_LinearSolverGet(linkedLinearSolver,linearSolver,err,error,*999)
      CALL SolverLinear_CreateFinish(linearSolver,err,error,*999)
      !Associate linear solver's KSP to nonlinear solver's SNES
      SELECT CASE(linearSolver%linearSolveType)
      CASE(SOLVER_LINEAR_DIRECT_SOLVE_TYPE)
        NULLIFY(directSolver)
        CALL SolverLinear_DirectSolverGet(linearSolver,directSolver,err,error,*999)
        CALL PETSc_SNESSetKsp(linesearchSolver%snes,directSolver%ksp,err,error,*999)
      CASE(SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE)
        NULLIFY(iterativeSolver)
        CALL SolverLinear_IterativeSolverGet(linearSolver,iterativeSolver,err,error,*999)
        CALL PETSc_SNESSetKSP(linesearchSolver%snes,iterativeSolver%ksp,err,error,*999)
      CASE DEFAULT
        localError="The linear solver type of "//TRIM(NumberToVString(linearSolver%linearSolveType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the nonlinear function
      NULLIFY(residualVector)
      CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualVector,err,error,*999)
      NULLIFY(residualPETScVector)
      CALL DistributedVector_PETScVectorGet(residualVector,residualPETScVector,err,error,*999)
      !Set the solver as a context for the SNES object
      CALL PETSc_SNESSetApplicationContext(linesearchSolver%snes,solver, err,error,*999)
      !Pass the linesearch solver object rather than the temporary solver
      CALL PETSc_SNESSetFunction(linesearchSolver%snes,residualPETScVector%vector,Problem_SolverResidualEvaluatePetsc, &
        & solver,err,error,*999)
      SELECT CASE(newtonSolver%convergenceTestType)
      CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
        !Default convergence test, do nothing
      CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM,SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
        CALL PETSc_SNESSetConvergenceTest(linesearchSolver%snes,Problem_SolverConvergenceTestPetsc,solver,err,error,*999)
      CASE DEFAULT
        localError="The specified convergence test type of "// &
          & TRIM(NumberToVString(newtonSolver%convergenceTestType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT      
      !Set the Jacobian
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
      IF(numberOfSolverMatrices/=1) THEN
        localError="Invalid number of solver matrices. The number of solver matrices is "// &
          & TRIM(NumberToVString(numberOfSolverMatrices,"*",err,error))//" and it should be 1."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(solverJacobian)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverJacobian,err,error,*999)
      NULLIFY(jacobianMatrix)
      CALL SolverMatrix_SolverDistributedMatrixGet(solverJacobian,jacobianMatrix,err,error,*999)
      NULLIFY(jacobianPETScMatrix)
      CALL DistributedMatrix_PETScMatrixGet(jacobianMatrix,jacobianPETScMatrix,err,error,*999)
      SELECT CASE(newtonSolver%jacobianCalculationType)
      CASE(SOLVER_NEWTON_JACOBIAN_NOT_CALCULATED)
        CALL FlagError("Cannot have no Jacobian calculation for a PETSc nonlinear linesearch solver.",err,error,*999)
      CASE(SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED)
        solverJacobian%updateMatrix=.TRUE. !CMISS will fill in the Jacobian values
        !Pass the linesearch solver object rather than the temporary solver
        CALL PETSc_SNESSetJacobian(linesearchSolver%snes,jacobianPETScMatrix%matrix,jacobianPETScMatrix%matrix, &
          & Problem_SolverJacobianEvaluatePetsc,linesearchSolver%newtonSolver%nonlinearSolver%solver,err,error,*999)
      CASE(SOLVER_NEWTON_JACOBIAN_FD_CALCULATED)
        solverJacobian%updateMatrix=.FALSE. !Petsc will fill in the Jacobian values
        CALL DistributedMatrix_Form(jacobianMatrix,err,error,*999)        
        SELECT CASE(sparsityType)
        CASE(SOLVER_SPARSE_MATRICES)
          CALL PETSc_MatColoringCreate(jacobianPETScMatrix%matrix,linesearchSolver%jacobianMatColoring,err,error,*999)
          CALL PETSc_MatColoringSetType(linesearchSolver%jacobianMatColoring,PETSC_MATCOLORING_SL,err,error,*999)
          CALL PETSc_MatColoringSetFromOptions(linesearchSolver%jacobianMatColoring,err,error,*999)
          CALL PETSc_MatColoringApply(linesearchSolver%jacobianMatColoring,linesearchSolver%jacobianISColoring,err,error,*999)
          CALL PETSc_MatColoringDestroy(linesearchSolver%jacobianMatColoring,err,error,*999)
          !Compute SNESComputeJacobianDefaultColor data structure
          CALL PETSc_MatFDColoringCreate(jacobianMatrix%petsc%matrix,linesearchSolver%jacobianISColoring, &
            & linesearchSolver%jacobianMatFDColoring,err,error,*999)
          !Pass the linesearch solver object rather than the temporary solver
          CALL PETSc_MatFDColoringSetFunction(linesearchSolver%jacobianMatFDColoring,Problem_SolverResidualEvaluatePetsc, &
            & linesearchSolver%newtonSolver%nonlinearSolver%solver,err,error,*999)
          CALL PETSc_MatFDColoringSetFromOptions(linesearchSolver%jacobianMatFDColoring,err,error,*999)
          CALL PETSc_MatFDColoringSetup(jacobianMatrix%petsc%matrix,linesearchSolver%jacobianISColoring, &
            & linesearchSolver%jacobianMatFDColoring,err,error,*999)
          CALL PETSc_ISColoringDestroy(linesearchSolver%jacobianISColoring,err,error,*999)
        CASE(SOLVER_FULL_MATRICES)
          !Do nothing
        CASE DEFAULT
          localError="The specified solver equations sparsity type of "// &
            & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL PETSc_SNESSetJacobian(linesearchSolver%snes,jacobianPETScMatrix%matrix,jacobianPETScMatrix%matrix, &
          & Problem_SolverJacobianFDCalculatePetsc,solver,err,error,*999)
      CASE DEFAULT
        localError="The Jacobian calculation type of "// &
          & TRIM(NumberToVString(newtonSolver%jacobianCalculationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      IF(solver%outputType>=SOLVER_MONITOR_OUTPUT) THEN
        !Set the monitor
        !Pass the linesearch solver object rather than the temporary solver
        CALL PETSc_SNESMonitorSet(linesearchSolver%snes,Problem_SolverNonlinearMonitorPETSC,solver,err,error,*999)
      ENDIF
      CALL PETSc_SNESGetLineSearch(linesearchSolver%snes,linesearchSolver%snesLineSearch,err,error,*999)
      !Set the line search type and order where applicable
      SELECT CASE(linesearchSolver%linesearchType)
      CASE(SOLVER_NEWTON_LINESEARCH_NONORMS)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_BASIC,err,error,*999)
        CALL PETSc_SNESLineSearchSetComputeNorms(linesearchSolver%snesLineSearch,.FALSE.,err,error,*999)
      CASE(SOLVER_NEWTON_LINESEARCH_LINEAR)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_CP,err,error,*999)
        CALL PETSc_SNESLineSearchSetOrder(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_ORDER_LINEAR,err,error,*999)
      CASE(SOLVER_NEWTON_LINESEARCH_QUADRATIC)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_BT,err,error,*999)
        CALL PETSc_SNESLineSearchSetOrder(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_ORDER_QUADRATIC,err,error,*999)
      CASE(SOLVER_NEWTON_LINESEARCH_CUBIC)
        CALL PETSc_SNESLineSearchSetType(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_BT,err,error,*999)
        CALL PETSc_SNESLineSearchSetOrder(linesearchSolver%snesLineSearch,PETSC_SNES_LINESEARCH_ORDER_CUBIC,err,error,*999)
      CASE DEFAULT
        localError="The nonlinear Newton line search type of "// &
          & TRIM(NumberToVString(linesearchSolver%linesearchType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(linesearchSolver%linesearchType)
      CASE(SOLVER_NEWTON_LINESEARCH_QUADRATIC,SOLVER_NEWTON_LINESEARCH_CUBIC)
        !Alpha parameter only applicable for back-tracking linesearch
        CALL PETSc_SNESLineSearchBTSetAlpha(linesearchSolver%snesLineSearch,linesearchSolver%linesearchAlpha,err,error,*999)
      END SELECT
      !Set step tolerances, leave iterative line search options as defaults.
!!TODO: set the rtol, atol, ltol and maxits properly.
      CALL PETSc_SNESLineSearchSetTolerances(linesearchSolver%snesLineSearch,linesearchSolver%linesearchStepTolerance, &
        & linesearchSolver%linesearchMaxstep,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER, &
        & err,error,*999)
      IF(linesearchSolver%linesearchMonitorOutput) THEN
        CALL PETSc_SNESLineSearchSetMonitor(linesearchSolver%snesLineSearch,PETSC_TRUE,err,error,*999)
      ELSE
        CALL PETSc_SNESLineSearchSetMonitor(linesearchSolver%snesLineSearch,PETSC_FALSE,err,error,*999)
      ENDIF
      !Set the tolerances for the SNES solver
      CALL PETSc_SNESSetTolerances(linesearchSolver%snes,newtonSolver%absoluteTolerance,newtonSolver%relativeTolerance, &
        & newtonSolver%solutionTolerance,newtonSolver%maximumNumberOfIterations,newtonSolver%maximumNumberOfFunctionEvaluations, &
        & err,error,*999)            
      !Set any further SNES options from the command line options
      CALL PETSc_SNESSetFromOptions(linesearchSolver%snes,err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "// &
        & TRIM(NumberToVString(linesearchSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinearNewtonLinesearch_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewtonLinesearch_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewtonLinesearch_CreateFinish

  !
  !================================================================================================================================
  !

  !>Finalise a nonlinear Newton line search solver and deallocate all memory
  SUBROUTINE SolverNonlinearNewton_LinesearchFinalise(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer the nonlinear Newton line search solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("SolverNonlinearNewton_LinesearchFinalise",err,error,*999)

    IF(ASSOCIATED(linesearchSolver)) THEN
      CALL PETSc_MatColoringFinalise(linesearchSolver%jacobianMatColoring,err,error,*999)
      CALL PETSc_ISColoringFinalise(linesearchSolver%jacobianISColoring,err,error,*999)
      CALL PETSc_MatFDColoringFinalise(linesearchSolver%jacobianMatFDColoring,err,error,*999)
      CALL PETSc_SNESLineSearchFinalise(linesearchSolver%snesLineSearch,err,error,*999)
      CALL PETSc_SNESFinalise(linesearchSolver%snes,err,error,*999)
      DEALLOCATE(linesearchSolver)
    ENDIF
        
    EXITS("SolverNonlinearNewton_LinesearchFinalise")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_LinesearchFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_LinesearchFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a nonlinear Newton line search solver for a Newton solver
  SUBROUTINE SolverNonlinearNewton_LinesearchInitialise(newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer the nonlinear Newton solver to initialise the Newton line search solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
  
    ENTERS("SolverNonlinearNewton_LinesearchInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*998)
    IF(ASSOCIATED(newtonSolver%linesearchSolver)) &
      & CALL FlagError("Netwon line search solver is already associated for this Newton solver.",err,error,*998)
      
    !Allocate and initialise the Newton linesearch solver
    ALLOCATE(newtonSolver%linesearchSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate nonlinear solver Newton line search solver.",err,error,*999)
    newtonSolver%linesearchSolver%newtonSolver=>newtonSolver
    newtonSolver%linesearchSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    newtonSolver%linesearchSolver%linesearchType=SOLVER_NEWTON_LINESEARCH_CUBIC
    newtonSolver%linesearchSolver%linesearchAlpha=0.0001_DP
    newtonSolver%linesearchSolver%linesearchMaxstep=1.0E8_DP
    newtonSolver%linesearchSolver%linesearchStepTolerance=CONVERGENCE_TOLERANCE
    CALL PETSc_MatColoringInitialise(newtonSolver%linesearchSolver%jacobianMatColoring,err,error,*999)
    CALL PETSc_ISColoringInitialise(newtonSolver%linesearchSolver%jacobianISColoring,err,error,*999)
    CALL PETSc_MatFDColoringInitialise(newtonSolver%linesearchSolver%jacobianMatFDColoring,err,error,*999)
    CALL PETSc_SNESInitialise(newtonSolver%linesearchSolver%snes,err,error,*999)
    CALL PETSc_SNESLineSearchInitialise(newtonSolver%linesearchSolver%snesLineSearch,err,error,*999)
    newtonSolver%linesearchSolver%linesearchMonitorOutput=.FALSE.
        
    EXITS("SolverNonlinearNewton_LinesearchInitialise")
    RETURN
999 CALL SolverNonlinearNewton_LinesearchFinalise(newtonSolver%linesearchSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverNonlinearNewton_LinesearchInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewton_LinesearchInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the line search maximum step for a nonlinear Newton linesearch solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonLineSearchMaxStepSet
  SUBROUTINE Solver_NewtonLinesearchMaxStepSet(solver,linesearchMaximumStep,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search maximum step for
    REAL(DP), INTENT(IN) :: linesearchMaximumStep !<The line search maximum step to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonLinesearchMaxStepSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
    IF(linesearchMaximumStep<=ZERO_TOLERANCE) THEN
      localError="The specified line search maximum step of "// &
        & TRIM(NumberToVString(linesearchMaximumStep,"*",err,error))// &
        & " is invalid. The line search maximum step must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    linesearchSolver%linesearchMaxstep=linesearchMaximumStep
    
    EXITS("Solver_NewtonLinesearchMaxStepSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonLinesearchMaxStepSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonLinesearchMaxStepSet
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton line search solver 
  SUBROUTINE SolverNonlinearNewtonLinesearch_Solve(linesearchSolver,err,error,*)

    !Argument variables
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver !<A pointer to the nonlinear Newton line search solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: convergedReason,numberOfIterations,numberOfSolverMatrices
    REAL(DP) :: functionNorm
    TYPE(DistributedVectorType), POINTER :: rhsVector,solverVector
    TYPE(DistributedVectorPETscType), POINTER :: rhsPETScVector,solverPETScVector
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(PetscVecType) :: functionPETScVector
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewtonLinesearch_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(linesearchSolver)) CALL FlagError("Linesearch solver is not associated.",err,error,*999)
    
    NULLIFY(newtonSolver)
    CALL SolverNonlinearNewtonLinesearch_NewtonSolverGet(linesearchSolver,newtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearNewton_NonlinearSolverGet(newtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    IF(numberOfSolverMatrices/=1) THEN
      localError="The number of solver matrices of "// &
        & TRIM(NumberToVString(numberOfSolverMatrices,"*",err,error))// &
        & " is invalid. There should only be one solver matrix for a Newton linesearch solver."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(solverMatrix)
    CALL SolverMatrices_SolverMatrixGet(solverMatrices,1,solverMatrix,err,error,*999)
    NULLIFY(rhsVector)
    CALL SolverMatrices_RHSDistributedVectorGet(solverMatrices,rhsVector,err,error,*999)
    NULLIFY(solverVector)
    CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
    SELECT CASE(linesearchSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      SELECT CASE(newtonSolver%solutionInitialiseType)
      CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
        !Zero the solution vector
        CALL DistributedVector_AllValuesSet(solverVector,0.0_DP,err,error,*999)
      CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
        !Make sure the solver vector contains the current dependent field values
        CALL Solver_SolutionUpdate(solver,err,error,*999)
      CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
        !Do nothing
      CASE DEFAULT
        localError="The Newton solver solution initialise type of "// &
          & TRIM(NumberToVString(newtonSolver%solutionInitialiseType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      NULLIFY(rhsPETScVector)
      CALL DistributedVector_PETScVectorGet(rhsVector,rhsPETScVector,err,error,*999)
      NULLIFY(solverPETScVector)
      CALL DistributedVector_PETScVectorGet(solverVector,solverPETScVector,err,error,*999)
      !Set step tolerances, leave iterative line search options as defaults in case the user has changed them
      CALL PETSc_SNESLineSearchSetTolerances(linesearchSolver%snesLineSearch,linesearchSolver%linesearchStepTolerance, &
        & linesearchSolver%linesearchMaxstep,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER, &
        & err,error,*999)
      !Set the tolerances for the SNES solver in case the user has changed them
      CALL PETSc_SNESSetTolerances(linesearchSolver%snes,newtonSolver%absoluteTolerance,newtonSolver%relativeTolerance, &
        & newtonSolver%solutionTolerance,newtonSolver%maximumNumberOfIterations,newtonSolver%maximumNumberOfFunctionEvaluations, &
        & err,error,*999)     
      !Solve the nonlinear equations
      CALL PETSc_SNESSolve(linesearchSolver%snes,rhsPETScVector%vector,solverPETScVector%vector,err,error,*999)
      !Check for convergence
      CALL PETSc_SNESGetConvergedReason(linesearchSolver%snes,convergedReason,err,error,*999)
      SELECT CASE(convergedReason)
      CASE(PETSC_SNES_DIVERGED_FUNCTION_DOMAIN)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged function domain.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_FUNCTION_COUNT)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged function count.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LINEAR_SOLVE)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged linear solve.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_FNORM_NAN)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged F Norm NaN.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_MAX_IT)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged maximum iterations.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LINE_SEARCH)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged line search.",err,error,*999)
      CASE(PETSC_SNES_DIVERGED_LOCAL_MIN)
        CALL FlagError("Nonlinear line search solver did not converge. PETSc diverged local minimum.",err,error,*999)
      END SELECT
      IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
        !Output solution characteristics
        CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Newton linesearch solver parameters:",err,error,*999)
        CALL PETSc_SNESGetIterationNumber(linesearchSolver%snes,numberOfIterations,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final number of iterations = ",numberOfIterations,err,error,*999)
        CALL PETSc_SNESGetFunction(linesearchSolver%snes,functionPETScVector,err,error,*999)
        CALL PETSc_VecNorm(functionPETScVector,PETSC_NORM_2,functionNorm,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Final function norm = ",functionNorm,err,error,*999)
        SELECT CASE(convergedReason)
        CASE(PETSC_SNES_CONVERGED_FNORM_ABS)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm absolute.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_FNORM_RELATIVE)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged F Norm relative.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_SNORM_RELATIVE)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged S Norm relative.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_ITS)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged its.",err,error,*999)
        CASE(PETSC_SNES_CONVERGED_ITERATING)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Converged Reason = PETSc converged iterating.",err,error,*999)
        END SELECT
      ENDIF
    CASE DEFAULT
      localError="The Newton line search solver library type of "// &
        & TRIM(NumberToVString(linesearchSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinearNewtonLinesearch_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewtonLinesearch_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewtonLinesearch_Solve
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search step tolerance for a nonlinear Newton line search solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonLineSearchStepTolSet
  SUBROUTINE Solver_NewtonLinesearchStepToleranceSet(solver,linesearchStepTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search step tolerance for
    REAL(DP), INTENT(IN) :: linesearchStepTolerance !<The line search step tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonLinesearchStepToleranceSet",err,error,*999)
    
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
    IF(linesearchStepTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified line search step tolerance of "// &
        & TRIM(NumberToVString(linesearchStepTolerance,"*",err,error))// &
        & " is invalid. The line search step tolerance must be > 0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    linesearchSolver%linesearchStepTolerance=linesearchStepTolerance
    
    EXITS("Solver_NewtonLinesearchStepToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonLinesearchStepToleranceSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonLinesearchStepToleranceSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the line search type for a nonlinear Newton linesearch solver \see OpenCMISS::Iron::cmfe_Solver_NewtonLineSearchTypeSet
  SUBROUTINE Solver_NewtonLinesearchTypeSet(solver,linesearchType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the line search type for
    INTEGER(INTG), INTENT(IN) :: linesearchType !<The line search type to set \see SolverRoutines_NewtonLineSearchTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: linesearchSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonLinesearchTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)    
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsLinesearch(newtonSolver,err,error,*999)
    NULLIFY(linesearchSolver)
    CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,linesearchSolver,err,error,*999)
    SELECT CASE(linesearchType)
    CASE(SOLVER_NEWTON_LINESEARCH_NONORMS)
      linesearchSolver%linesearchType=SOLVER_NEWTON_LINESEARCH_NONORMS
    CASE(SOLVER_NEWTON_LINESEARCH_LINEAR)
      linesearchSolver%linesearchType=SOLVER_NEWTON_LINESEARCH_LINEAR
    CASE(SOLVER_NEWTON_LINESEARCH_QUADRATIC)
      linesearchSolver%linesearchType=SOLVER_NEWTON_LINESEARCH_QUADRATIC
    CASE(SOLVER_NEWTON_LINESEARCH_CUBIC)
      linesearchSolver%linesearchType=SOLVER_NEWTON_LINESEARCH_CUBIC
    CASE DEFAULT
      localError="The specified line search type of "//TRIM(NumberToVString(linesearchType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_NewtonLinesearchTypeSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonLinesearchTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonLinesearchTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of function evaluations for a nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonMaximumFunctionEvaluationsSet
  SUBROUTINE Solver_NewtonMaximumFunctionEvaluationsSet(solver,maximumFunctionEvaluations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the maximum function evaluations for
    INTEGER(INTG), INTENT(IN) :: maximumFunctionEvaluations !<The maximum function evaluations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonMaximumFunctionEvaluationsSet",err,error,*999)
    
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(maximumFunctionEvaluations<1) THEN
      localError="The specified maximum number of function evaluations of "// &
        & TRIM(NumberToVString(maximumFunctionEvaluations,"*",err,error))// &
        & " is invalid. The maximum number of function evaluations must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    newtonSolver%maximumNumberOfFunctionEvaluations=maximumFunctionEvaluations
    
    EXITS("Solver_NewtonMaximumFunctionEvaluationsSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonMaximumFunctionEvaluationsSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonMaximumFunctionEvaluationsSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the maximum number of iterations for a nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonMaximumIterationsSet
  SUBROUTINE Solver_NewtonMaxNumberOfIterationsSet(solver,maxNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the maximum iterations for
    INTEGER(INTG), INTENT(IN) :: maxNumberOfIterations !<The maximum iterations to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonMaxNumberOfIterationsSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(maxNumberOfIterations<1) THEN
      localError="The specified maximum iterations of "//TRIM(NumberToVString(maxNumberOfIterations,"*",err,error))// &
        & " is invalid. The maximum number of iterations must be >= 1."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    newtonSolver%maximumNumberOfIterations=maxNumberOfIterations
    
    EXITS("Solver_NewtonMaxNumberOfIterationsSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonMaxNumberOfIterationsSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonMaxNumberOfIterationsSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the relative tolerance for a nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonRelativeToleranceSet
  SUBROUTINE Solver_NewtonRelativeToleranceSet(solver,relativeTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the relative tolerance for
    REAL(DP), INTENT(IN) :: relativeTolerance !<The relative tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonRelativeToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(relativeTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified relative tolerance of "//TRIM(NumberToVString(relativeTolerance,"*",err,error))// &
        & " is invalid. The relative tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    newtonSolver%relativeTolerance=relativeTolerance
    
    EXITS("Solver_NewtonRelativeToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonRelativeToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonRelativeToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution initialisation for a nonlinear Newton solver
  SUBROUTINE Solver_NewtonSolutionInitialiseTypeSet(solver,solutionInitialiseType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the solution tolerance for
    INTEGER(INTG), INTENT(IN) :: solutionInitialiseType !<The solution initialise type to set \see SolverRoutines_SolutionInitialiseTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonSolutionInitialiseTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    SELECT CASE(solutionInitialiseType)
    CASE(SOLVER_SOLUTION_INITIALISE_ZERO)
      newtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_ZERO
    CASE(SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD)
      newtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_CURRENT_FIELD
    CASE(SOLVER_SOLUTION_INITIALISE_NO_CHANGE)
      newtonSolver%solutionInitialiseType=SOLVER_SOLUTION_INITIALISE_NO_CHANGE
    CASE DEFAULT
      localError="The specified solution initialise type  of "// &
        & TRIM(NumberToVString(solutionInitialiseType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Solver_NewtonSolutionInitialiseTypeSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonSolutionInitialiseTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonSolutionInitialiseTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the solution tolerance for a nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonSolutionToleranceSet
  SUBROUTINE Solver_NewtonSolutionToleranceSet(solver,solutionTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the solution tolerance for
    REAL(DP), INTENT(IN) :: solutionTolerance !<The solution tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonSolutionToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(solutionTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified solution tolerance of "//TRIM(NumberToVString(solutionTolerance,"*",err,error))// &
        & " is invalid. The relative tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    newtonSolver%solutionTolerance=solutionTolerance
    
    EXITS("Solver_NewtonSolutionToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonSolutionToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonSolutionToleranceSet
        
  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton solver 
  SUBROUTINE SolverNonlinearNewton_Solve(newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer to the nonlinear Newton solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewton_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*999)
    
    SELECT CASE(newtonSolver%newtonSolveType)
    CASE(SOLVER_NEWTON_LINESEARCH)
      CALL SolverNonlinearNewtonLinesearch_Solve(newtonSolver%linesearchSolver,err,error,*999)
    CASE(SOLVER_NEWTON_TRUSTREGION)
      CALL SolverNonlinearNewtonTrustregion_Solve(newtonSolver%trustregionSolver,err,error,*999)
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearNewton_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_Solve
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating nonlinear Newton trust region solver
  SUBROUTINE SolverNonlinearNewtonTrustregion_CreateFinish(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer the nonlinear Newton trust region solver to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    EXTERNAL :: Problem_SolverResidualEvaluatePetsc
    INTEGER(INTG) :: equationsMatrixIdx,equationsSetIdx,groupCommunicator
    TYPE(DistributedVectorType), POINTER :: residualVector
    TYPE(DistributedVectorPETScType), POINTER :: residualPETScVector
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: linearVariable
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(WorkGroupType), POINTER :: workGroup
    TYPE(VARYING_STRING) :: localError
  
    ENTERS("SolverNonlinearNewtonTrustregion_CreateFinish",err,error,*999)
    
    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trust region solver is not associated.",err,error,*999)
    
    NULLIFY(newtonSolver)
    CALL SolverNonlinearNewtonTrustregion_NewtonSolverGet(trustregionSolver,newtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearNewton_NonlinearSolverGet(newtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(workGroup)
    CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    SELECT CASE(trustregionSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Loop over the equations set in the solver equations
      DO equationsSetIdx=1,solverMapping%numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          !If there are any linear matrices create temporary vector for matrix-vector products
          NULLIFY(vectorMatrices)
          CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          DO equationsMatrixIdx=1,linearMatrices%numberOfLinearMatrices
            NULLIFY(equationsMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
            IF(.NOT.ASSOCIATED(equationsMatrix%tempVector)) THEN
              NULLIFY(linearVariable)
              CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,linearVariable,err,error,*999)
              NULLIFY(domainMapping)
              CALL FieldVariable_DomainMappingGet(linearVariable,domainMapping,err,error,*999)
              CALL DistributedVector_CreateStart(domainMapping,equationsMatrix%tempVector,err,error,*999)
              CALL DistributedVector_DataTypeSet(equationsMatrix%tempVector,DISTRIBUTED_MATRIX_VECTOR_DP_TYPE,err,error,*999)
              CALL DistributedVector_CreateFinish(equationsMatrix%tempVector,err,error,*999)
            ENDIF
          ENDDO !equationsMatrixIdx
        ENDIF
      ENDDO !equationsSetIdx                 
      !Create the solver matrices and vectors
      NULLIFY(solverMatrices)
      CALL SolverMatrices_CreateStart(solverEquations,solverMatrices,err,error,*999)
      CALL SolverMatrices_LibraryTypeSet(solverMatrices,SOLVER_PETSC_LIBRARY,err,error,*999)
!!TODO: set up the matrix structure if using an analytic Jacobian
      SELECT CASE(solverEquations%symmetryType)
      CASE(SOLVER_SYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_SYMMETRIC_TYPE,err,error,*999)
      CASE(SOLVER_UNSYMMETRIC_MATRICES)
        CALL SolverMatrices_SymmetryTypesSet(solverMatrices,DISTRIBUTED_MATRIX_UNSYMMETRIC_TYPE,err,error,*999)
      CASE DEFAULT
        localError="The specified solver equations symmetry type of "// &
          & TRIM(NumberToVString(solverEquations%symmetryType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL SolverMatrices_CreateFinish(solverMatrices,err,error,*999)
      !Create the PETSc SNES solver
      CALL PETSc_SNESCreate(groupCommunicator,trustregionSolver%snes,err,error,*999)
      !Set the nonlinear solver type to be a Newton trust region solver
      CALL PETSc_SNESSetType(trustregionSolver%snes,PETSC_SNESNEWTONTR,err,error,*999)
      !Set the solver as the SNES application context
      CALL PETSc_SNESSetApplicationContext(trustregionSolver%snes, &
        & trustregionSolver%newtonSolver%nonlinearSolver%solver,err,error,*999)
      !Set the nonlinear function
      NULLIFY(residualVector)
      CALL SolverMatrices_ResidualDistributedVectorGet(solverMatrices,residualVector,err,error,*999)
      NULLIFY(residualPETScVector)
      CALL DistributedVector_PETScVectorGet(residualVector,residualPETScVector,err,error,*999)
      CALL PETSc_SNESSetFunction(trustregionSolver%snes,residualPETScVector%vector, &
        & Problem_SolverResidualEvaluatePetsc,trustregionSolver%newtonSolver%nonlinearSolver%solver,err,error,*999)
      !Set the Jacobian if necessary
      !Set the trust region delta ???
      
      !Set the trust region tolerance
      CALL PETSc_SNESSetTrustRegionTolerance(trustregionSolver%snes,trustregionSolver%trustregionTolerance, &
        & err,error,*999)
      !Set the tolerances for the SNES solver
      CALL PETSc_SNESSetTolerances(trustregionSolver%snes,newtonSolver%absoluteTolerance, &
        & newtonSolver%relativeTolerance,newtonSolver%solutionTolerance, &
        & newtonSolver%maximumNumberOfIterations,newtonSolver%maximumNumberOfFunctionEvaluations, &
        & err,error,*999)
      !Set any further SNES options from the command line options
      CALL PETSc_SNESSetFromOptions(trustregionSolver%snes,err,error,*999)
    CASE DEFAULT
      localError="The solver library type of "// &
        & TRIM(NumberToVString(trustregionSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinearNewtonTrustregion_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewtonTrustregion_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewtonTrustregion_CreateFinish

  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region delta0 for a nonlinear Newton trust region solver solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonTrustRegionDelta0Set
  SUBROUTINE Solver_NewtonTrustregionDelta0Set(solver,trustregionDelta0,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the trust region delta0 for
    REAL(DP), INTENT(IN) :: trustregionDelta0 !<The trust region delta0 to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonTrustregionDelta0Set",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsTrustregion(newtonSolver,err,error,*999)
    NULLIFY(trustregionSolver)
    CALL SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,trustregionSolver,err,error,*999)
    IF(trustregionDelta0<=ZERO_TOLERANCE) THEN
      localError="The specified trust region delta0 of "// &
        & TRIM(NumberToVString(trustregionDelta0,"*",err,error))// &
        & " is invalid. The trust region delta0 must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    trustregionSolver%trustregionDelta0=trustregionDelta0
    
    EXITS("Solver_NewtonTrustregionDelta0Set")
    RETURN
999 ERRORSEXITS("Solver_NewtonTrustregionDelta0Set",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NewtonTrustregionDelta0Set
        
  !
  !================================================================================================================================
  !
  
  !>Finalise a nonlinear Newton trust region solver and deallocate all memory
  SUBROUTINE SolverNonlinearNewton_TrustregionFinalise(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer the non linear trust region solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
  
    ENTERS("SolverNonlinearNewton_TrustregionFinalise",err,error,*999)

    IF(ASSOCIATED(trustregionSolver)) THEN      
      CALL PETSc_SNESFinalise(trustregionSolver%snes,err,error,*999)
      DEALLOCATE(trustregionSolver)
    ENDIF
    
    EXITS("SolverNonlinearNewton_TrustregionFinalise")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewton_TrustregionFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewton_TrustregionFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a Newton trust region solver for a nonlinear solver
  SUBROUTINE SolverNonlinearNewton_TrustregionInitialise(newtonSolver,err,error,*)

    !Argument variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver !<A pointer the Newton solver to initialise the trust region solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
  
    ENTERS("SolverNonlinearNewton_TrustregionInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(newtonSolver)) CALL FlagError("Newton solver is not associated.",err,error,*998)
    IF(ASSOCIATED(newtonSolver%trustregionSolver)) &
      & CALL FlagError("Trust region solver is already associated for this nonlinear solver.",err,error,*998)
    
    ALLOCATE(newtonSolver%trustregionSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate Newton solver trust region solver.",err,error,*999)
    newtonSolver%trustregionSolver%newtonSolver=>newtonSolver
    newtonSolver%trustregionSolver%solverLibrary=SOLVER_PETSC_LIBRARY
!!TODO: set this properly
    newtonSolver%trustregionSolver%trustregionDelta0=0.01_DP
    CALL PETSc_SNESInitialise(newtonSolver%trustregionSolver%snes,err,error,*999)
        
    EXITS("SolverNonlinearNewton_TrustregionInitialise")
    RETURN
999 CALL SolverNonlinearNewton_TrustregionFinalise(newtonSolver%trustregionSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("SolverNonlinearNewton_TrustregionInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinearNewton_TrustregionInitialise

  !
  !================================================================================================================================
  !

  !Solves a nonlinear Newton trust region solver 
  SUBROUTINE SolverNonlinearNewtonTrustregion_Solve(trustregionSolver,err,error,*)

    !Argument variables
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver !<A pointer to the nonlinear Newton trust region solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinearNewtonTrustregion_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(trustregionSolver)) CALL FlagError("Trustregion solver is not associated.",err,error,*999)
    
    NULLIFY(newtonSolver)
    CALL SolverNonlinearNewtonTrustregion_NewtonSolverGet(trustregionSolver,newtonSolver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL SolverNonlinearNewton_NonlinearSolverGet(newtonSolver,nonlinearSolver,err,error,*999)
    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    SELECT CASE(trustregionSolver%solverLibrary)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The nonlinear Newton trust region solver library type of "// &
        & TRIM(NumberToVString(trustregionSolver%solverLibrary,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
         
    EXITS("SolverNonlinearNewtonTrustregion_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinearNewtonTrustregion_Solve",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinearNewtonTrustregion_Solve
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the trust region tolerance for a nonlinear Newton trust region solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonTrustRegionToleranceSet
  SUBROUTINE Solver_NewtonTrustregionToleranceSet(solver,trustregionTolerance,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the trust region tolerance for
    REAL(DP), INTENT(IN) :: trustregionTolerance !<The trust region tolerance to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonTrustregionSolverType), POINTER :: trustregionSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_NewtonTrustregionToleranceSet",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    CALL SolverNonlinearNewton_AssertIsTrustregion(newtonSolver,err,error,*999)
    NULLIFY(trustregionSolver)
    CALL SolverNonlinearNewton_TrustregionSolverGet(newtonSolver,trustregionSolver,err,error,*999)
    IF(trustregionTolerance<=ZERO_TOLERANCE) THEN
      localError="The specified trust region tolerance of "// &
        & TRIM(NumberToVString(trustregionTolerance,"*",err,error))// &
        & " is invalid. The trust region tolerance must be > 0.0."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    trustregionSolver%trustregionTolerance=trustregionTolerance
    
    EXITS("Solver_NewtonTrustregionToleranceSet")
    RETURN
999 ERRORSEXITS("Solver_NewtonTrustregionToleranceSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonTrustregionToleranceSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of nonlinear Newton solver. \see OpenCMISS::Iron::cmfe_Solver_NewtonTypeSet
  SUBROUTINE Solver_NewtonTypeSet(solver,newtonSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear Newton solver type
    INTEGER(INTG), INTENT(IN) :: newtonSolveType !<The type of nonlinear solver to set \see SolverRoutines_NewtonSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solver_NewtonTypeSet",err,error,*998)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    CALL SolverNonlinear_AssertIsNewton(nonlinearSolver,err,error,*999)
    NULLIFY(newtonSolver)
    CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
    IF(newtonSolveType/=newtonSolver%newtonSolveType) THEN
      !Intialise the new solver type
      SELECT CASE(newtonSolveType)
      CASE(SOLVER_NEWTON_LINESEARCH)
        CALL SolverNonlinearNewton_LinesearchInitialise(newtonSolver,err,error,*999)
      CASE(SOLVER_NEWTON_TRUSTREGION)
        CALL SolverNonlinearNewton_TrustregionInitialise(newtonSolver,err,error,*999)
      CASE DEFAULT
        localError="The Newton solver type of "//TRIM(NumberToVString(newtonSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old solver type
      SELECT CASE(newtonSolver%newtonSolveType)
      CASE(SOLVER_NEWTON_LINESEARCH)
        CALL SolverNonlinearNewton_LinesearchFinalise(newtonSolver%linesearchSolver,err,error,*999)
      CASE(SOLVER_NEWTON_TRUSTREGION)
        CALL SolverNonlinearNewton_TrustregionFinalise(newtonSolver%trustregionSolver,err,error,*999)
      CASE DEFAULT
        localError="The Newton solver type of "// &
          & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      newtonSolver%newtonSolveType=newtonSolveType
    ENDIF
    
    EXITS("Solver_NewtonTypeSet")
    RETURN
999 SELECT CASE(newtonSolveType)
    CASE(SOLVER_NEWTON_LINESEARCH)
      CALL SolverNonlinearNewton_LinesearchFinalise(newtonSolver%linesearchSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_NEWTON_TRUSTREGION)
      CALL SolverNonlinearNewton_TrustregionFinalise(newtonSolver%trustregionSolver,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("Solver_NewtonTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NewtonTypeSet
        
  !
  !================================================================================================================================
  !

  !>Finishes the process of creating a nonlinear solver 
  SUBROUTINE SolverNonlinear_CreateFinish(nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinear_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
    
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      CALL SolverNonlinearNewton_CreateFinish(nonlinearSolver%newtonSolver,err,error,*999)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      CALL SolverNonlinearQuasiNewton_CreateFinish(nonlinearSolver%quasiNewtonSolver,err,error,*999)
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverNonlinear_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverNonlinear_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_CreateFinish

  !
  !================================================================================================================================
  !

  !>Instead of warning on nonlinear divergence, exit with error
  SUBROUTINE Solver_NonlinearDivergenceExit(solver,err,error,*)
    
    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local variables
    INTEGER(INTG) :: convergedReason
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonLinesearchSolverType), POINTER :: newtonLinesearchSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: quasiNewtonLinesearchSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_NonlinearDivergenceExit",err,error,*999)

    CALL Solver_AssertIsNonlinear(solver,err,error,*999)    
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      SELECT CASE (newtonSolver%newtonSolveType)
      CASE(SOLVER_NEWTON_LINESEARCH)
        NULLIFY(newtonLinesearchSolver)
        CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,newtonLinesearchSolver,err,error,*999)
        CALL PETSc_SNESGetConvergedReason(newtonLinesearchSolver%snes,convergedReason,err,error,*999)
        SELECT CASE(convergedReason)
        CASE(PETSC_SNES_DIVERGED_FUNCTION_COUNT)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged function count.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LINEAR_SOLVE)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged linear solve.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_FNORM_NAN)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged F Norm NaN.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_MAX_IT)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged maximum iterations.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LINE_SEARCH)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged line search.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LOCAL_MIN)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged local minimum.", &
            & err,error,*999)
        END SELECT
      CASE(SOLVER_NEWTON_TRUSTREGION)
        !Not yet implemented. Don't kick up a fuss, just exit
      END SELECT
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      !Not yet implemented. Don't kick up a fuss, just exit
    CASE(SOLVER_NONLINEAR_SQP)
      !Not yet implemented. Don't kick up a fuss, just exit
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      SELECT CASE (quasiNewtonSolver%quasiNewtonSolveType)
      CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
        NULLIFY(quasiNewtonLinesearchSolver)
        CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,quasiNewtonLinesearchSolver,err,error,*999)
        CALL PETSc_SNESGetConvergedReason(quasiNewtonLinesearchSolver%snes,convergedReason,err,error,*999)
        SELECT CASE(convergedReason)
        CASE(PETSC_SNES_DIVERGED_FUNCTION_COUNT)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged function count.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LINEAR_SOLVE)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged linear solve.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_FNORM_NAN)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged F Norm NaN.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_MAX_IT)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged maximum iterations.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LINE_SEARCH)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged line search.", &
            & err,error,*999)
        CASE(PETSC_SNES_DIVERGED_LOCAL_MIN)
          CALL FlagError("Nonlinear line search solver did not converge. Exit due to PETSc diverged local minimum.", &
            & err,error,*999)
        END SELECT
      CASE(SOLVER_NEWTON_TRUSTREGION)
        !Not yet implemented. Don't kick up a fuss, just exit
      END SELECT
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("Solver_NonlinearDivergenceExit")
    RETURN
999 ERRORSEXITS("Solver_NonlinearDivergenceExit",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_NonlinearDivergenceExit

  !
  !================================================================================================================================
  !

  !>Finalise a nonlinear solver for a solver.
  RECURSIVE SUBROUTINE Solver_NonlinearFinalise(nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer the nonlinear solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_NonlinearFinalise",err,error,*999)

    IF(ASSOCIATED(nonlinearSolver)) THEN
      SELECT CASE(nonlinearSolver%nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        CALL SolverNonlinear_NewtonFinalise(nonlinearSolver%newtonSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_NONLINEAR_SQP)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        CALL SolverNonlinear_QuasiNewtonFinalise(nonlinearSolver%quasiNewtonSolver,err,error,*999)
      CASE DEFAULT
        localError="The nonlinear solver type of "// &
          & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      DEALLOCATE(nonlinearSolver)
    ENDIF
         
    EXITS("Solver_NonlinearFinalise")
    RETURN
999 ERRORSEXITS("Solver_NonlinearFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NonlinearFinalise

  !
  !================================================================================================================================
  !

  !>Initialise a nonlinear solver for a solver.
  SUBROUTINE Solver_NonlinearInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the nonlinear solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("Solver_NonlinearInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%nonlinearSolver)) &
        & CALL FlagError("Nonlinear solver is already associated for this solver.",err,error,*998)
      
    !Allocate and initialise a Nonlinear solver
    ALLOCATE(solver%nonlinearSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver nonlinear solver.",err,error,*999)
    solver%nonlinearSolver%solver=>solver
    NULLIFY(solver%nonlinearSolver%newtonSolver)
    !Default to a nonlinear Newton solver
    solver%nonlinearSolver%nonlinearSolveType=SOLVER_NONLINEAR_NEWTON
    CALL SolverNonlinear_NewtonInitialise(solver%nonlinearSolver,err,error,*999)
        
    EXITS("Solver_NonlinearInitialise")
    RETURN
999 CALL Solver_NonlinearFinalise(solver%nonlinearSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_NonlinearInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NonlinearInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for a nonlinear solver.
  SUBROUTINE SolverNonlinear_LibraryTypeSet(nonlinearSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer the nonlinear solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the nonlinear solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("SolverNonlinear_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
    
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      CALL SolverNonlinearNewton_LibraryTypeSet(newtonSolver,solverLibraryType,err,error,*999)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      CALL SolverNonlinearQuasiNewton_LibraryTypeSet(quasiNewtonSolver,solverLibraryType,err,error,*999)
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("SolverNonlinear_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverNonlinear_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Monitors the nonlinear solve.
  SUBROUTINE SolverNonlinear_Monitor(nonlinearSolver,its,norm,err,error,*)

   !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to monitor
    INTEGER(INTG), INTENT(IN) :: its !<The number of iterations
    REAL(DP), INTENT(IN) :: norm !<The residual norm
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(DP) :: xnorm !<The norm of the current solution 
    REAL(DP) :: fnorm !<The norm of the current function 
    REAL(DP) :: ynorm !<The norm of the current update
    TYPE(NewtonLinesearchSolverType), POINTER :: newtonLinesearchSolver
    TYPE(NewtonSolverType), POINTER :: newtonSolver
    TYPE(NewtonSolverConvergenceTestType), POINTER :: convergenceTest
    TYPE(QuasiNewtonLinesearchSolverType), POINTER :: quasiNewtonlinesearchSolver
    TYPE(QuasiNewtonSolverType), POINTER :: quasiNewtonSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("SolverNonlinear_Monitor",err,error,*999)

    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)
        
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"Nonlinear solve monitor: ",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,ERROR,*999)
    CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Iteration number = ",its,err,error,*999)
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      NULLIFY(newtonSolver)
      CALL SolverNonlinear_NewtonSolverGet(nonlinearSolver,newtonSolver,err,error,*999)
      SELECT CASE(newtonSolver%convergenceTestType)
      CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Function Norm = ",norm,err,error,*999)
      CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
        NULLIFY(convergenceTest)
        CALL SolverNonlinearNewton_ConvergenceTestGet(newtonSolver,convergenceTest,err,error,*999)
        SELECT CASE(newtonSolver%newtonSolveType)
        CASE(SOLVER_NEWTON_LINESEARCH)
          NULLIFY(newtonLinesearchSolver)
          CALL SolverNonlinearNewton_LinesearchSolverGet(newtonSolver,newtonLinesearchSolver,err,error,*999)
          CALL PETSc_SNESLineSearchGetNorms(newtonLinesearchSolver%sneslinesearch,xnorm,fnorm,ynorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Solution Norm          = ",xnorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Solution Update Norm   = ",ynorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Function Norm          = ",fnorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Normalised Energy Norm = ",convergenceTest%normalisedEnergy,err,error,*999)
        CASE(SOLVER_NEWTON_TRUSTREGION)
          CALL FlagError("The Newton Trust region solver is not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The Newton solve type of "// &
            & TRIM(NumberToVString(newtonSolver%newtonSolveType,"*",err,error))//"is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
        CALL FlagError("The Sum of differentiated ratios of unconstrained to constrained residuals"// &
          &  "convergence test type is not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The convergence test type of "//TRIM(NumberToVString(newtonSolver%convergenceTestType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      CALL WriteString(GENERAL_OUTPUT_TYPE,"  Newton solver information: ",err,error,*999)          
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"    Number of function evaluations = ",newtonSolver% &
        & totalNumberOfFunctionEvaluations,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"    Number of Jacobian evaluations = ",newtonSolver% &
        & totalNumberOfJacobianEvaluations,err,error,*999)            
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      !Do nothing
    CASE(SOLVER_NONLINEAR_SQP)
      !Do nothing
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      NULLIFY(quasiNewtonSolver)
      CALL SolverNonlinear_QuasiNewtonSolverGet(nonlinearSolver,quasiNewtonSolver,err,error,*999)
      SELECT CASE(quasiNewtonSolver%convergenceTestType)
      CASE(SOLVER_NEWTON_CONVERGENCE_PETSC_DEFAULT)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Function Norm    = ",norm,err,error,*999)
      CASE(SOLVER_NEWTON_CONVERGENCE_ENERGY_NORM)
        CALL SolverNonlinearQuasiNewton_ConvergenceTestGet(quasiNewtonSolver,convergenceTest,err,error,*999)
        SELECT CASE(quasiNewtonSolver%quasiNewtonSolveType)
        CASE(SOLVER_QUASI_NEWTON_LINESEARCH)
          NULLIFY(quasiNewtonLinesearchSolver)
          CALL SolverNonlinearQuasiNewton_LinesearchSolverGet(quasiNewtonSolver,quasiNewtonLinesearchSolver,err,error,*999)
          NULLIFY(convergenceTest)
          CALL PETSc_SNESLineSearchGetNorms(quasiNewtonLinesearchSolver%sneslinesearch,xnorm,fnorm,ynorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Solution Norm          = ",xnorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Solution Update Norm   = ",ynorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Function Norm          = ",fnorm,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Normalised Energy Norm = ",convergenceTest%normalisedEnergy,err,error,*999)
        CASE(SOLVER_NEWTON_TRUSTREGION)
          CALL FlagError("The Quasi-Newton trustregion solver is not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The Quasi-Newton solve type of "// &
            & TRIM(NumberToVString(quasiNewtonSolver%quasiNewtonSolveType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(SOLVER_NEWTON_CONVERGENCE_DIFFERENTIATED_RATIO)
        CALL FlagError("The Sum of differentiated ratios of unconstrained to constrained residuals "// &
          &  "convergence test type is not implemented.",err,error,*999)
      END SELECT
      CALL WriteString(GENERAL_OUTPUT_TYPE,"  Quasi-Newton solver information: ",err,error,*999)          
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"    Number of function evaluations = ",quasiNewtonSolver% &
        & totalNumberOfFunctionEvaluations,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"    Number of Jacobian evaluations = ",quasiNewtonSolver% &
        & totalNumberOfJacobianEvaluations,err,error,*999)            
    CASE DEFAULT
      localError="The nonlinear solver type of "// &
        & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("SolverNonlinear_Monitor")
    RETURN
999 ERRORSEXITS("SolverNonlinear_Monitor",err,error)
    RETURN 1
    
  END SUBROUTINE SolverNonlinear_Monitor

  !
  !================================================================================================================================
  !

  !Solves a nonlinear solver 
  SUBROUTINE SolverNonlinear_Solve(nonlinearSolver,err,error,*)

    !Argument variables
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver !<A pointer to the nonlinear solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfSolverMatrices,solverMatrixIdx
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverNonlinear_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(nonlinearSolver)) CALL FlagError("Nonlinear solver is not associated.",err,error,*999)

    NULLIFY(solver)
    CALL SolverNonlinear_SolverGet(nonlinearSolver,solver,err,error,*999)
    SELECT CASE(nonlinearSolver%nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      CALL SolverNonlinearNewton_Solve(nonlinearSolver%newtonSolver,err,error,*999)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      CALL SolverNonlinearQuasiNewton_Solve(nonlinearSolver%quasiNewtonSolver,err,error,*999)
    CASE DEFAULT
      localError="The nonlinear solver type of "//TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    IF(solver%outputType>=SOLVER_SOLVER_OUTPUT) THEN
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Solver solution vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of solution vectors = ",numberOfSolverMatrices,err,error,*999)
      DO solverMatrixIdx=1,numberOfSolverMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        NULLIFY(solverVector)
        CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Solution vector for solver matrix : ",solverMatrixIdx,err,error,*999)
        CALL DistributedVector_Output(GENERAL_OUTPUT_TYPE,solverVector,err,error,*999)
      ENDDO !solverMatrixIdx
    ENDIF
        
    EXITS("SolverNonlinear_Solve")
    RETURN
999 ERRORSEXITS("SolverNonlinear_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverNonlinear_Solve
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of nonlinear solver. \see OpenCMISS::Iron::cmfe_Solver_NonlinearTypeSet
  SUBROUTINE Solver_NonlinearTypeSet(solver,nonlinearSolveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the nonlinear solver type
    INTEGER(INTG), INTENT(IN) :: nonlinearSolveType !<The type of nonlinear solver to set \see SolverRoutines_NonlinearSolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solver_NonlinearTypeSet",err,error,*998)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsNonlinear(solver,err,error,*999)
    NULLIFY(nonlinearSolver)
    CALL Solver_NonlinearSolverGet(solver,nonlinearSolver,err,error,*999)
    IF(nonlinearSolveType/=nonlinearSolver%nonlinearSolveType) THEN
      CALL Solver_LinkedSolverRemove(solver,SOLVER_LINEAR_TYPE,err,error,*999)
      !Finalise the old solver type
      SELECT CASE(nonlinearSolver%nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        CALL SolverNonlinear_NewtonFinalise(nonlinearSolver%newtonSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
        CALL FlagError("Not implemented.",err,error,*999)                
      CASE(SOLVER_NONLINEAR_SQP)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        CALL SolverNonlinear_QuasiNewtonFinalise(nonlinearSolver%quasiNewtonSolver,err,error,*999)
      CASE DEFAULT
        localError="The nonlinear solver type of "// &
          & TRIM(NumberToVString(nonlinearSolver%nonlinearSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      nonlinearSolver%nonlinearSolveType=nonlinearSolveType
      !Intialise the new solver type
      SELECT CASE(nonlinearSolveType)
      CASE(SOLVER_NONLINEAR_NEWTON)
        NULLIFY(nonlinearSolver%newtonSolver)
        CALL SolverNonlinear_NewtonInitialise(nonlinearSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_NONLINEAR_SQP)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
        CALL SolverNonlinear_QuasiNewtonInitialise(nonlinearSolver,err,error,*999)
      CASE DEFAULT
        localError="The specified nonlinear solver type of "// &
          & TRIM(NumberToVString(nonlinearSolveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_NonlinearTypeSet")
    RETURN
999 SELECT CASE(nonlinearSolveType)
    CASE(SOLVER_NONLINEAR_NEWTON)
      CALL SolverNonlinear_NewtonFinalise(nonlinearSolver%newtonSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_NONLINEAR_BFGS_INVERSE)
      CALL FlagError("Not implemented.",err,error,*998)                
    CASE(SOLVER_NONLINEAR_SQP)
      CALL FlagError("Not implemented.",err,error,*998)      
    CASE(SOLVER_NONLINEAR_QUASI_NEWTON)
      CALL SolverNonlinear_QuasiNewtonFinalise(nonlinearSolver%quasiNewtonSolver,dummyErr,dummyError,*998)
    END SELECT
998 ERRORSEXITS("Solver_NonlinearTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_NonlinearTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of certainty for an optimisation solver.
  SUBROUTINE SolverOptimiser_CertaintyTypeSet(optimiserSolver,solverCertaintyType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to set the certainty type for.
    INTEGER(INTG), INTENT(IN) :: solverCertaintyType !<The type of certainty for the optimiser solver to set. \see SolverRoutines_OptimiserCertaintyType,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverOptimiser_CertaintyTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    SELECT CASE(solverCertaintyType)
    CASE(SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY)
      optimiserSolver%certaintyType=SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY
    CASE(SOLVER_OPTIMISER_STOCHASTIC_CERTAINTY)
      optimiserSolver%certaintyType=SOLVER_OPTIMISER_STOCHASTIC_CERTAINTY
    CASE DEFAULT
      localError="The specified certainty type of "//TRIM(NumberToVString(solverCertaintyType,"*",err,error))// &
        & " is invalid for an optimiser solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverOptimiser_CertaintyTypeSet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_CertaintyTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_CertaintyTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of constraint for an optimisation solver.
  SUBROUTINE SolverOptimiser_ConstraintTypeSet(optimiserSolver,solverConstraintType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to set the constraint type for.
    INTEGER(INTG), INTENT(IN) :: solverConstraintType !<The type of constraint for the optimiser solver to set. \see SolverRoutines_OptimiserConstraintType,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverOptimiser_ConstraintTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    SELECT CASE(solverConstraintType)
    CASE(SOLVER_OPTIMISER_UNCONSTRAINED)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_OPTIMISER_BOUND_CONSTRAINED)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_OPTIMISER_LINEAR_CONSTRAINTS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_OPTIMISER_NONLINEAR_CONSTRAINTS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_OPTIMISER_PDE_CONSTRAINTS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The specified constraint type of "//TRIM(NumberToVString(solverConstraintType,"*",err,error))// &
        & " is invalid for an optimiser solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverOptimiser_ConstraintTypeSet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_ConstraintTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_ConstraintTypeSet

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an optimiser solver 
  SUBROUTINE SolverOptimiser_CreateFinish(optimiserSolver,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer to the optimiser solver to finish the creation of.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverOptimiser_CreateFinish",err,error,*999)

    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverOptimiser_CreateFinish")
    RETURN
999 ERRORSEXITS("SolverOptimiser_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Finalise a optimiser solver.
  SUBROUTINE Solver_OptimiserFinalise(optimiserSolver,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solver_OptimiserFinalise",err,error,*999)

    IF(ASSOCIATED(optimiserSolver)) THEN        
      DEALLOCATE(optimiserSolver)
    ENDIF
         
    EXITS("Solver_OptimiserFinalise")
    RETURN
999 ERRORSEXITS("Solver_OptimiserFinalise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_OptimiserFinalise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of gradient calculation type for an optimisation solver.
  SUBROUTINE Solver_OptimiserGradientCalculationTypeSet(solver,gradientCalculationType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the gradient calculation type
    INTEGER(INTG), INTENT(IN) :: gradientCalculationType !<The type of gradient calculation type to set for an optimisation solver \see SolverRoutines_OptimiserGradientCalculationTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_OptimiserGradientCalculationTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsOptimiser(solver,err,error,*999)
    NULLIFY(optimiserSolver)
    CALL Solver_OptimiserSolverGet(solver,optimiserSolver,err,error,*999)
    IF(gradientCalculationType/=optimiserSolver%gradientCalculationType) THEN
      SELECT CASE(gradientCalculationType)
      CASE(SOLVER_OPTIMISER_GRADIENT_NOT_CALCULATED)
        optimiserSolver%gradientCalculationType=SOLVER_OPTIMISER_GRADIENT_NOT_CALCULATED
      CASE(SOLVER_OPTIMISER_GRADIENT_EQUATIONS_CALCULATED)
        optimiserSolver%gradientCalculationType=SOLVER_OPTIMISER_GRADIENT_EQUATIONS_CALCULATED
      CASE(SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED)
        optimiserSolver%gradientCalculationType=SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED
      CASE DEFAULT
        localError="The gradient calculation type of "//TRIM(NumberToVString(gradientCalculationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_OptimiserGradientCalculationTypeSet")
    RETURN
999 ERRORSEXITS("Solver_OptimiserGradientCalculationTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OptimiserGradientCalculationTypeSet
        
  !
  !================================================================================================================================
  !

  !>Sets/changes the type of Hessian calculation type for an optimisation solver.
  SUBROUTINE Solver_OptimiserHessianCalculationTypeSet(solver,hessianCalculationType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the Hessian calculation type
    INTEGER(INTG), INTENT(IN) :: hessianCalculationType !<The type of Hessian calculation type to set for an optimisation solver \see SolverRoutines_OptimiserHessianCalculationTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_OptimiserHessianCalculationTypeSet",err,error,*999)

    CALL Solver_AssertNotFinished(solver,err,error,*999)
    CALL Solver_AssertIsOptimiser(solver,err,error,*999)
    NULLIFY(optimiserSolver)
    CALL Solver_OptimiserSolverGet(solver,optimiserSolver,err,error,*999)
    IF(hessianCalculationType/=optimiserSolver%hessianCalculationType) THEN
      SELECT CASE(hessianCalculationType)
      CASE(SOLVER_OPTIMISER_HESSIAN_NOT_CALCULATED)
        optimiserSolver%hessianCalculationType=SOLVER_OPTIMISER_HESSIAN_NOT_CALCULATED
      CASE(SOLVER_OPTIMISER_HESSIAN_EQUATIONS_CALCULATED)
        optimiserSolver%hessianCalculationType=SOLVER_OPTIMISER_HESSIAN_EQUATIONS_CALCULATED
      CASE(SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED)
        optimiserSolver%hessianCalculationType=SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED
      CASE DEFAULT
        localError="The Hessian calculation type of "//TRIM(NumberToVString(hessianCalculationType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF
    
    EXITS("Solver_OptimiserHessianCalculationTypeSet")
    RETURN
999 ERRORSEXITS("Solver_OptimiserHessianCalculationTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_OptimiserHessianCalculationTypeSet
        
  !
  !================================================================================================================================
  !

  !>Initialise an optimiser solver for a solver.
  SUBROUTINE Solver_OptimiserInitialise(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to initialise the optimiser solver for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solver_OptimiserInitialise",err,error,*998)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*998)
    IF(ASSOCIATED(solver%optimiserSolver)) CALL FlagError("Optimiser solver is already associated for this solver.",err,error,*998)
     
    ALLOCATE(solver%optimiserSolver,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solver optimiser solver.",err,error,*999)
    solver%optimiserSolver%solver=>solver
    solver%optimiserSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    solver%optimiserSolver%variableType=SOLVER_OPTIMISER_CONTINUOUS_VARIABLES
    solver%optimiserSolver%objectiveType=SOLVER_OPTIMISER_ONE_OBJECTIVE
    solver%optimiserSolver%constraintType=SOLVER_OPTIMISER_UNCONSTRAINED
    solver%optimiserSolver%certaintyType=SOLVER_OPTIMISER_DETERMINISTIC_CERTAINTY
    solver%optimiserSolver%gradientCalculationType=SOLVER_OPTIMISER_GRADIENT_FD_CALCULATED
    solver%optimiserSolver%hessianCalculationType=SOLVER_OPTIMISER_HESSIAN_FD_CALCULATED
        
    EXITS("Solver_OptimiserInitialise")
    RETURN
999 CALL Solver_OptimiserFinalise(solver%optimiserSolver,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_OptimiserInitialise",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_OptimiserInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of library to use for an optimisation solver.
  SUBROUTINE SolverOptimiser_LibraryTypeSet(optimiserSolver,solverLibraryType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to get the library type for.
    INTEGER(INTG), INTENT(IN) :: solverLibraryType !<The type of library for the optimiser solver to set. \see SolverRoutines_SolverLibraries,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverOptimiser_LibraryTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    SELECT CASE(solverLibraryType)
    CASE(SOLVER_CMISS_LIBRARY)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(SOLVER_PETSC_LIBRARY)
      optimiserSolver%solverLibrary=SOLVER_PETSC_LIBRARY
    CASE DEFAULT
      localError="The specified solver library type of "//TRIM(NumberToVString(solverLibraryType,"*",err,error))// &
        & " is invalid for an optimiser solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverOptimiser_LibraryTypeSet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_LibraryTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_LibraryTypeSet

  !
  !================================================================================================================================
  !

  !>Monitors the optimiser solve.
  SUBROUTINE SolverOptimiser_Monitor(optimiserSolver,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer to the optimiser solver to monitor
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    
    ENTERS("SolverOptimiser_Monitor",err,error,*999)

    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
        
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"Optimiser solve monitor: ",err,error,*999)
    CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
     
    EXITS("SolverOptimiser_Monitor")
    RETURN
999 ERRORSEXITS("SolverOptimiser_Monitor",err,error)
    RETURN 1
    
  END SUBROUTINE SolverOptimiser_Monitor

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of objective for an optimisation solver.
  SUBROUTINE SolverOptimiser_ObjectiveTypeSet(optimiserSolver,solverObjectiveType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to set the objective type for.
    INTEGER(INTG), INTENT(IN) :: solverObjectiveType !<The type of objective for the optimiser solver to set. \see SolverRoutines_OptimiserObjectiveType,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverOptimiser_ObjectiveTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    SELECT CASE(solverObjectiveType)
    CASE(SOLVER_OPTIMISER_NO_OBJECTIVE)
      optimiserSolver%objectiveType=SOLVER_OPTIMISER_NO_OBJECTIVE
    CASE(SOLVER_OPTIMISER_ONE_OBJECTIVE)
      optimiserSolver%objectiveType=SOLVER_OPTIMISER_ONE_OBJECTIVE
    CASE(SOLVER_OPTIMISER_MANY_OBJECTIVE)
      optimiserSolver%objectiveType=SOLVER_OPTIMISER_MANY_OBJECTIVE
    CASE DEFAULT
      localError="The specified objective type of "//TRIM(NumberToVString(solverObjectiveType,"*",err,error))// &
        & " is invalid for an optimiser solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverOptimiser_ObjectiveTypeSet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_ObjectiveTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_ObjectiveTypeSet

  !
  !================================================================================================================================
  !

  !>Solve an optimiser solver
  SUBROUTINE SolverOptimiser_Solve(optimiserSolver,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("SolverOptimiser_Solve",err,error,*999)

    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    CALL FlagError("Not implemented.",err,error,*999)
         
    EXITS("SolverOptimiser_Solve")
    RETURN
999 ERRORSEXITS("SolverOptimiser_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_Solve

  !
  !================================================================================================================================
  !

  !>Sets/changes the type of variable type for an optimisation solver.
  SUBROUTINE SolverOptimiser_VariableTypeSet(optimiserSolver,solverVariableType,err,error,*)

    !Argument variables
    TYPE(OptimiserSolverType), POINTER :: optimiserSolver !<A pointer the optimiser solver to set the variable type for.
    INTEGER(INTG), INTENT(IN) :: solverVariableType !<The type of variable for the optimiser solver to set. \see SolverRoutines_OptimiserVariableType,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolverOptimiser_VariableTypeSet",err,error,*999)
    
    IF(.NOT.ASSOCIATED(optimiserSolver)) CALL FlagError("Optimiser solver is not associated.",err,error,*999)
    
    SELECT CASE(solverVariableType)
    CASE(SOLVER_OPTIMISER_CONTINUOUS_VARIABLES)
      optimiserSolver%variableType=SOLVER_OPTIMISER_CONTINUOUS_VARIABLES
    CASE(SOLVER_OPTIMISER_DISCRETE_VARIABLES)
      optimiserSolver%variableType=SOLVER_OPTIMISER_DISCRETE_VARIABLES
    CASE DEFAULT
      localError="The specified variable type of "//TRIM(NumberToVString(solverVariableType,"*",err,error))// &
        & " is invalid for an optimiser solver."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("SolverOptimiser_VariableTypeSet")
    RETURN
999 ERRORSEXITS("SolverOptimiser_VariableTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE SolverOptimiser_VariableTypeSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for a solver. \see OpenCMISS::Iron::cmfe_Solver_OutputTypeSet
  SUBROUTINE Solver_OutputTypeSet(solver,outputType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The type of solver output to be set \see SolverRoutines_OutputTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Solver_OutputTypeSet",err,error,*999)

    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    SELECT CASE(outputType)
    CASE(SOLVER_NO_OUTPUT)
      solver%outputType=SOLVER_NO_OUTPUT
    CASE(SOLVER_MONITOR_OUTPUT)
      solver%outputType=SOLVER_MONITOR_OUTPUT
    CASE(SOLVER_PROGRESS_OUTPUT)
      solver%outputType=SOLVER_PROGRESS_OUTPUT
    CASE(SOLVER_TIMING_OUTPUT)
      solver%outputType=SOLVER_TIMING_OUTPUT
    CASE(SOLVER_SOLVER_OUTPUT)
      solver%outputType=SOLVER_SOLVER_OUTPUT
    CASE(SOLVER_MATRIX_OUTPUT)
      solver%outputType=SOLVER_MATRIX_OUTPUT         
    CASE DEFAULT
      localError="The specified solver output type of "// &
        & TRIM(NumberToVString(outputType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("Solver_OutputTypeSet")
    RETURN
999 ERRORSEXITS("Solver_OutputTypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_OutputTypeSet
        
  !
  !================================================================================================================================
  !

  !>Updates the solver solution from the field variables
  SUBROUTINE Solver_SolutionUpdate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to update the solution from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnNumber,equationsSetIdx,interfaceConditionIdx,localDOFIdx,localNumber,numberOfDOFs, &
      & numberOfEquationsSets,numberOfInterfaceConditions,numberOfSolverMatrices,numberOfVariables,solverDOFIdx, &
      & solverMatrixIdx,variableDOFIdx,variableIdx,variableType
    REAL(DP) :: additiveConstant,dofValue,couplingCoefficient
    REAL(DP), POINTER :: variableData(:)
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(DomainMappingType), POINTER :: columnDOFsMapping,domainMapping
    TYPE(EquationsMatricesToSolverMatrixMapType), POINTER :: equationsMatricesToSolverMatrixMap
    TYPE(EquationsSetToSolverMatricesMapType), POINTER :: equationsSetToSolverMatricesMap
    TYPE(FieldType), POINTER :: dependentField,lagrangeField
    TYPE(FieldVariableType), POINTER :: dependentVariable,lagrangeVariable
    TYPE(InterfaceConditionToSolverMatricesMapType), POINTER :: interfaceConditionToSolverMatricesMap
    TYPE(InterfaceMatricesToSolverMatrixMapType), POINTER :: interfaceMatricesToSolverMatrixMap
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VariableDOFToSolverDOFsMapType), POINTER :: varDOFToSolverDOFsMap
 
    NULLIFY(variableData)
    
    ENTERS("Solver_SolutionUpdate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
      NULLIFY(solverMatrixToEquationsMap)
      CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
      NULLIFY(columnDOFsMapping)
      CALL SolverMappingSMToEQSMap_ColumnDOFsMappingGet(solverMatrixToEquationsMap,columnDOFsMapping,err,error,*999)
      NULLIFY(solverVector)
      CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSetToSolverMatricesMap)
        CALL SolverMapping_EquationsSetToSolverMatricesMapGet(solverMapping,equationsSetIdx,equationsSetToSolverMatricesMap, &
          & err,error,*999)
        NULLIFY(equationsMatricesToSolverMatrixMap)
        CALL SolverMappingESToSMSMap_EquationsMatricesToSolverMatrixMapGet(equationsSetToSolverMatricesMap,solverMatrixIdx, &
          & equationsMatricesToSolverMatrixMap,err,error,*999)
        CALL SolverMappingEMSToSMMap_NumberOfVariablesGet(equationsMatricesToSolverMatrixMap,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL SolverMappingEMSToSMMap_VariableGet(equationsMatricesToSolverMatrixMap,variableIdx,dependentVariable, &
            & err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(dependentVariable,domainMapping,err,error,*999)
          NULLIFY(varDOFToSolverDOFsMap)
          CALL SolverMappingEMSToSMMap_VariableDOFToSolverDOFsMapGet(equationsMatricesToSolverMatrixMap,variableIdx, &
            & varDOFToSolverDOFsMap,err,error,*999)
          NULLIFY(variableData)
          CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_VALUES_SET_TYPE,variableData,err,error,*999)
          CALL FieldVariable_NumberOfDOFsGet(dependentVariable,numberOfDOFs,err,error,*999)
          DO variableDOFIdx=1,numberOfDOFs
            CALL SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet(varDOFToSolverDOFsMap,variableDOFIdx,solverDOFIdx, &
              & couplingCoefficient,additiveConstant,err,error,*999)
            IF(solverDOFIdx/=0) THEN
              dofValue=variableData(variableDOFIdx)*couplingCoefficient+additiveConstant
              CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,solverDOFIdx,1,localDOFIdx,err,error,*999)
              CALL DistributedVector_ValuesSet(solverVector,localDOFIdx,dofValue,err,error,*999)
            ENDIF
          ENDDO !variableDOFIdx
          CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_VALUES_SET_TYPE,variableData,err,error,*999)
        ENDDO !variableIdx
      ENDDO !equationsSetIdx
      !Loop over interface conditions
      CALL SolverMapping_NumberOfInterfaceConditionsGet(solverMapping,numberOfInterfaceConditions,err,error,*999)
      DO interfaceConditionIdx=1,numberOfInterfaceConditions
        NULLIFY(interfaceConditionToSolverMatricesMap)
        CALL SolverMapping_InterfaceConditionToSolverMatricesMapGet(solverMapping,interfaceConditionIdx, &
          & interfaceConditionToSolverMatricesMap,err,error,*999)
        NULLIFY(interfaceMatricesToSolverMatrixMap)
        CALL SolverMappingICToSMSMap_InterfaceMatricesToSolverMatrixMapGet(interfaceConditionToSolverMatricesMap,solverMatrixIdx, &
          & interfaceMatricesToSolverMatrixMap,err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL SolverMappingIMSToSMMap_LagrangeVariableGet(interfaceMatricesToSolverMatrixMap,lagrangeVariable,err,error,*999)        
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(lagrangeVariable,domainMapping,err,error,*999)
        NULLIFY(varDOFToSolverDOFsMap)
        CALL SolverMappingIMSToSMMap_LagrangeVarDOFToSolverDOFsMapGet(interfaceMatricesToSolverMatrixMap,varDOFToSolverDOFsMap, &
          & err,error,*999)
        NULLIFY(variableData)
        CALL FieldVariable_ParameterSetDataGet(lagrangeVariable,FIELD_VALUES_SET_TYPE,variableData,err,error,*999)
        CALL FieldVariable_NumberOfDOFsGet(lagrangeVariable,numberOfDOFs,err,error,*999)
        DO variableDOFIdx=1,numberOfDOFs
          CALL SolverMappingVDOFToSDOFsMap_SolverDOFCouplingGet(varDOFToSolverDOFsMap,variableDOFIdx,solverDOFIdx, &
            & couplingCoefficient,additiveConstant,err,error,*999)
          IF(solverDOFIdx/=0) THEN
            dofValue=variableData(variableDOFIdx)*couplingCoefficient+additiveConstant
            CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,solverDOFIdx,1,localDOFIdx,err,error,*999)
            CALL DistributedVector_ValuesSet(solverVector,localDOFIdx,dofValue,err,error,*999)
          ENDIF
        ENDDO !variableDOFIdx
        CALL FieldVariable_ParameterSetDataRestore(lagrangeVariable,FIELD_VALUES_SET_TYPE,variableData,err,error,*999)
      ENDDO !interfaceConditionIdx
      CALL DistributedVector_UpdateStart(solverVector,err,error,*999)
      CALL DistributedVector_UpdateFinish(solverVector,err,error,*999)
    ENDDO !solverMatrixIdx
    
    EXITS("Solver_SolutionUpdate")
    RETURN
999 ERRORSEXITS("Solver_SolutionUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_SolutionUpdate
  
  !
  !================================================================================================================================
  !

  !>Solve the problem. 
  RECURSIVE SUBROUTINE Solver_Solve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    REAL(SP) :: systemElapsed,systemTime1(1),systemTime2(1),userElapsed,userTime1(1),userTime2(1)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_Solve",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)

    IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)          
    ENDIF
    
    !Solve the system depending on the solver type
    SELECT CASE(solver%solveType)
    CASE(SOLVER_LINEAR_TYPE)
      !Solve linear equations
      CALL SolverLinear_Solve(solver%linearSolver,err,error,*999)
    CASE(SOLVER_NONLINEAR_TYPE)
      !Solve nonlinear equations
      CALL SolverNonlinear_Solve(solver%nonlinearSolver,err,error,*999)
    CASE(SOLVER_DYNAMIC_TYPE)
      !Solve dynamic equations
      CALL SolverDynamic_Solve(solver%dynamicSolver,err,error,*999)
    CASE(SOLVER_DAE_TYPE)
      !Solve differential-algebraic equations
      CALL SolverDAE_Solve(solver%DAESolver,err,error,*999)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      !Solve eigenproblem
      CALL SolverEigenproblem_Solve(solver%eigenproblemSolver,err,error,*999)
    CASE(SOLVER_OPTIMISER_TYPE)
      !Solve an optimisation problem
      CALL SolverOptimiser_Solve(solver%optimiserSolver,err,error,*999)
    CASE(SOLVER_CELLML_EVALUATOR_TYPE)
      !Solve a CellML evaluator
      CALL SolverCellMLEvaluator_Solve(solver%cellMLEvaluatorSolver,err,error,*999)
    CASE DEFAULT
      localError="The solver type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !If necessary output the timing information
    IF(solver%outputType>=SOLVER_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      IF(solver%outputType>=SOLVER_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Total time for solve",userElapsed,systemElapsed,err,error,*999)
    ENDIF
       
    EXITS("Solver_Solve")
    RETURN
999 ERRORSEXITS("Solver_Solve",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_Solve

  !
  !================================================================================================================================
  !

  !>Sets/changes the type for a solver.
  SUBROUTINE Solver_TypeSet(solver,solveType,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to set the solver type for.
    INTEGER(INTG), INTENT(IN) :: solveType !<The type of solver to be set \see SolverRoutines_SolverTypes,SolverRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("Solver_TypeSet",err,error,*998)

    CALL Solver_AssertNotFinished(solver,err,error,*998)
    IF(ASSOCIATED(solver%linkingSolver)) &
      & CALL FlagError("Can not changed the solver type for a solve that has been linked.",err,error,*998)
        
    IF(solveType/=solver%solveType) THEN
      !Initialise the new solver type 
      SELECT CASE(solveType)
      CASE(SOLVER_LINEAR_TYPE)
        CALL Solver_LinearInitialise(solver,err,error,*999)
      CASE(SOLVER_NONLINEAR_TYPE)
        CALL Solver_NonlinearInitialise(solver,err,error,*999)
      CASE(SOLVER_DYNAMIC_TYPE)
        CALL Solver_DynamicInitialise(solver,err,error,*999)
      CASE(SOLVER_DAE_TYPE)
        CALL Solver_DAEInitialise(solver,err,error,*999)
      CASE(SOLVER_EIGENPROBLEM_TYPE)
        CALL Solver_EigenproblemInitialise(solver,err,error,*999)
      CASE(SOLVER_OPTIMISER_TYPE)
        CALL Solver_OptimiserInitialise(solver,err,error,*999)
      CASE(SOLVER_CELLML_EVALUATOR_TYPE)
        CALL Solver_CellMLEvaluatorInitialise(solver,err,error,*999)
      CASE(SOLVER_GEOMETRIC_TRANSFORMATION_TYPE)
        CALL Solver_GeometricTransformationInitialise(solver,err,error,*999)
      CASE DEFAULT
        localError="The specified solve type of "//TRIM(NumberToVString(solveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Finalise the old solve type
      SELECT CASE(solver%solveType)
      CASE(SOLVER_LINEAR_TYPE)
        CALL Solver_LinearFinalise(solver%linearSolver,err,error,*999)
      CASE(SOLVER_NONLINEAR_TYPE)
        CALL Solver_NonlinearFinalise(solver%nonlinearSolver,err,error,*999)
      CASE(SOLVER_DYNAMIC_TYPE)
        CALL Solver_DynamicFinalise(solver%dynamicSolver,err,error,*999)
      CASE(SOLVER_DAE_TYPE)
        CALL Solver_DAEFinalise(solver%DAESolver,err,error,*999)
      CASE(SOLVER_EIGENPROBLEM_TYPE)
        CALL Solver_EigenproblemFinalise(solver%eigenproblemSolver,err,error,*999)
      CASE(SOLVER_OPTIMISER_TYPE)
        CALL Solver_OptimiserFinalise(solver%optimiserSolver,err,error,*999)
      CASE(SOLVER_CELLML_EVALUATOR_TYPE)
        CALL Solver_CellMLEvaluatorFinalise(solver%cellMLEvaluatorSolver,err,error,*999)
      CASE(SOLVER_GEOMETRIC_TRANSFORMATION_TYPE)
        CALL Solver_GeometricTransformationFinalise(solver%geometricTransformationSolver,err,error,*999)
      CASE DEFAULT
        localError="The solver solve type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set the solve type
      solver%solveType=solveType
    ENDIF
    
    EXITS("Solver_TypeSet")
    RETURN
999 SELECT CASE(solveType)
    CASE(SOLVER_LINEAR_TYPE)
      CALL Solver_LinearFinalise(solver%linearSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_NONLINEAR_TYPE)
      CALL Solver_NonlinearFinalise(solver%nonlinearSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_DYNAMIC_TYPE)
      CALL Solver_DynamicFinalise(solver%dynamicSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_DAE_TYPE)
      CALL Solver_DAEFinalise(solver%DAESolver,dummyErr,dummyError,*998)
    CASE(SOLVER_EIGENPROBLEM_TYPE)
      CALL Solver_EigenproblemFinalise(solver%eigenproblemSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_OPTIMISER_TYPE)
      CALL Solver_OptimiserFinalise(solver%optimiserSolver,dummyErr,dummyError,*998)
    CASE(SOLVER_GEOMETRIC_TRANSFORMATION_TYPE)
      CALL Solver_GeometricTransformationFinalise(solver%geometricTransformationSolver,err,error,*999)
    END SELECT
998 ERRORSEXITS("Solver_TypeSet",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_TypeSet
        
  !
  !================================================================================================================================
  !

  !>Updates the dependent variables from the solver solution for dynamic solvers
  SUBROUTINE Solver_VariablesDynamicFieldUpdate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,dynamicVariableType,equationsDOFIdx,equationIdx,equationsSetIdx,equationType,numberOfEquationDOFs, &
      & numberOfSolverDOFs,numberOfSolverMatrices,numberOfVariables,solverDOFIdx,solverIdx,solverMatrixIdx,variableDOF,variableIdx
    REAL(DP) :: currentAcceleration,additiveConstant,deltaT,currentDisplacement,previousAcceleration, &
      & previousDisplacement,previousVelocity,alphaValue,variableCoefficient,currentVelocity
    REAL(DP), POINTER :: solverData(:)
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(solverData)
    
    ENTERS("Solver_VariablesDynamicFieldUpdate",err,error,*998)

    CALL Solver_AssertIsFinished(solver,err,error,*998)
    CALL Solver_AssertIsDynamic(solver,err,error,*998)    
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    deltaT=dynamicSolver%timeIncrement
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
      NULLIFY(solverVector)
      CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
      !Get the solver variables data
      CALL DistributedVector_DataGet(solverVector,solverData,err,error,*999)
      !Loop over the solver variable dofs
      NULLIFY(solverMatrixToEquationsMap)
      CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverIdx,solverMatrixToEquationsMap,err,error,*999)
      CALL SolverMappingSMToEQSMap_NumberOfDOFsGet(solverMatrixToEquationsMap,numberOfSolverDOFs,err,error,*999)
      DO solverDOFIdx=1,numberOfSolverDOFs
        NULLIFY(solverDOFToVariableDOFsMap)
        CALL SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,solverDOFIdx, &
          & solverDOFToVariableDOFsMap,err,error,*999)
        !Loop over the equations sets associated with this dof
        CALL SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs,err,error,*999)
        DO equationsDOFIdx=1,numberOfEquationDOFs
          CALL SolverMappingSDOFToVDOFsMap_EquationsInfoGet(solverDOFToVariableDOFsMap,equationsDOFIdx,equationType,equationIdx, &
            & err,error,*999)
          SELECT CASE(equationType)
          CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET, &
            & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION)
            NULLIFY(dependentVariable)
            CALL SolverMappingSDOFToVDOFsMap_VariableGet(solverDOFToVariableDOFsMap,equationsDOFIdx,dependentVariable, &
              & err,error,*999)
            !Get the dependent field DOF the solver DOF is mapped to
            CALL SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet(solverDOFToVariableDOFsMap,equationsDOFIdx,variableDOF, &
              & variableCoefficient,additiveConstant,err,error,*999)
            alphaValue=solverData(solverDOFIdx)*variableCoefficient+additiveConstant
            !Set the dependent field DOF
            IF(dynamicSolver%solverInitialised) THEN
              SELECT CASE(dynamicSolver%degree)
              CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,variableDOF, &
                  & previousDisplacement,err,error,*999)
                currentDisplacement=previousDisplacement+deltaT*alphaValue
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,variableDOF, &
                  & currentDisplacement,err,error,*999)
              CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,variableDOF, &
                  & previousDisplacement,err,error,*999)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,variableDOF, &
                  & previousVelocity,err,error,*999)
                currentDisplacement=previousDisplacement+deltaT*previousVelocity+(deltaT*deltaT/2.0_DP)*alphaValue
                currentVelocity=previousVelocity+deltaT*alphaValue
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,variableDOF, &
                  & currentDisplacement,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,variableDOF, &
                  & currentVelocity,err,error,*999)
              CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,variableDOF, &
                  & previousDisplacement,err,error,*999)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VELOCITY_SET_TYPE,variableDOF, &
                  & previousVelocity,err,error,*999)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_ACCELERATION_SET_TYPE,variableDOF, &
                  & previousAcceleration,err,error,*999)
                currentDisplacement=previousDisplacement+deltaT*previousVelocity+ &
                  & (deltaT*deltaT/2.0_DP)*previousAcceleration+ &
                  & (deltaT*deltaT*deltaT/6.0_DP)*alphaValue
                currentVelocity=previousVelocity+deltaT*previousAcceleration+(deltaT*deltaT/2.0_DP)*alphaValue
                currentAcceleration=previousAcceleration+deltaT*alphaValue
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,variableDOF, &
                  & currentDisplacement,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,variableDOF, &
                  & currentVelocity,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,variableDOF, &
                  & currentAcceleration,err,error,*999)
              CASE DEFAULT
                localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              SELECT CASE(dynamicSolver%order)
              CASE(SOLVER_DYNAMIC_FIRST_ORDER)
                SELECT CASE(dynamicSolver%degree)
                CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                  !Do nothing
                CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,variableDOF, &
                    & alphaValue,err,error,*999)
                CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,variableDOF, &
                    & alphaValue,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE, &
                    & variableDOF,0.0_DP,err,error,*999)
                CASE DEFAULT
                  localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(SOLVER_DYNAMIC_SECOND_ORDER)
                IF(dynamicSolver%degree==SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE, &
                    & variableDOF,alphaValue,err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="The dynamic solver order of "//TRIM(NumberToVString(dynamicSolver%order,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE DEFAULT
            localError="The equations type of "//TRIM(NumberToVString(equationType,"*",err,error))//" of equations DOF index "// &
              & TRIM(NumberToVString(equationsDOFIdx,"*",err,error))//" for solver DOF index "// &
              & TRIM(NumberToVString(solverDOFIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !equationsDOFIdx
      ENDDO !solverDOFIdx
      !Restore the solver DOF data
      CALL DistributedVector_DataRestore(solverVector,solverData,err,error,*999)
      !Start the transfer of the field DOFs
      NULLIFY(solverMappingVariables)
      CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
      CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(dependentVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,dependentVariable,err,error,*999)
        IF(dynamicSolver%solverInitialised) THEN
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
            CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
            IF(dynamicSolver%degree>SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          SELECT CASE(dynamicSolver%order)
          CASE(SOLVER_DYNAMIC_FIRST_ORDER)
            SELECT CASE(dynamicSolver%degree)
            CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
              !Do nothing
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE)                                  
              CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,err,error,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
              CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE,err,error,*999)
            CASE DEFAULT
              localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(SOLVER_DYNAMIC_SECOND_ORDER)
            IF(dynamicSolver%degree==SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The dynamic solver order of "//TRIM(NumberToVString(dynamicSolver%order,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ENDDO !variableIdx
      !Finish the transfer of the field DOFs
      DO variableIdx=1,numberOfVariables
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(dependentVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,dependentVariable,err,error,*999)
        IF(dynamicSolver%solverInitialised) THEN
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(dynamicSolver%degree>SOLVER_DYNAMIC_FIRST_DEGREE) THEN
            CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,err,error,*999)
            IF(dynamicSolver%degree>SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,err,error,*999)
            ENDIF
          ENDIF
        ELSE
          SELECT CASE(dynamicSolver%order)
          CASE(SOLVER_DYNAMIC_FIRST_ORDER)
            SELECT CASE(dynamicSolver%degree)
            CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
              !Do nothing
            CASE(SOLVER_DYNAMIC_SECOND_DEGREE)                                  
              CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,err,error,*999)
            CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
              CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_INITIAL_VELOCITY_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE,err,error,*999)
            CASE DEFAULT
              localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(SOLVER_DYNAMIC_SECOND_ORDER)
            IF(dynamicSolver%degree==SOLVER_DYNAMIC_THIRD_DEGREE) THEN
              CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_INITIAL_ACCELERATION_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The dynamic solver order of "//TRIM(NumberToVString(dynamicSolver%order,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      ENDDO !variableIdx
    ENDDO !solverMatrixIdx
    
    EXITS("Solver_VariablesDynamicFieldUpdate")
    RETURN
999 IF(ASSOCIATED(solverData)) CALL DistributedVector_DataRestore(solverVector,solverData,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_VariablesDynamicFieldUpdate",err,error)
    RETURN 1
   
  END SUBROUTINE Solver_VariablesDynamicFieldUpdate

  !
  !================================================================================================================================
  !

  !>Updates the previous values from the solver solution for dynamic solvers
  SUBROUTINE Solver_VariablesDynamicFieldPreviousValuesUpdate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,linearityType,numberOfEquationsSets,numberOfResiduals,numberOfSolverMatrices, &
      & numberOfVariables,residualIdx,residualVariableIdx,solverMatrixIdx,variableIdx
    TYPE(DistributedVectorType), POINTER :: currentVector,previousVector,previous2Vector,previous3Vector
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: field
    TYPE(FieldVariableType), POINTER :: fieldVariable,residualVariable
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solver_VariablesDynamicFieldPreviousValuesUpdate",err,error,*999)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    CALL Solver_AssertIsDynamic(solver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    !Loop over the solver matrices
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverMatrixToEquationsMap)
      CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
      !Loop over the variables involved in the solver matrix.
      NULLIFY(solverMappingVariables)
      CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
      CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(fieldVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
        !Copy the displacements
        CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP, &
          & err,error,*999)
        IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
          !Copy velocity
          CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_VELOCITY_VALUES_SET_TYPE,FIELD_PREVIOUS_VELOCITY_SET_TYPE, &
            & 1.0_DP,err,error,*999)
          IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
            !Copy acceleration
            CALL FieldVariable_ParameterSetsCopy(fieldVariable,FIELD_ACCELERATION_VALUES_SET_TYPE, &
              & FIELD_PREVIOUS_ACCELERATION_SET_TYPE,1.0_DP,err,error,*999)
          ENDIF
        ENDIF
      ENDDO !variableIdx
      IF(dynamicSolver%linearity==SOLVER_DYNAMIC_NONLINEAR) THEN
        !Loop over the equations sets and copy any residuals
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeGet(equations,linearityType,err,error,*999)
          IF(linearityType==EQUATIONS_NONLINEAR) THEN
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            NULLIFY(vectorMapping)
            CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
            NULLIFY(vectorMatrices)
            CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
            NULLIFY(nonlinearMapping)
            CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
            NULLIFY(nonlinearMatrices)
            CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
            CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
            DO residualIdx=1,numberOfResiduals
              NULLIFY(residualVector)
              CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
              NULLIFY(currentVector)
              CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                & currentVector,err,error,*999)
              NULLIFY(previousVector)
              CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS_VECTOR, &
                & previousVector,err,error,*999)
              IF(dynamicSolver%degree>=SOLVER_DYNAMIC_SECOND_DEGREE) THEN
                NULLIFY(previous2Vector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                  & previous2Vector,err,error,*999)
                IF(dynamicSolver%degree>=SOLVER_DYNAMIC_THIRD_DEGREE) THEN
                  NULLIFY(previous3Vector)
                  CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_PREVIOUS2_VECTOR, &
                    & previous3Vector,err,error,*999)
                  CALL DistributedVector_VectorCopy(previous3Vector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,1.0_DP, &
                    & previous2Vector,err,error,*999)
                ENDIF
                CALL DistributedVector_VectorCopy(previous2Vector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,1.0_DP, &
                  & previousVector,err,error,*999)                
              ENDIF
              CALL DistributedVector_VectorCopy(previousVector,DISTRIBUTED_MATRIX_VECTOR_INCLUDE_GHOSTS_TYPE,1.0_DP, &
                & currentVector,err,error,*999)
            ENDDO !residualIdx
          ENDIF
        ENDDO !equationsSetIdx
      ENDIF
    ENDDO !solverMatrixIdx

    EXITS("Solver_VariablesDynamicFieldPreviousValuesUpdate")
    RETURN
999 ERRORS("Solver_VariablesDynamicFieldPreviousValuesUpdate",err,error)
    EXITS("Solver_VariablesDynamicFieldPreviousValuesUpdate")
    RETURN 1

  END SUBROUTINE Solver_VariablesDynamicFieldPreviousValuesUpdate

  !
  !================================================================================================================================
  !

  !>Update the field values form the dynamic factor * current solver values AND add in previous values
  SUBROUTINE Solver_VariablesDynamicNonlinearUpdate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    INTEGER(INTG) :: dummyErr,dynamicVariableType,equationsDOFIdx,equationIdx,equationsSetIdx,equationType,numberOfDOFs, &
      & numberOfEquationDOFs,numberOfSolverMatrices,numberOfVariables,solverDOFIdx,solverMatrixIdx,variableDOF
    REAL(DP) :: additiveConstant,alphaValue,alphaDOFValue,currentDisplacement,deltaT,previousDisplacement, &
      & previousVelocity,previousAcceleration,variableCoefficient
    INTEGER(INTG) :: variableIdx,variableType,interfaceConditionIdx
    REAL(DP), POINTER :: solverData(:)
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(SolverType), POINTER :: linkingSolver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(solverData)
    
    ENTERS("Solver_VariablesDynamicNonlinearUpdate",err,error,*998)

    CALL Solver_AssertIsFinished(solver,err,error,*998)
    NULLIFY(linkingSolver)
    CALL Solver_LinkingSolverGet(solver,linkingSolver,err,error,*999)
    CALL Solver_AssertIsDynamic(linkingSolver,err,error,*999)
    NULLIFY(dynamicSolver)
    CALL Solver_DynamicSolverGet(linkingSolver,dynamicSolver,err,error,*999)
    !Define the dynamic alpha factor
    IF(dynamicSolver%solverInitialised) THEN
      deltaT=dynamicSolver%timeIncrement
      !Set the dependent field for calculating the nonlinear residual and Jacobian values
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMatrices)
      CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
      DO solverMatrixIdx=1,numberOfSolverMatrices
        NULLIFY(solverMatrix)
        CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
        NULLIFY(solverVector)
        CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
        !Get the solver variables data                  
        CALL DistributedVector_DataGet(solverVector,solverData,err,error,*999)
        !Loop over the solver variable dofs
        NULLIFY(solverMatrixToEquationsMap)
        CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
        CALL SolverMappingSMToEQSMap_NumberOfDOFsGet(solverMatrixToEquationsMap,numberOfDOFS,err,error,*999)
        DO solverDOFIdx=1,numberOfDOFs
          NULLIFY(solverDOFToVariableDOFsMap)
          CALL SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,solverDOFIdx, &
            & solverDOFToVariableDOFsMap,err,error,*999)
          !Loop over the equations associated with this dof
          CALL SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs,err,error,*999)
          DO equationsDOFIdx=1,numberOfEquationDOFs
            CALL SolverMappingSDOFToVDOFsMap_EquationsInfoGet(solverDOFToVariableDOFsMap,equationsDOFIdx,equationType, &
              & equationIdx,err,error,*999)
            SELECT CASE(equationType)
            CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET, &
              & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION)
              !Get the field variable
              NULLIFY(dependentVariable)
              CALL SolverMappingSDOFToVDOFsMap_VariableGet(solverDOFToVariableDOFsMap,equationsDOFIdx,dependentVariable, &
                & err,error,*999)
              !Get the dependent field variable DOF the solver DOF is mapped to
              CALL SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet(solverDOFToVariableDOFsMap,equationsDOFIdx,variableDOF, &
                & variableCoefficient,additiveConstant,err,error,*999)
              alphaValue=solverData(solverDOFIdx)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,variableDOF, &
                & alphaValue,err,error,*999)
              !Calculate solver data only
              alphaDOFValue=alphaValue*variableCoefficient+additiveConstant
              !Get the predicted displacement data       
              CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,variableDOF, &
                & previousDisplacement,err,error,*999)
              SELECT CASE(dynamicSolver%degree)
              CASE(SOLVER_DYNAMIC_FIRST_DEGREE)
                currentDisplacement=previousDisplacement+deltaT*alphaDOFValue
              CASE(SOLVER_DYNAMIC_SECOND_DEGREE)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,variableDOF, &
                  & previousVelocity,err,error,*999)
                currentDisplacement=previousDisplacement+deltaT*previousVelocity+deltaT*deltaT*alphaDOFValue/2.0_DP
              CASE(SOLVER_DYNAMIC_THIRD_DEGREE)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_VELOCITY_VALUES_SET_TYPE,variableDOF, &
                  & previousVelocity,err,error,*999)
                CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,variableDOF, &
                  & previousAcceleration,err,error,*999)
                currentDisplacement=previousDisplacement+deltaT*previousVelocity+ &
                  & deltaT*deltaT*previousAcceleration/2.0_DP+deltaT+deltaT*deltaT*alphaDOFValue/6.0_DP                  
              CASE DEFAULT
                localError="The dynamic solver degree of "//TRIM(NumberToVString(dynamicSolver%degree,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,variableDOF, &
                & currentDisplacement,err,error,*999)
            CASE DEFAULT
              localError="The equations type of "//TRIM(NumberToVString(equationType,"*",err,error))//" of equations index "// &
                & TRIM(NumberToVString(equationsDOFIdx,"*",err,error))//" for solver degree-of-freedom "// &
                & TRIM(NumberToVString(solverDOFIdx,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsDOFIdx
        ENDDO !solverDOFIdx
        !Restore the solver dof data
        CALL DistributedVector_DataRestore(solverVector,solverData,err,error,*999)
        !Start the transfer of the field dofs
        NULLIFY(solverMappingVariables)
        CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
        CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(solverMappingVariable)
          CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
          NULLIFY(dependentVariable)
          CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
        ENDDO !variableIdx
        !Finish the transfer of the field DOFs
        DO variableIdx=1,numberOfVariables
          NULLIFY(solverMappingVariable)
          CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
          NULLIFY(dependentVariable)
          CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,err,error,*999)
        ENDDO !variableIdx

      ENDDO !solverMatrixIdx     
    ENDIF !solver initialised
    
    EXITS("Solver_VariablesDynamicNonlinearUpdate")
    RETURN
999 IF(ASSOCIATED(solverData)) CALL DistributedVector_DataRestore(solverVector,solverData,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_VariablesDynamicNonlinearUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_VariablesDynamicNonlinearUpdate


  !
  !================================================================================================================================
  !

  !>Updates the dependent variables from the solver solution for static solvers
  SUBROUTINE Solver_VariablesFieldUpdate(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer the solver to update the variables from
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsDOFIdx,equationIdx,equationsSetIdx,equationType,numberOfEquationDOFs, &
      & numberOfSolverDOFs,numberOfSolverMatrices,numberOfVariables,solverDOFIdx,solverMatrixIdx,variableDOF, &
      & variableIdx,variableType
    REAL(DP) :: additiveConstant,solverDOFValue,variableCoefficient
    REAL(DP), POINTER :: solverData(:)
    TYPE(DistributedVectorType), POINTER :: solverVector
    TYPE(FieldVariableType), POINTER :: dependentVariable,fieldVariable,lagrangeVariable
    TYPE(SolverDOFToVariableDOFsMapType), POINTER :: solverDOFToVariableDOFsMap
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMappingVariableType), POINTER :: solverMappingVariable
    TYPE(SolverMappingVariablesType), POINTER :: solverMappingVariables
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverMatrixToEquationsMapType), POINTER :: solverMatrixToEquationsMap
    TYPE(VARYING_STRING) :: dummyError,localError

    NULLIFY(solverData)
    
    ENTERS("Solver_VariablesFieldUpdate",err,error,*998)

    CALL Solver_AssertIsFinished(solver,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
    NULLIFY(solverMatrices)
    CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
    DO solverMatrixIdx=1,numberOfSolverMatrices
      NULLIFY(solverMatrix)
      CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
      NULLIFY(solverVector)
      CALL SolverMatrix_SolverDistributedVectorGet(solverMatrix,solverVector,err,error,*999)
      !Get the solver variables data7
      CALL DistributedVector_DataGet(solverVector,solverData,err,error,*999)
      !Loop over the solver variable dofs
      NULLIFY(solverMatrixToEquationsMap)
      CALL SolverMapping_SolverMatrixToEquationsMapGet(solverMapping,solverMatrixIdx,solverMatrixToEquationsMap,err,error,*999)
      CALL SolverMappingSMToEQSMap_NumberOfDOFsGet(solverMatrixToEquationsMap,numberOfSolverDOFs,err,error,*999)            
      DO solverDOFIdx=1,numberOfSolverDOFs
        NULLIFY(solverDOFToVariableDOFsMap)
        CALL SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,solverDOFIdx, &
          & solverDOFToVariableDOFsMap,err,error,*999)
        !Loop over the equations associated with this dof
        CALL SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs,err,error,*999)
        DO equationsDOFIdx=1,numberOfEquationDOFs
          !Get the field variable DOF the solver DOF is mapped to
          CALL SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet(solverDOFToVariableDOFsMap,equationsDOFIdx,variableDOF, &
            & variableCoefficient,additiveConstant,err,error,*999)
          !Calculate the field variable DOF value
          solverDOFValue=solverData(solverDOFIdx)*variableCoefficient+additiveConstant
          CALL SolverMappingSDOFToVDOFsMap_EquationsInfoGet(solverDOFToVariableDOFsMap,equationsDOFIdx,equationType,equationIdx, &
            & err,error,*999)
          SELECT CASE(equationType)
          CASE(SOLVER_MAPPING_EQUATIONS_EQUATIONS_SET, &
            & SOLVER_MAPPING_EQUATIONS_INTERFACE_CONDITION)
            !Equations set DOF.
            NULLIFY(fieldVariable)
            CALL SolverMappingSDOFToVDOFsMap_VariableGet(solverDOFToVariableDOFsMap,equationsDOFIdx,fieldVariable, &
              & err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(fieldVariable,FIELD_VALUES_SET_TYPE,variableDOF,solverDOFValue, &
              & err,error,*999)
          CASE DEFAULT
            localError="The equations type of "//TRIM(NumberToVString(equationType,"*",err,error))// &
              & " of equations DOF index "//TRIM(NumberToVString(equationsDOFIdx,"*",err,error))//" for solver DOF index "// &
              & TRIM(NumberToVString(solverDOFIdx,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !equationsDOFIdx
      ENDDO !solverDOFIdx
      IF(diagnostics2) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  Solver matrix index = ",solverMatrixIdx,err,error,*999)
        DO solverDOFIdx=1,numberOfSolverDOFs
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Solver DOF index = ",solverDOFIdx,err,error,*999)
          NULLIFY(solverDOFToVariableDOFsMap)
          CALL SolverMappingSMToEQSMap_SolverDOFToVariableDOFsMapGet(solverMatrixToEquationsMap,solverDOFIdx, &
            & solverDOFToVariableDOFsMap,err,error,*999)
          !Loop over the equations DOFs associated with this DOF
          CALL SolverMappingSDOFToVDOFsMap_NumberOfEquationDOFsGet(solverDOFToVariableDOFsMap,numberOfEquationDOFs,err,error,*999)
          DO equationsDOFIdx=1,numberOfEquationDOFs
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"      Equations DOF index = ",equationsDOFIdx,err,error,*999)
            CALL SolverMappingSDOFToVDOFsMap_VariableDOFCouplingGet(solverDOFToVariableDOFsMap,equationsDOFIdx,variableDOF, &
              & variableCoefficient,additiveConstant,err,error,*999)
            solverDOFValue=solverData(solverDOFIdx)*variableCoefficient+additiveConstant
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable dof = ",variableDOF,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Variable coefficient = ",variableCoefficient,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Additive constant = ",additiveConstant,err,error,*999)
            CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"        Value = ",solverDOFValue,err,error,*999)
          ENDDO
        ENDDO
      ENDIF
      !Restore the solver DOF data
      CALL DistributedVector_DataRestore(solverVector,solverData,err,error,*999)
      !Start the transfer of the field dofs
      NULLIFY(solverMappingVariables)
      CALL SolverMappingSMToEQSMap_VariablesListGet(solverMatrixToEquationsMap,solverMappingVariables,err,error,*999)
      CALL SolverMappingVariables_NumberOfVariablesGet(solverMappingVariables,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(fieldVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
      !Finish the transfer of the field dofs
      DO variableIdx=1,numberOfVariables
        NULLIFY(solverMappingVariable)
        CALL SolverMappingVariables_VariableGet(solverMappingVariables,variableIdx,solverMappingVariable,err,error,*999)
        NULLIFY(fieldVariable)
        CALL SolverMappingVariable_FieldVariableGet(solverMappingVariable,fieldVariable,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fieldVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
    ENDDO !solverMatrixIdx
  
    EXITS("Solver_VariablesFieldUpdate")
    RETURN
999 IF(ASSOCIATED(solverData)) CALL DistributedVector_DataRestore(solverVector,solverData,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solver_VariablesFieldUpdate",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_VariablesFieldUpdate
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of solvers.
  SUBROUTINE Solvers_CreateFinish(solvers,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverIdx
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SolverType), POINTER :: SOLVER
   
    ENTERS("Solvers_CreateFinish",err,error,*999)

    CALL Solvers_AssertNotFinished(solvers,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    !Finish the solver creation
    IF(.NOT.ALLOCATED(solvers%solvers)) CALL FlagError("Solvers solvers is not allocated.",err,error,*999)
    DO solverIdx=1,solvers%numberOfSolvers
      NULLIFY(solver)
      CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
      CALL Solver_CreateFinish(solver,err,error,*999)
    ENDDO !solverIdx            
    solvers%solversFinished=.TRUE.
     
    EXITS("Solvers_CreateFinish")
    RETURN
999 ERRORSEXITS("Solvers_CreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_CreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of a solvers for the control loop. 
  SUBROUTINE Solvers_CreateStart(controlLoop,solvers,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to create the solvers for
    TYPE(SolversType), POINTER :: solvers !<On exit, a pointer to the solvers. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfSubLoops
    TYPE(VARYING_STRING) :: localError

    ENTERS("Solvers_CreateStart",err,error,*999)

    IF(ASSOCIATED(solvers)) CALL FlagError("Solvers is already associated.",err,error,*999)
    CALL ControlLoop_AssertIsFinished(controlLoop,err,error,*999)
    CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
    
    IF(numberOfSubLoops/=0) THEN
      localError="Invalid control loop setup. The specified control loop has "// &
        & TRIM(NumberToVString(controlLoop%numberOfSubLoops,"*",err,error))// &
        & " sub loops. To create solvers the control loop must have 0 sub loops."          
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    NULLIFY(solvers)
    !Initialise the solvers
    CALL Solvers_Initialise(controlLoop,err,error,*999)
    !Return the pointer
    solvers=>controlLoop%solvers
    
    EXITS("Solvers_CreateStart")
    RETURN
999 ERRORSEXITS("Solvers_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Destroys the solvers
  SUBROUTINE Solvers_Destroy(solvers,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Solvers_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(solvers)) CALL FlagError("Solvers is not associated.",err,error,*999)
    
    CALL Solvers_Finalise(solvers,err,error,*999)
      
    EXITS("Solvers_Destroy")
    RETURN
999 ERRORSEXITS("Solvers_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_Destroy

  !
  !================================================================================================================================
  !

  !>Finalises the solvers and deallocates all memory
  SUBROUTINE Solvers_Finalise(solvers,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverIdx
 
    ENTERS("Solvers_Finalise",err,error,*999)

    IF(ASSOCIATED(solvers)) THEN
      IF(ALLOCATED(solvers%solvers)) THEN
        DO solverIdx=1,SIZE(solvers%solvers,1)
          CALL Solver_Finalise(solvers%solvers(solverIdx)%ptr,err,error,*999)
        ENDDO !solverIdx
        DEALLOCATE(solvers%solvers)
      ENDIF
      DEALLOCATE(solvers)
    ENDIF
       
    EXITS("Solvers_Finalise")
    RETURN
999 ERRORSEXITS("Solvers_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_Finalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the solvers for a control loop.
  SUBROUTINE Solvers_Initialise(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to initialise the solvers for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,solverIdx
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("Solvers_Initialise",err,error,*998)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*998)
    IF(ASSOCIATED(controlLoop%solvers)) CALL FlagError("Solvers is already allocated for this control loop.",err,error,*998)
    
    ALLOCATE(controlLoop%solvers,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate control loop solvers.",err,error,*999)
    controlLoop%solvers%controlLoop=>controlLoop
    controlLoop%solvers%solversFinished=.FALSE.
    controlLoop%solvers%numberOfSolvers=1
    ALLOCATE(controlLoop%solvers%solvers(controlLoop%solvers%numberOfSolvers),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate solvers solvers.",err,error,*999)
    DO solverIdx=1,controlLoop%solvers%numberOfSolvers
      NULLIFY(controlLoop%solvers%solvers(solverIdx)%ptr)
      CALL Solvers_SolverInitialise(controlLoop%solvers,solverIdx,err,error,*999)
    ENDDO !solverIdx
       
    EXITS("Solvers_Initialise")
    RETURN
999 CALL Solvers_Finalise(controlLoop%solvers,dummyErr,dummyError,*998)
998 ERRORSEXITS("Solvers_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_Initialise
  
 
  !
  !================================================================================================================================
  !

  !>Sets/changes the number of solvers.
  SUBROUTINE Solvers_NumberOfSolversSet(solvers,numberOfSolvers,err,error,*)

    !Argument variables
    TYPE(SolversType), POINTER :: solvers !<A pointer to the solvers to set the number for
    INTEGER(INTG), INTENT(IN) :: numberOfSolvers !<The number of solvers to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: solverIdx, oldNumberOfSolvers
    TYPE(SolverPtrType), ALLOCATABLE :: oldSolvers(:)
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Solvers_NumberOfSolversSet",err,error,*998)

    CALL Solvers_AssertNotFinished(solvers,err,error,*999)
    IF(numberOfSolvers<=0) THEN
      localError="The specified number of solvers of "//TRIM(NumberToVString(numberOfSolvers,"*",err,error))// &
        & " is invalid. The number of solvers must be > 0."
      CALL FlagError(localError,err,error,*998)
    ENDIF
    
    oldNumberOfSolvers=solvers%numberOfSolvers
    IF(numberOfSolvers/=oldNumberOfSolvers) THEN
      ALLOCATE(oldSolvers(oldNumberOfSolvers),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate old solvers.",err,error,*999)
      DO solverIdx=1,oldNumberOfSolvers
        oldSolvers(solverIdx)%ptr=>solvers%solvers(solverIdx)%ptr
      ENDDO !solverIdx
      IF(ALLOCATED(solvers%solvers)) DEALLOCATE(solvers%solvers)
      ALLOCATE(solvers%solvers(numberOfSolvers),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate solvers.",err,error,*999)
      IF(numberOfSolvers>oldNumberOfSolvers) THEN
        DO solverIdx=1,oldNumberOfSolvers
          solvers%solvers(solverIdx)%ptr=>oldSolvers(solverIdx)%ptr
        ENDDO !solverIdx
        solvers%numberOfSolvers=numberOfSolvers
        DO solverIdx=oldNumberOfSolvers+1,numberOfSolvers
          NULLIFY(solvers%solvers(solverIdx)%ptr)
          CALL Solvers_SolverInitialise(solvers,solverIdx,err,error,*999)
        ENDDO !solution_idx
      ELSE
        DO solverIdx=1,numberOfSolvers
          solvers%solvers(solverIdx)%ptr=>oldSolvers(solverIdx)%ptr
        ENDDO !solverIdx
        DO solverIdx=numberOfSolvers+1,oldNumberOfSolvers
          CALL Solver_Finalise(oldSolvers(solverIdx)%ptr,err,error,*999)
        ENDDO !solverIdx
        solvers%numberOfSolvers=numberOfSolvers
      ENDIF
    ENDIF
    
    EXITS("Solvers_NumberOfSolversSet")
    RETURN
999 IF(ALLOCATED(oldSolvers)) DEALLOCATE(oldSolvers)
998 ERRORSEXITS("Solvers_NumberOfSolversSet",err,error)
    RETURN 1
    
  END SUBROUTINE Solvers_NumberOfSolversSet
  
  !
  !================================================================================================================================
  !

  !>Adds a linked solver to the solver. Also sets the solver type for the linked solver, als well as its linking solver.
  SUBROUTINE Solver_LinkedSolverAdd(solver,solverToLink,solverTypeToLink,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to add the linked solver to.
    TYPE(SolverType), POINTER :: solverToLink !<A pointer the the solver to be linked. 
    INTEGER(INTG), INTENT(IN) :: solverTypeToLink !<The solver type of the solver to be linked.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    TYPE(SolverPtrType), ALLOCATABLE, TARGET :: oldLinkedSolvers(:)
    INTEGER(INTG) :: solverIdx

    ENTERS("Solver_LinkedSolverAdd",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solverToLink)) CALL FlagError("The solver to link is not associated.",err,error,*999)
    IF(solverTypeToLink<1.OR.solverTypeToLink<=SOLVER_NUMBER_OF_SOLVER_TYPES) THEN
      localError="The specified solver type to link of "//TRIM(NumberToVString(solverTypeToLink,"*",err,error))//&
        & " is invalid. The solver type should be >= 1 and <= "// &
        & TRIM(NumberToVString(SOLVER_NUMBER_OF_SOLVER_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Does the solver have already linked solvers?
    IF(solver%numberOfLinkedSolvers==0) THEN
      !no - then start the creation of linked solvers
      ALLOCATE(solver%linkedSolvers(1),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate linked solvers.",err,error,*999)
      DO solverIdx=1,SOLVER_NUMBER_OF_SOLVER_TYPES
        NULLIFY(solver%linkedSolverTypeMap(solverIdx)%ptr)
      ENDDO !solverIdx
      solver%linkedSolverTypeMap(solverTypeToLink)%ptr=>solverToLink
      solver%linkedSolvers(1)%ptr=>solverToLink
      solver%numberOfLinkedSolvers=solver%numberOfLinkedSolvers+1
    ELSE IF(solver%numberOfLinkedSolvers>0.AND.solver%numberOfLinkedSolvers<=SOLVER_NUMBER_OF_SOLVER_TYPES) THEN
      !yes, there are already linked solvers
      !check if a solver of the same type has already been linked
      DO solverIdx=1,solver%numberOfLinkedSolvers
        IF(solver%linkedSolvers(solverIdx)%ptr%solveType==solverTypeToLink) THEN
          localError="The solver has already a linked solver of type "// &
            & TRIM(NumberToVString(solverTypeToLink,"*",err,error))//" attached to it."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !solverIdx
      ALLOCATE(oldLinkedSolvers(solver%numberOfLinkedSolvers),STAT=err)
      IF(err/=0) CALL FlagError("Could not old linked solvers.",err,error,*999)
      DO solverIdx=1,solver%numberOfLinkedSolvers
        oldLinkedSolvers(solverIdx)%ptr=>solver%linkedSolvers(solverIdx)%ptr
      ENDDO
      DEALLOCATE(solver%linkedSolvers)
      ALLOCATE(solver%linkedSolvers(solver%numberOfLinkedSolvers+1),STAT=err)
      IF(err/=0) CALL FlagError("Could not new linked solvers.",err,error,*999)
      DO solverIdx=1,solver%numberOfLinkedSolvers
        solver%linkedSolvers(solverIdx)%ptr=>oldLinkedSolvers(solverIdx)%ptr
      ENDDO
      solver%linkedSolvers(solver%numberOfLinkedSolvers+1)%ptr=>solverToLink
      solver%linkedSolverTypeMap(solverTypeToLink)%ptr=>solverToLink
      solver%numberOfLinkedSolvers=solver%numberOfLinkedSolvers+1
      DEALLOCATE(oldLinkedSolvers)
    ELSE
      localError="The number of linked solvers is "//TRIM(NumberToVString(solver%numberOfLinkedSolvers,"*",err,error))// &
        & " but should be between 0 and "//TRIM(NumberToVString(SOLVER_NUMBER_OF_SOLVER_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    !set the solver type for the linked solver
    solver%linkedSolverTypeMap(solverTypeToLink)%ptr%solveType=solverTypeToLink
    !set the linking solver for the linked solver
    solver%linkedSolverTypeMap(solverTypeToLink)%ptr%linkingSolver=>solver

    EXITS("Solver_LinkedSolverAdd")
    RETURN
999 ERRORSEXITS("Solver_LinkedSolverAdd",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinkedSolverAdd

  !
  !================================================================================================================================
  !

  !>Adds a linked solver to the solver. Also sets the solver type for the linked solver, als well as its linking solver.
  SUBROUTINE Solver_LinkedSolverRemove(solver,solveTypeToRemove,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver to add the linked solver to.
    INTEGER(INTG), INTENT(IN) :: solveTypeToRemove !<The solver type of the solver to be linked.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: solverIdx

    ENTERS("Solver_LinkedSolverRemove",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    IF(solveTypeToRemove<1.OR.solveTypeToRemove>SOLVER_NUMBER_OF_SOLVER_TYPES) THEN
      localError="The specified solver type is "//TRIM(NumberToVString(solveTypeToRemove,"*",err,error))//&
        & " but should be between 1 and "//TRIM(NumberToVString(SOLVER_NUMBER_OF_SOLVER_TYPES,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Check if there is any linked solvers
    IF(solver%numberOfLinkedSolvers>0.AND.solver%numberOfLinkedSolvers<=SOLVER_NUMBER_OF_SOLVER_TYPES) THEN
      !Check if a solver of the same type has already been linked
      DO solverIdx=1,solver%numberOfLinkedSolvers
        IF(solver%linkedSolvers(solverIdx)%ptr%solveType==solveTypeToRemove) THEN
          DEALLOCATE(solver%linkedSolvers)
          solver%numberOfLinkedSolvers=solver%numberOfLinkedSolvers-1
        ENDIF
      ENDDO !solverIdx
    ENDIF

    EXITS("Solver_LinkedSolverRemove")
    RETURN
999 ERRORSEXITS("Solver_LinkedSolverRemove",err,error)
    RETURN 1
    
  END SUBROUTINE Solver_LinkedSolverRemove

  !
  !================================================================================================================================
  !

END MODULE SolverRoutines

!
!================================================================================================================================
!

!>Called from the PETSc TS solvers to monitor the dynamic solver
SUBROUTINE Solver_TimeSteppingMonitorPETSc(ts,steps,time,X,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE ISO_VARYING_STRING
  USE Kinds
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE

  !Argument variables
  TYPE(PetscTSType), INTENT(INOUT) :: ts !<The PETSc ts type
  INTEGER(INTG), INTENT(INOUT) :: steps !<The iteration number
  REAL(DP), INTENT(INOUT) :: time !<The current time
  TYPE(PetscVecType), INTENT(INOUT) :: X !<The current iterate
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(DAESolverType), POINTER :: daeSolver
  TYPE(VARYING_STRING) :: error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*999)
  IF(ctx%solveType/=SOLVER_DAE_TYPE) THEN
    localError="Invalid solve type. The solve type of "//TRIM(NumberToVString(ctx%solveType,"*",err,error))// &
      & " does not correspond to a differntial-algebraic equations solver."
    CALL FlagError(localError,err,error,*999)
  ENDIF
  
  daeSolver=>ctx%DAESolver

  CALL SolverDAE_TimeSteppingMonitor(daeSolver,steps,time,err,error,*999)

  RETURN
999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error monitoring differential-algebraic equations solve.",err,error,*997)
997 RETURN
  
END SUBROUTINE Solver_TimeSteppingMonitorPETSc

!
!================================================================================================================================
!

!>Called from the PETSc SNES solvers to monitor the Newton nonlinear solver
SUBROUTINE Solver_NonlinearMonitorPETSc(snes,its,norm,ctx,err)

  USE BaseRoutines
  USE CmissPetscTypes
  USE ISO_VARYING_STRING
  USE Kinds
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Strings
  USE Types

  IMPLICIT NONE

  !Argument variables
  TYPE(PetscSnesType), INTENT(INOUT) :: snes !<The PETSc SNES type
  INTEGER(INTG), INTENT(INOUT) :: its !<The iteration number
  REAL(DP), INTENT(INOUT) :: norm !<The residual norm
  TYPE(SolverType), POINTER :: ctx !<The passed through context
  INTEGER(INTG), INTENT(INOUT) :: err !<The error code
  !Local Variables
  TYPE(NonlinearSolverType), POINTER :: nonlinearSolver
  TYPE(VARYING_STRING) :: error,localError

  IF(.NOT.ASSOCIATED(ctx)) CALL FlagError("Solver context is not associated.",err,error,*999)
  IF(ctx%solveType/=SOLVER_NONLINEAR_TYPE) THEN
    localError="Invalid solve type. The solve type of "//TRIM(NumberToVString(ctx%solveType,"*",err,error))// &
      & " does not correspond to a nonlinear solver."
    CALL FlagError(localError,err,error,*999)
  ENDIF
  
  nonlinearSolver=>ctx%nonlinearSolver

  CALL SolverNonlinear_Monitor(nonlinearSolver,its,norm,err,error,*999)

  RETURN
999 CALL WriteError(err,error,*998)
998 CALL FlagWarning("Error monitoring nonlinear solve.",err,error,*997)
997 RETURN
  
END SUBROUTINE Solver_NonlinearMonitorPETSc
