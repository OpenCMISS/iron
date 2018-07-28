!>This module handles all Navier-Stokes fluid routines.
MODULE NAVIER_STOKES_EQUATIONS_ROUTINES

  USE ADVECTION_EQUATION_ROUTINES
  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CHARACTERISTIC_EQUATION_ROUTINES
  USE CmissMPI
  USE CmissPetsc
  USE CmissPetscTypes
  USE ComputationEnvironment
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE COORDINATE_ROUTINES
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lapack
  USE Maths
  USE MatrixVector
  USE MESH_ROUTINES
  USE MeshAccessRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE STREE_EQUATION_ROUTINES
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PUBLIC NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE

  PUBLIC NavierStokes_EquationsSetSpecificationSet

  PUBLIC NavierStokes_EquationsSetSolutionMethodSet

  PUBLIC NAVIER_STOKES_EQUATIONS_SET_SETUP

  PUBLIC NavierStokes_PreSolveALEUpdateParameters,NavierStokes_PreSolveUpdateBoundaryConditions, &
    & NavierStokes_PreSolveALEUpdateMesh

  PUBLIC NAVIER_STOKES_PRE_SOLVE, NAVIER_STOKES_POST_SOLVE

  PUBLIC NavierStokes_ProblemSpecificationSet

  PUBLIC NAVIER_STOKES_PROBLEM_SETUP

  PUBLIC NavierStokes_FiniteElementResidualEvaluate,NavierStokes_FiniteElementJacobianEvaluate

  PUBLIC NavierStokes_BoundaryConditionsAnalyticCalculate

  PUBLIC NavierStokes_ResidualBasedStabilisation

  PUBLIC NavierStokes_Couple1D0D

  PUBLIC NavierStokes_CoupleCharacteristics

  PUBLIC NavierStokes_FiniteElementPreResidualEvaluate

  PUBLIC NavierStokes_ControlLoopPostLoop

  PUBLIC NavierStokes_UpdateMultiscaleBoundary

  PUBLIC NavierStokes_WallShearStressCalculate

CONTAINS

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Navier-Stokes flow equation type of an fluid mechanics equations set class.
  SUBROUTINE NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_NODAL_SOLUTION_METHOD
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
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes flow equation type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("NavierStokes_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_EquationsSetSolutionMethodSet

!
!================================================================================================================================
!

  !>Sets the equation specification for a Navier-Stokes fluid type of a fluid mechanics equations set class.
  SUBROUTINE NavierStokes_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("NavierStokes_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      ENDIF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_Coupled1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE)
        CALL FlagError("Not implemented yet.",err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      ENDIF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("NavierStokes_EquationsSetSpecificationSet",err,error)
    EXITS("NavierStokes_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE NavierStokes_EquationsSetSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Navier-Stokes fluid setup.
  SUBROUTINE NAVIER_STOKES_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD,ANALYTIC_FIELD,dependentField,geometricField
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: GEOMETRIC_SCALING_TYPE,GEOMETRIC_MESH_COMPONENT,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_ANALYTIC_COMPONENTS,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,GEOMETRIC_COMPONENT_NUMBER,I,componentIdx,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: MATERIAL_FIELD_NUMBER_OF_VARIABLES,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
    INTEGER(INTG) :: elementBasedComponents,nodeBasedComponents,constantBasedComponents
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS

    ENTERS("NAVIER_STOKES_EQUATIONS_SET_SETUP",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)
    NULLIFY(EQUATIONS_SET_FIELD_FIELD)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              CALL EquationsSet_LabelSet(EQUATIONS_SET,"Navier-Stokes equations set",err,error,*999)
           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              CALL EquationsSet_LabelSet(EQUATIONS_SET,"Navier-Stokes equations set",err,error,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 1
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Penalty Coefficient",err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                !Default the penalty coefficient value to 1E4
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                 & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1.0E4_DP,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL NavierStokes_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
              CALL EquationsSet_LabelSet(EQUATIONS_SET,"Navier-Stokes equations set",err,error,*999)
              EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 3
              nodeBasedComponents = 2  ! boundary flow, pressure
              elementBasedComponents = 11  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1, coupledNodeNumber
              constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field field for SUPG element metrics
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL FIELD_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                 CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_SET_FIELD_FIELD, &
                   & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                 CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & "BoundaryFlow",err,error,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_V_VARIABLE_TYPE, &
                   & "ElementMetrics",err,error,*999)
                 CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET_FIELD_FIELD,FIELD_U1_VARIABLE_TYPE, &
                   & "EquationsConstants",err,error,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_V_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,FIELD_U1_VARIABLE_TYPE, &
                   & FIELD_DP_TYPE,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_U_VARIABLE_TYPE,nodeBasedComponents,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_V_VARIABLE_TYPE,elementBasedComponents,err,error,*999)
                 CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                   & FIELD_U1_VARIABLE_TYPE,constantBasedComponents,err,error,*999)
              ELSE
                localError="User-specified fields are not yet implemented for an equations set field field &
                  & setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                  & SETUP_TYPE,"*",err,error))// " for a Navier-Stokes fluid."
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                !Default the Element Metrics parameter values 0.0
                nodeBasedComponents = 2  ! boundary flow,pressure
                elementBasedComponents = 11  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1, coupledNodeNumber
                constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
                ! Init boundary flux to 0
                DO componentIdx=1,nodeBasedComponents
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,ERR,ERROR,*999)
                END DO
                ! Init Element Metrics to 0 (except C1)
                DO componentIdx=1,elementBasedComponents-1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                END DO
                ! Default C1 to -1 for now, will be calculated in ResidualBasedStabilisation if not specified by user
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elementBasedComponents,-1.0_DP,err,error,*999)
                ! Boundary stabilisation scale factor (beta): default to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,0.0_DP,ERR,ERROR,*999)
                ! Max Courant (CFL) number: default to 0.0 (do not use)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,0.0_DP,ERR,ERROR,*999)
                ! Init Time increment to 0
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
                ! Stabilisation type: default to 1 for RBS (0=none, 2=RBVM)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)
              ELSE
                localError="User-specified fields are not yet implemented for an equations set field field &
                  & setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                  & SETUP_TYPE,"*",err,error))// " for a Navier-Stokes fluid."
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE, &
                & "*",err,error))// " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP% &
                & SETUP_TYPE,"*",err,error))// " is not implemented for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            !Do nothing???
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 1
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                DO componentIdx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              nodeBasedComponents = 2  ! boundary flow, pressure
              elementBasedComponents = 11  ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, C1, coupledNodeNumber
              constantBasedComponents = 4  ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET_FIELD_FIELD,&
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                ! Node-based fields
                DO componentIdx = 1, nodeBasedComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                END DO
                ! Element-based fields
                DO componentIdx = 1, elementBasedComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,componentIdx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                END DO
                ! Constant fields: boundary stabilisation scale factor and max courant #
                DO componentIdx = 1, constantBasedComponents
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
              ELSE
                !Do nothing
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              ! do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !-----------------------------------------------------------------
          ! D e p e n d e n t   f i e l d
          !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  DEPENDENT_FIELD_NUMBER_OF_VARIABLES=3
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_W_VARIABLE_TYPE],err,error,*999)
                ELSE
                  DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                ENDIF
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "del U/del n",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension and one for pressure
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO componentIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO !componentIdx
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_W_VARIABLE_TYPE, &
                    & "Wss",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_W_VARIABLE_TYPE, &
                    & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_W_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_W_VARIABLE_TYPE,1,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_W_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_W_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                !Other solutions not defined yet
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
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  DEPENDENT_FIELD_NUMBER_OF_VARIABLES=3
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_W_VARIABLE_TYPE],err,error,*999)
                ELSE
                  DEPENDENT_FIELD_NUMBER_OF_VARIABLES=2
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                ENDIF
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension and one for pressure
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_W_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_W_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_W_VARIABLE_TYPE,1,err,error,*999)
                ENDIF
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                      & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_W_VARIABLE_TYPE, &
                      & 1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
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
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                 & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                 & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                 & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
                  DO componentIdx=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_PRESSURE_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                  END DO
                END IF
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE, &
                & "*",err,error))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
               EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !set number of variables to 5 (U,DELUDELN,V,U1,U2)
              DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
              !calculate number of components (Q,A) for U and dUdN
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 6 (U,DELUDELN,V,U1,U2)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                !set data type
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)

                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
               !calculate number of components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                ! 2 component (W1,W2) for V
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & FIELD_U1_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                        & FIELD_U2_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  END DO
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
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
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field- Characteristic equations
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                  & err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)

                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components (Q,A) for U and dUdN
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                ! 2 component (W1,W2) for V
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE, &
                & "*",err,error))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                !start field creation with name 'DEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                !define new created field to be dependent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 5 (U,DELUDELN,V,U1,U2)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE, &
                    & FIELD_U3_VARIABLE_TYPE],err,error,*999)
                ELSE
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                    & err,error,*999)
                END IF
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components (Q,A)
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U1_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U1_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U1_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U2_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U2_VARIABLE_TYPE,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
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
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !set number of variables to 5 (U,DELUDELN,V,U1,U2)
                DEPENDENT_FIELD_NUMBER_OF_VARIABLES=5
                !Check the user specified field- Characteristic equations
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,DEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components (Q,A) for U and dUdN
                DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE, &
                  & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                    & "*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE, &
                & "*",err,error))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% &
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO !componentIdx
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                !Other solutions not defined yet
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !set number of variables to 1
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              !normalDirection for wave relative to node for W1,W2
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              !Create the auto created independent field
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error, &
                  & *999)
                !define new created field to be independent
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% &
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                !calculate number of components with one component for each dimension
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                  DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field- Characteristic equation
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                  & err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(localError,err,error,*999)
            END SELECT
         CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
           & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
           SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
           !Set start action
           CASE(EQUATIONS_SET_SETUP_START_ACTION)
             !set number of variables to 1
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              !normalDirection for wave relative to node for W1,W2
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              !Create the auto created independent field
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                ! Do nothing? independent field should be set up by characteristic equation routines
              ELSE
                !Check the user specified field- Characteristic equation
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES, &
                 & ERR,ERROR,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & ERR,ERROR,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,ERR,ERROR,*999)
              END IF
            !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,ERR,ERROR,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,EQUATIONS_SET% &
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO !componentIdx
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                  !Default geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                !Other solutions not defined yet
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                !calculate number of components with one component for each dimension
                INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                   DO componentIdx=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO !componentIdx
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END IF
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Navier-Stokes fluid"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
              !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
              IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                  dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                  IF(ASSOCIATED(dependentField)) THEN
                    EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                    IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                      IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                        geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                        IF(ASSOCIATED(geometricField)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 2 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            END IF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
                             & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                            !Check that this is a 1D equations set
                            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
                              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                              !Set analytic function type
                              EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                              !Set numbrer of components- Q,A (same as N-S depenedent field)
                              NUMBER_OF_ANALYTIC_COMPONENTS=2
                            ELSE
                              localError="The third equations set specification must by a TRANSIENT1D or COUPLED1D0D "// &
                                & "to use an analytic function of type "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))//"."
                              CALL FlagError(localError,err,error,*999)
                            END IF
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                            !Check that domain is 2D/3D
                            IF(NUMBER_OF_DIMENSIONS<2 .OR. NUMBER_OF_DIMENSIONS>3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 2 or 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            END IF
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            !Set numbrer of components
                            NUMBER_OF_ANALYTIC_COMPONENTS=10
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 2 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            END IF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5
                          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
                            !Set analtyic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for an analytic Navier-Stokes problem."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                          !Create analytic field if required
                          IF(NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                            IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                              !Create the auto created analytic field
                              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                                & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                              CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",err,error,*999)
                              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                                & err,error,*999)
                              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                                & err,error,*999)
                              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                                & GEOMETRIC_DECOMPOSITION,err,error,*999)
                              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                                & GEOMETRIC_FIELD,err,error,*999)
                              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,err,error,*999)
                              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                                & err,error,*999)
                              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & "Analytic",err,error,*999)
                              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_DP_TYPE,err,error,*999)
                              !Set the number of analytic components
                              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                              !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                              DO componentIdx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & componentIdx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                                IF(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE == &
                                 & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                    & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                                ELSE
                                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                    & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                                END IF
                              END DO !componentIdx
                              !Default the field scaling to that of the geometric field
                              CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                                & err,error,*999)
                              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                            ELSE
                              !Check the user specified field
                              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                              IF(NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                              ELSE
                                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                              END IF
                              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                                & err,error,*999)
                              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                            END IF
                          END IF
                        ELSE
                          CALL FlagError("Equations set materials is not finished.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set materials is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations analytic is not associated.",err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
              IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
                IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                  IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                    !Finish creating the analytic field
                    CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                    !Set the default values for the analytic field
                    SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                        !Default the analytic parameter values (L, H, U_mean, Pout) to 0.0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,4,0.0_DP,err,error,*999)
                      CASE DEFAULT
                        localError="The analytic function type of "// &
                          & TRIM(NumberToVString(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                          & " is invalid for an analytical static Navier-Stokes equation."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
                        !Default the analytic parameter values (U_characteristic, L) to 0.0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                        !Default the analytic parameter values to 0
                        NUMBER_OF_ANALYTIC_COMPONENTS = 10
                        DO componentIdx = 1,NUMBER_OF_ANALYTIC_COMPONENTS
                          CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                        END DO
                      CASE DEFAULT
                        localError="The analytic function type of "// &
                          & TRIM(NumberToVString(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                          & " is invalid for an analytical transient Navier-Stokes equation."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                      SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                        !Default the analytic parameter period values to 0
                        CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                      CASE DEFAULT
                        localError="The analytic function type of "// &
                          & TRIM(NumberToVString(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                          & " is invalid for a 1D Navier-Stokes equation."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT
                    CASE DEFAULT
                      localError="The third equations set specification of "// &
                        & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                        & " is invalid for an analytical Navier-Stokes equation set."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  END IF
                END IF
              ELSE
                CALL FlagError("Equations set analytic is not associated.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for an analytic Navier-Stokes problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=1
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=2! viscosity, density
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    &[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & "Materials",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials field
                  ! viscosity,density=1
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", &
                & err,error))//" is invalid for Navier-Stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=2
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=2! U_var (constant)  : viscosity scale, density
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS2=2! V_var (gaussBased): viscosity, shear rate
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"MaterialsField",err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    &[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                  ! Set up U_VARIABLE (constants)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                     & "MaterialsConstants",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  DO componentIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  END DO
                  ! Set up V_VARIABLE (gauss-point based, CellML in/out parameters)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                     & "ConstitutiveValues",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                  DO componentIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & componentIdx,GEOMETRIC_COMPONENT_NUMBER,ERR,ERROR,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  END DO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD, &
                    & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                  ! Check the U_VARIABLE
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,err,error,*999)
                  ! Check the U_VARIABLE
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials constants (viscosity scale, density)
                  DO componentIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  END DO
                  !Set the default values for the materials consitutive parameters (viscosity scale, density)
                  DO componentIdx=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  END DO
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", &
                & err,error))//" is invalid for Navier-Stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! 1 variables for the 1D Navier-Stokes materials
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=2
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS1=8
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS2=3
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            !Specify start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  !start field creation with name 'MATERIAL_FIELD'
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                    & EQUATIONS_SET%MATERIALS%MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  !label the field
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD, &
                    & GEOMETRIC_DECOMPOSITION,err,error,*999)
                  !point new field to geometric field
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  ! 2 U,V materials field
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    &[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  ! Set up Navier-Stokes materials parameters
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                    & FIELD_V_VARIABLE_TYPE,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  DO I=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS1 !(MU,RHO,alpha,pressureExternal,LengthScale,TimeScale,MassScale)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & I,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  END DO
                  DO I=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS2 !(A0,E,H0)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                  ! Set up coupling materials parameters
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                    & err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE], &
                    & err,error,*999)
                  ! Check N-S field variable
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS1,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & MATERIAL_FIELD_NUMBER_OF_COMPONENTS2,err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
              !Specify start action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", &
                & err,error))//" is invalid for Navier-Stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            !\todo: Think about gravity
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !Do nothing
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
              !? Maybe set finished flag????
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*", &
                & err,error))//" for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*", &
                & err,error))//" is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              &  " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              &  " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                ! Use the analytic Jacobian calculation
                !CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                !  & err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                ! Use the analytic Jacobian calculation
                !CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                !  & err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Navier-stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                ! Use the analytic Jacobian calculation
                !CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED],err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE, &
                    & MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Navier-Stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
              IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                  CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
                  CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
                  CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
                ELSE
                  CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Equations materials is not associated.",err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the creation of the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                ! Use the analytic Jacobian calculation
                !CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                !  & err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD, &
                  & "*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Navier-Stokes equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a Navier-Stokes fluid subtype."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("NAVIER_STOKES_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_EQUATIONS_SET_SETUP",err,error)
    RETURN 1

  END SUBROUTINE NAVIER_STOKES_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem pre solve.
  SUBROUTINE NAVIER_STOKES_PRE_SOLVE(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(NONLINEAR_SOLVER_TYPE), POINTER :: nonlinearSolver
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVER_MATRICES_TYPE), POINTER :: SOLVER_MATRICES
    TYPE(SOLVER_MATRIX_TYPE), POINTER :: SOLVER_MATRIX
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2,cellmlSolver
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: solver_matrix_idx,iteration,equationsSetIdx
    REAL(DP) :: timeIncrement,currentTime

    NULLIFY(SOLVER2)

    ENTERS("NAVIER_STOKES_PRE_SOLVE",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          !Since we can have a fluid mechanics navier stokes equations set in a coupled problem setup we do not necessarily
          !have PROBLEM%SPECIFICATION(1)==FLUID_MECHANICS_CLASS
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  ! TODO: Set up for multiple equations sets
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
                    IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
                      !Update boundary conditions and any analytic values
                      CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver mapping is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations is not associated.",err,error,*999)
              END IF
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
              !Update transient boundary conditions and any analytic values
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
            CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
               & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
              SELECT CASE(SOLVER%SOLVE_TYPE)
              ! --- D y n a m i c    S o l v e r s ---
              CASE(SOLVER_DYNAMIC_TYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equationsSetIdx = 1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetIdx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        ! --- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
                        CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                          IF (CONTROL_LOOP%PROBLEM%SPECIFICATION(3) == PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
                            IF (CONTROL_LOOP%PARENT_LOOP%WHILE_LOOP%ITERATION_NUMBER == 1) THEN
                              ! Only update fixed BCs once per timestep
                              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                            END IF
                          ELSE
                            !Update boundary conditions and any analytic values
                            CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                          END IF
                        ! --- 1 D    N a v i e r - S t o k e s   E q u a t i o n s ---
                        CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                           & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
                             & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
                             & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
                          ELSE
                            CALL FlagError("Dependent field is not associated.",err,error,*999)
                          END IF
                          !Update boundary conditions and any analytic values
                          CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,ERR,ERROR,*999)
                        ! --- A d v e c t i o n   S o l v e r ---
                        CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
                          CALL Advection_PreSolve(solver,err,error,*999)
                        CASE DEFAULT
                          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*", &
                            & err,error))//" is not valid for a nonlinear Navier-Stokes solver."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    END DO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                END IF
              ! --- N o n l i n e a r    S o l v e r s ---
              CASE(SOLVER_NONLINEAR_TYPE)
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,ERR,ERROR,*999)
                iteration = CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equationsSetIdx = 1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetIdx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        ! --- C h a r a c t e r i s t i c   E q u a t i o n s ---
                        CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          IF(ASSOCIATED(dependentField)) THEN
                            CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                             & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,err,error,*999)
                            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
                             & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,err,error,*999)
                            IF(iteration == 1) THEN
                              CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                              ! Extrapolate new W from Q,A if this is the first timestep
                              ! (otherwise will be calculated based on Navier-Stokes values)
                              CALL Characteristic_Extrapolate(SOLVER,currentTime,timeIncrement,err,error,*999)
                            END IF
                          ELSE
                            CALL FlagError("Dependent field is not associated.",err,error,*999)
                          END IF
                        CASE DEFAULT
                          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*", &
                            & err,error))//" is not valid for a nonlinear Navier-Stokes solver."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    END DO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                END IF
              ! --- C e l l M L   S o l v e r ---
              CASE(SOLVER_DAE_TYPE)
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,ERR,ERROR,*999)
                CALL SOLVER_DAE_TIMES_SET(SOLVER,currentTime,currentTime + timeIncrement,err,error,*999)
                CALL SOLVER_DAE_TIME_STEP_SET(SOLVER,timeIncrement/1000.0_DP,err,error,*999)
              ! --- L i n e a r   S o l v e r s ---
              CASE(SOLVER_LINEAR_TYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equationsSetIdx = 1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetIdx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        ! --- S t r u c t u r e d   T r e e  E q u a t i o n s ---
                        CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE, &
                           & EQUATIONS_SET_STREE1D0D_ADV_SUBTYPE)
                          CALL Stree_PRE_SOLVE(SOLVER,err,error,*999)
                        CASE DEFAULT
                          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*", &
                            & err,error))//" is not valid for a linear Navier-Stokes solver."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    END DO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="The solve type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",err,error))// &
                  & " is invalid for a multiscale Navier-Stokes problem type."
                CALL FlagError(localError,err,error,*999)
              END SELECT

            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
              &  PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              &  PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
              &  PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
              &  PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
              &  PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)

              SELECT CASE(SOLVER%SOLVE_TYPE)
              ! This switch takes advantage of the uniqueness of the solver types to do pre-solve operations
              ! for each of solvers in the various possible 1D subloops

              ! --- C h a r a c t e r i s t i c   S o l v e r ---
              CASE(SOLVER_NONLINEAR_TYPE)
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,err,error,*999)
                iteration = CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
                EQUATIONS_SET=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                ! Characteristic solver effectively solves for the mass/momentum conserving fluxes at the
                ! *NEXT* timestep by extrapolating current field values and then solving a system of nonlinear
                ! equations: cons mass, continuity of pressure, and the characteristics.
                NULLIFY(fieldVariable)
                CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_INPUT_DATA1_SET_TYPE)%ptr)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
                END IF
                CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                 & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,err,error,*999)
                CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
                 & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,err,error,*999)

                IF(iteration == 1) THEN
                  NULLIFY(fieldVariable)
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                  END IF
                  ! Extrapolate new W from Q,A if this is the first timestep (otherwise will be calculated based on Navier-Stokes
                  ! values)
                  CALL Characteristic_Extrapolate(SOLVER,currentTime,timeIncrement,err,error,*999)
                END IF

              ! --- 1 D   N a v i e r - S t o k e s   S o l v e r ---
              CASE(SOLVER_DYNAMIC_TYPE)
                IF(SOLVER%global_number==2) THEN
                  ! update solver matrix
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                      SOLVER_MATRICES=>SOLVER_equations%SOLVER_MATRICES
                      IF(ASSOCIATED(SOLVER_MATRICES)) THEN
                        DO solver_matrix_idx=1,SOLVER_MAPPING%NUMBER_OF_SOLVER_MATRICES
                          SOLVER_MATRIX=>SOLVER_MATRICES%MATRICES(solver_matrix_idx)%ptr
                          IF(ASSOCIATED(SOLVER_MATRIX)) THEN
                            SOLVER_MATRIX%UPDATE_MATRIX=.TRUE.
                          ELSE
                            CALL FlagError("Solver Matrix is not associated.",err,error,*999)
                          END IF
                        END DO
                      ELSE
                        CALL FlagError("Solver Matrices is not associated.",err,error,*999)
                      END IF
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                        IF(ASSOCIATED(dependentField)) THEN
                          CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
                           & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                          CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
                           & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
                        ELSE
                          CALL FlagError("Dependent field is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Solver mapping is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations is not associated.",err,error,*999)
                  END IF
                ELSE
                  ! --- A d v e c t i o n   S o l v e r ---
                  CALL Advection_PreSolve(solver,err,error,*999)
                END IF
                ! Update boundary conditions
                CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)

              ! --- C e l l M L    S o l v e r ---
              CASE(SOLVER_DAE_TYPE)
                ! DAE solver-set time
                CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,currentTime,timeIncrement,err,error,*999)
                CALL SOLVER_DAE_TIMES_SET(SOLVER,currentTime,currentTime + timeIncrement,err,error,*999)
                CALL SOLVER_DAE_TIME_STEP_SET(SOLVER,timeIncrement/1000.0_DP,err,error,*999)

              ! --- S T R E E    S o l v e r ---
              CASE(SOLVER_LINEAR_TYPE)
                CALL Stree_PRE_SOLVE(SOLVER,err,error,*999)

              CASE DEFAULT
                localError="The solve type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
                  & " is invalid for a 1D Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT

            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
              !First update mesh and calculates boundary velocity values
              CALL NavierStokes_PreSolveALEUpdateMesh(SOLVER,err,error,*999)
              !Then apply both normal and moving mesh boundary conditions
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                !Update boundary conditions for mesh-movement
                CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,err,error,*999)
                IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                  SOLVER2%DYNAMIC_SOLVER%ALE=.FALSE.
                ELSE
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                END IF
                !Update material properties for Laplace mesh movement
                CALL NavierStokes_PreSolveALEUpdateParameters(SOLVER,err,error,*999)
                !Pre solve for the linear solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                IF(SOLVER%DYNAMIC_SOLVER%ALE) THEN
                  !First update mesh and calculates boundary velocity values
                  CALL NavierStokes_PreSolveALEUpdateMesh(SOLVER,err,error,*999)
                  !Then apply both normal and moving mesh boundary conditions
                  CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
                ELSE
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver type is not associated for ALE problem.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
                !Pre solve for the linear solver
                IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  !TODO if first time step smooth imported mesh with respect to absolute nodal position?
                  !Update boundary conditions for mesh-movement
                  CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
                  IF(CONTROL_LOOP%problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                    & CONTROL_LOOP%problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                    CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,err,error,*999)
                  ELSE
                    CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,4,SOLVER2,err,error,*999)
                  ENDIF
                  IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                    SOLVER2%DYNAMIC_SOLVER%ALE=.FALSE.
                  ELSE
                    CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                  END IF
                  !Update material properties for Laplace mesh movement
                  CALL NavierStokes_PreSolveALEUpdateParameters(SOLVER,err,error,*999)
                  !Pre solve for the dynamic solver which deals with the coupled FiniteElasticity-NavierStokes problem
                ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  IF(SOLVER%DYNAMIC_SOLVER%ALE) THEN
                    !Apply both normal and moving mesh boundary conditions
                    CALL NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*999)
                  ELSE
                    CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver type is not associated for ALE problem.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
                CALL FlagError(localError,Err,Error,*999)
              END SELECT
            CASE DEFAULT
              localError="Problem type "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",err,error))// &
                & " is not valid for NAVIER_STOKES_PRE_SOLVE of a multi physics problem class."
              CALL FlagError(localError,Err,Error,*999)
            END SELECT
          CASE DEFAULT
            localError="Problem class "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(1),"*",err,error))// &
              & " is not valid for Navier-Stokes fluid types."
            CALL FlagError(localError,Err,Error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers are not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    END IF

    EXITS("NAVIER_STOKES_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_PRE_SOLVE",err,error)
    RETURN 1

  END SUBROUTINE NAVIER_STOKES_PRE_SOLVE

!
!================================================================================================================================
!

  !>Sets/changes the problem subtype for a Navier-Stokes fluid type.
  SUBROUTINE NavierStokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("NavierStokes_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
          !All ok
        CASE(PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE)
          CALL FlagError("Not implemented yet.",err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Navier-Stokes fluid mechanics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_NAVIER_STOKES_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Navier-Stokes problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_ProblemSpecificationSet")
    RETURN
999 ERRORS("NavierStokes_ProblemSpecificationSet",err,error)
    EXITS("NavierStokes_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE NavierStokes_ProblemSpecificationSet

!
!================================================================================================================================
!

  !>Sets up the Navier-Stokes problem.
  SUBROUTINE NAVIER_STOKES_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Navier-Stokes fluid on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CELLML_EQUATIONS_TYPE), POINTER :: CELLML_EQUATIONS,cellMLEquations
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: iterativeWhileLoop,iterativeWhileLoop2,iterativeWhileLoop3,simpleLoop
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,MESH_SOLVER_EQUATIONS,BIF_SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(SOLVER_TYPE), POINTER :: SOLVER, MESH_SOLVER,BIF_SOLVER,cellmlSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("NAVIER_STOKES_PROBLEM_SETUP",err,error,*999)

    NULLIFY(BIF_SOLVER)
    NULLIFY(BIF_SOLVER_EQUATIONS)
    NULLIFY(cellmlSolver)
    NULLIFY(CELLML_EQUATIONS)
    NULLIFY(CONTROL_LOOP)
    NULLIFY(CONTROL_LOOP_ROOT)
    NULLIFY(MESH_SOLVER)
    NULLIFY(MESH_SOLVER_EQUATIONS)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
        !All steady state cases of Navier-Stokes
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Transient cases and moving mesh
      CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !Set the first solver to be an CellML Evaluator for time varying boundary conditions
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes boundary condition CellML evaluation solver",err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !Set the second solver to be a first order dynamic solver
            NULLIFY(solver)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes dynamic nonlinear solver",err,error,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            !setup CellML evaluator for constitutive law
            IF(PROBLEM%SPECIFICATION(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              !Create the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_EVALUATOR_CREATE(SOLVER,cellmlSolver,ERR,ERROR,*999)
              !Link the CellML evaluator solver to the solver
              CALL SOLVER_LINKED_SOLVER_ADD(SOLVER,cellmlSolver,SOLVER_CELLML_EVALUATOR_TYPE,ERR,ERROR,*999)
            END IF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
              & err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            NULLIFY(solver)
            !Get the boundary condition solver
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            !Create the CellML equations
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            IF(PROBLEM%specification(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              !Get the multiscale solver
              NULLIFY(solver)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,cellmlSolver,err,error,*999)
              !Create the CellML equations
              CALL CELLML_EQUATIONS_CREATE_START(cellmlSolver,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            ENDIF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Get the CellML boundary condition solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
            IF(PROBLEM%specification(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              !Get the multiscale CellML solver
              NULLIFY(solver)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              !Get the CellML evaluator solver
              CALL SOLVER_NEWTON_CELLML_SOLVER_GET(SOLVER,cellmlSolver,err,error,*999)
              !Get the CellML equations for the CellML evaluator solver
              CALL SOLVER_CELLML_EQUATIONS_GET(cellmlSolver,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a CellML setup for a  transient Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ! Multiscale: 3D/1D/0D coupled
      CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
            NULLIFY(iterativeWhileLoop)
            ! The 3D-1D boundary value iterative coupling loop
            CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,ERR,ERROR,*999)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
            CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,100,ERR,ERROR,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E-6_DP,err,error,*999)
            CALL ControlLoop_RelativeToleranceSet(iterativeWhileLoop,0.001_DP,err,error,*999)
            IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"3D-0D Iterative Loop",ERR,ERROR,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              ! The 3D Navier-Stokes simple loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"3D Navier-Stokes",ERR,ERROR,*999)
            ELSE
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"3D-1D Iterative Loop",ERR,ERROR,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,ERR,ERROR,*999)
              NULLIFY(iterativeWhileLoop2)
              ! The 1D-0D boundary value iterative coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,100,ERR,ERROR,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,0.1_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D-0D Iterative Coupling Convergence Loop",ERR,ERROR,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop2,2,ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",ERR,ERROR,*999)
              NULLIFY(iterativeWhileLoop3)
              ! The Characteristics branch solver/ Navier-Stokes coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,2,iterativeWhileLoop3,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop3,PROBLEM_CONTROL_WHILE_LOOP_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop3,1000,ERR,ERROR,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop3,1.0E10_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop3,"1D Iterative Loop",ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              ! The 3D Navier-Stokes simple loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,ERR,ERROR,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"3D Navier-Stokes",ERR,ERROR,*999)
            END IF

          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            ! The 3D-1D iterative coupling loop
            NULLIFY(iterativeWhileLoop)
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
            IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(solver,"DAE Solver",ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Navier-Stokes 3D Solver",ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            ELSE
              ! Iterative loop 2 couples 1D and 0D, checking convergence
              ! (this subloop holds 2 subloops)
              NULLIFY(iterativeWhileLoop2)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)

              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,1,simpleLoop,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,ERR,ERROR,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,ERR,ERROR,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(solver,"DAE Solver",ERR,ERROR,*999)

              ! Iterative loop 3 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop3)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,2,iterativeWhileLoop3,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop3,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_NUMBER_SET(SOLVERS,2,ERR,ERROR,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Characteristic Solver",ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,ERR,ERROR,*999)
!!!-- 1 D   N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes 1D Solver",ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Navier-Stokes 3D Solver",ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
            END IF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            NULLIFY(iterativeWhileLoop)
            NULLIFY(iterativeWhileLoop2)
            NULLIFY(SOLVERS)
            IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              !Finish the 0D solvers
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              !Finish the 3D solver
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
            ELSE
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)
              !Finish the 0D solvers
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              !Finish the 1D solvers
              NULLIFY(iterativeWhileLoop3)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,2,iterativeWhileLoop3,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop3,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              !Finish the 3D solver
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a transient Navier-Stokes fluid."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          !Get the control loop
          IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            NULLIFY(iterativeWhileLoop)
            NULLIFY(simpleLoop)
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            ! 3D0D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
            ! 3D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              !Create the solver equations
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                & ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              !Finish the solver equations creation
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          ELSE
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            NULLIFY(iterativeWhileLoop)
            NULLIFY(iterativeWhileLoop2)
            NULLIFY(iterativeWhileLoop3)
            NULLIFY(simpleLoop)
            NULLIFY(SOLVERS)
            NULLIFY(SOLVER)
            ! 3D1D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
            ! 1D0D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)
            ! 1D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,2,iterativeWhileLoop3,ERR,ERROR,*999)
            ! 3D subloop
            CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,simpleLoop,ERR,ERROR,*999)
            SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              ! 1D NS/C subloop
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop3,SOLVERS,ERR,ERROR,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 1 D  N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC, &
                & ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              !Create the solver equations
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
                & ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              ! 1D NS/C subloop
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop3,SOLVERS,ERR,ERROR,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 1 D  N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- 3 D  N A V I E R   S T O K E S --!!!
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
              !Finish the solver equations creation
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a Navier-Stokes fluid."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          END IF

          !Create the CELLML solver equations
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
              ! 3D-1D loop
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              ! 0D loop
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(cellMLSolver)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL CELLML_EQUATIONS_CREATE_START(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
            ELSE
              ! 3D-1D loop
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              ! 1D-0D loop
              NULLIFY(iterativeWhileLoop2)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)
              ! 0D loop
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              NULLIFY(cellMLSolver)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL CELLML_EQUATIONS_CREATE_START(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
            END IF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            IF (PROBLEM%SPECIFICATION(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              ! 3D-0D
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              ! 0D
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
            ELSE
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
              ! 3D-1D
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,ERR,ERROR,*999)
              ! 1D-0D
              NULLIFY(iterativeWhileLoop2)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,iterativeWhileLoop2,ERR,ERROR,*999)
              ! 0D
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,1,simpleLoop,ERR,ERROR,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,ERR,ERROR,*999)
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(solver,CELLML_EQUATIONS,ERR,ERROR,*999)
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
              & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT

        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",ERR,ERROR))// &
            & " is invalid for a 3D-1d-0D Navier-Stokes fluid type."
          CALL FlagError(localError,ERR,ERROR,*999)
        END SELECT

      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &    !1D Navier-Stokes
        & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &    !  with coupled 0D boundaries
        & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, & !  with coupled advection
        & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, & !  with coupled 0D boundaries and advection
        & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &      !  with stree
        & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)     !  with stree and advection

        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for Coupled1dDaeNavierStokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            NULLIFY(CONTROL_LOOP_ROOT)
            !Time Loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The 1D-0D boundary value iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                & err,error,*999)
            ELSE IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The 1D-0D boundary value iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                & err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple Advection solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",err,error,*999)
            ELSE IF(PROBLEM%specification(3) == PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E3_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
            ELSE IF(PROBLEM%specification(3) == PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E6_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple Advection solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",err,error,*999)
            ELSE IF(PROBLEM%specification(3) == PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The 1D-0D boundary value iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,1,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                & err,error,*999)
            ELSE IF(PROBLEM%specification(3) == PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              ! The 1D-0D boundary value iterative coupling loop
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(CONTROL_LOOP,2,err,error,*999)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
              CALL CONTROL_LOOP_NUMBER_OF_SUB_LOOPS_SET(iterativeWhileLoop,2,err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple CellML solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"0D CellML solver Loop",err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(iterativeWhileLoop2,PROBLEM_CONTROL_WHILE_LOOP_TYPE,err,error,*999)
              CALL CONTROL_LOOP_MAXIMUM_ITERATIONS_SET(iterativeWhileLoop2,1000,err,error,*999)
              CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
                & err,error,*999)
              NULLIFY(simpleLoop)
              ! The simple Advection solver loop
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_TYPE_SET(simpleLoop,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
              CALL CONTROL_LOOP_LABEL_SET(simpleLoop,"Advection",err,error,*999)
            END IF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a 1d transient Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Create the solvers
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          NULLIFY(CONTROL_LOOP)
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              ! Simple loop 1 contains the Advection solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- A D V E C T I O N --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Advection Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)

              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"DAE Solver",err,error,*999)

              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)

              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- A D V E C T I O N --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Advection Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Linear Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- A D V E C T I O N --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Advection Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_LINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)

              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"DAE Solver",err,error,*999)

              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)

              ! Simple loop 1 contains the 0D/CellML DAE solver
              ! (this subloop holds 1 solver)
              NULLIFY(simpleLoop)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL SOLVERS_CREATE_START(simpleLoop,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,1,err,error,*999)
!!!-- D A E --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(solvers,1,solver,err,error,*999)
              CALL SOLVER_TYPE_SET(solver,SOLVER_LINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Linear Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)

              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL SOLVERS_CREATE_START(iterativeWhileLoop2,solvers,err,error,*999)
              CALL SOLVERS_NUMBER_SET(solvers,2,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Characteristic Solver",err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
!!!-- N A V I E R   S T O K E S --!!!
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(solver,"Navier-Stokes Solver",err,error,*999)
              CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            IF(PROBLEM%specification(3)==PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            ELSE IF(PROBLEM%specification(3)==PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            ELSE IF(PROBLEM%specification(3)==PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            ELSE IF(PROBLEM%specification(3)==PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            ELSE IF(PROBLEM%specification(3)==PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            ELSE IF(PROBLEM%specification(3)==PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
              NULLIFY(simpleLoop)
              NULLIFY(SOLVERS)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
            END IF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a 1d transient Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Create the solver equations
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- D A E --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- D A E --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
              CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            NULLIFY(CONTROL_LOOP)
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- D A E --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(simpleLoop)
              ! Iterative loop couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,2,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
!!!-- A D V E C T I O N --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
              NULLIFY(iterativeWhileLoop)
              ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
              ! (this subloop holds 2 subloops)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- D A E --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER_EQUATIONS)
              NULLIFY(iterativeWhileLoop2)
              ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
              ! (this subloop holds 2 solvers)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(iterativeWhileLoop2,SOLVERS,err,error,*999)
!!!-- C H A R A C T E R I S T I C --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
              NULLIFY(SOLVER)
              NULLIFY(SOLVER_EQUATIONS)
!!!-- N A V I E R   S T O K E S --!!!
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
              CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          !Create the CELLML solver equations
        CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
            ELSE
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            END IF
            NULLIFY(SOLVER)
            NULLIFY(cellMLSolver)
            NULLIFY(CELLML_EQUATIONS)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_START(solver,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup Navier-Stokes equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            IF(PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & PROBLEM%specification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              NULLIFY(iterativeWhileLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(CONTROL_LOOP,1,iterativeWhileLoop,err,error,*999)
              NULLIFY(simpleLoop)
              CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop,1,simpleLoop,err,error,*999)
              CALL CONTROL_LOOP_SOLVERS_GET(simpleLoop,SOLVERS,err,error,*999)
            ELSE
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            END IF
            NULLIFY(SOLVER)
            SELECT CASE(PROBLEM%specification(3))
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(solver,CELLML_EQUATIONS,err,error,*999)
              CALL CELLML_EQUATIONS_CREATE_FINISH(CELLML_EQUATIONS,err,error,*999)
            CASE DEFAULT
              localError="The third problem specification of "// &
                & TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for cellML equations setup Navier-Stokes fluid mechanics problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a 1d transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
        !Quasi-static Navier-Stokes
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a quasistatic Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a quasistatic Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Navier-Stokes ALE cases
      CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a ALE Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a ALE Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,2,err,error,*999)
            !Set the first solver to be a linear solver for the Laplace mesh movement problem
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(MESH_SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(MESH_SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            !Set the solver to be a first order dynamic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a ALE Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(MESH_SOLVER,MESH_SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(MESH_SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
              & err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,MESH_SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(MESH_SOLVER,MESH_SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(MESH_SOLVER_EQUATIONS,err,error,*999)
            !Get the solver equations
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a ALE Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third problem specification of "//TRIM(NumberToVString(PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid mechanics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("NAVIER_STOKES_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_PROBLEM_SETUP",err,error)
    RETURN 1

  END SUBROUTINE NAVIER_STOKES_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    INTEGER(INTG) :: ng,mh,mhs,mi,ms,nh,nhs,ni,ns,nhs_max,mhs_max,nhs_min,mhs_min,xv,out
    INTEGER(INTG) :: FIELD_VAR_TYPE,MESH_COMPONENT1,MESH_COMPONENT2,MESH_COMPONENT_NUMBER
    INTEGER(INTG) :: nodeIdx,xiIdx,coordIdx,derivativeIdx,versionIdx,elementVersionNumber,componentIdx
    INTEGER(INTG) :: numberOfVersions,nodeNumber,numberOfElementNodes,numberOfParameters,firstNode,lastNode
    REAL(DP) :: JGW,SUM,X(3),DXI_DX(3,3),DPHIMS_DXI(3),DPHINS_DXI(3),PHIMS,PHINS,momentum,mass,QUpwind,AUpwind,pExternal
    REAL(DP) :: U_VALUE(3),W_VALUE(3),U_DERIV(3,3),Q_VALUE,A_VALUE,Q_DERIV,A_DERIV,area,pressure,normalWave,normal,Lref,Tref,Mref
    REAL(DP) :: MU_PARAM,RHO_PARAM,A0_PARAM,E_PARAM,H_PARAM,A0_DERIV,E_DERIV,H_DERIV,alpha,beta,kappa,G0_PARAM,muScale
    REAL(DP), POINTER :: dependentParameters(:),materialsParameters(:),materialsParameters1(:)
    LOGICAL :: updateStiffnessMatrix,updateDampingMatrix,updateRHSVector,updateNonlinearResidual
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_FiniteElementResidualEvaluate",err,error,*999)

    updateStiffnessMatrix=.FALSE.
    updateDampingMatrix=.FALSE.
    updateRHSVector=.FALSE.
    updateNonlinearResidual=.FALSE.
    X=0.0_DP
    out=0

    NULLIFY(DEPENDENT_BASIS,GEOMETRIC_BASIS)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(linearMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(dynamicMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(linearMatrices)
    NULLIFY(nonlinearMatrices)
    NULLIFY(dynamicMatrices)
    NULLIFY(rhsVector)
    NULLIFY(stiffnessMatrix, dampingMatrix)
    NULLIFY(dependentField,independentField,geometricField,materialsField)
    NULLIFY(dependentParameters,materialsParameters,materialsParameters1)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(QUADRATURE_SCHEME)
    NULLIFY(QUADRATURE_SCHEME1, QUADRATURE_SCHEME2)
    NULLIFY(DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          !Set general and specific pointers
          dependentField=>equations%interpolation%dependentField
          independentField=>equations%interpolation%independentField
          geometricField=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          vectorMatrices=>vectorEquations%vectorMatrices
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            linearMatrices=>vectorMatrices%linearMatrices
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            stiffnessMatrix=>linearMatrices%MATRICES(1)%ptr
            linearMapping=>vectorMapping%linearMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
            linearMatrices=>vectorMatrices%linearMatrices
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            stiffnessMatrix=>linearMatrices%MATRICES(1)%ptr
            linearMapping=>vectorMapping%linearMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%MATRICES(1)%ptr
            dampingMatrix=>dynamicMatrices%MATRICES(2)%ptr
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%MATRICES(1)%ptr
            dampingMatrix=>dynamicMatrices%MATRICES(2)%ptr
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)

          CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%MATRICES(1)%ptr
            dampingMatrix=>dynamicMatrices%MATRICES(2)%ptr
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            independentField=>equations%interpolation%independentField
            INDEPENDENT_BASIS=>independentField%DECOMPOSITION%DOMAIN(independentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
              & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            stiffnessMatrix=>dynamicMatrices%MATRICES(1)%ptr
            dampingMatrix=>dynamicMatrices%MATRICES(2)%ptr
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            dampingMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
            IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
            IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
            IF(ASSOCIATED(nonlinearMatrices)) updateNonlinearResidual=nonlinearMatrices%updateResidual
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          !Loop over Gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              W_VALUE(1)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              W_VALUE(2)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                W_VALUE(3)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
              END IF
            ELSE
              W_VALUE=0.0_DP
            END IF

            ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
              ! Note the constant from the U_VARIABLE is a scale factor
              muScale = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              ! Get the gauss point based value returned from the CellML solver
              CALL Field_ParameterSetGetLocalGaussPoint(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ng,ELEMENT_NUMBER,1,MU_PARAM,err,error,*999)
              MU_PARAM=MU_PARAM*muScale
            ELSE
              MU_PARAM = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
            END IF
            RHO_PARAM = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)

            !Start with matrix calculations
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              !Loop over field components
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF(updateStiffnessMatrix.OR.updateDampingMatrix) THEN
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        &(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      ! JGW=equations%interpolation%geometricInterpPointMetrics%JACOBIAN*QUADRATURE_SCHEME2%&
                      ! &GAUSS_WEIGHTS(ng)
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        !Calculate some general values
                        DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                          DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DXI_DX(mi,ni)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & DXI_DX(mi,ni)
                          END DO
                          DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        END DO !ni
                        PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        !Laplace only matrix
                        IF(updateStiffnessMatrix) THEN
                          !LAPLACE TYPE
                          IF(nh==mh) THEN
                            SUM=0.0_DP
                            !Calculate SUM
                            DO xv=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                  SUM=SUM+MU_PARAM*DPHINS_DXI(ni)*DXI_DX(ni,xv)*DPHIMS_DXI(mi)*DXI_DX(mi,xv)
                                END DO !ni
                              END DO !mi
                            END DO !x
                            !Calculate MATRIX
                            stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                        END IF
                        !General matrix
                        IF(updateStiffnessMatrix) THEN
                          !GRADIENT TRANSPOSE TYPE
                          IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE) THEN
                            IF(nh<FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                  !note mh/nh derivative in DXI_DX
                                  SUM=SUM+MU_PARAM*DPHINS_DXI(mi)*DXI_DX(mi,mh)*DPHIMS_DXI(ni)*DXI_DX(ni,nh)
                                END DO !ni
                              END DO !mi
                              !Calculate MATRIX
                              stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs) &
                               & +SUM*JGW
                           END IF
                          END IF
                        END IF
                        !Contribution through ALE
                        IF(updateStiffnessMatrix) THEN
                          !GRADIENT TRANSPOSE TYPE
                          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
                            IF(nh==mh) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                                  SUM=SUM-RHO_PARAM*W_VALUE(mi)*DPHINS_DXI(ni)*DXI_DX(ni,mi)*PHIMS
                                END DO !ni
                              END DO !mi
                              !Calculate MATRIX
                              stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+ &
                               & SUM*JGW
                            END IF
                          END IF
                        END IF
                        !Pressure contribution (B transpose)
                        IF(updateStiffnessMatrix) THEN
                          !LAPLACE TYPE
                          IF(nh==FIELD_VARIABLE%NUMBER_OF_COMPONENTS) THEN
                            SUM=0.0_DP
                            !Calculate SUM
                            DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              SUM=SUM-PHINS*DPHIMS_DXI(ni)*DXI_DX(ni,mh)
                            END DO !ni
                            !Calculate MATRIX
                            stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                        END IF
                        !Damping matrix
                        IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE .OR. &
                          & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                          IF(updateDampingMatrix) THEN
                            IF(nh==mh) THEN
                              SUM=0.0_DP
                              !Calculate SUM
                              SUM=PHIMS*PHINS*RHO_PARAM
                              !Calculate MATRIX
                              dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                            END IF
                          END IF
                        END IF
                      END DO !ns
                    END DO !nh
                  END IF
                END DO !ms
              END DO !mh
              !Analytic RHS vector
              IF(rhsVector%firstAssembly) THEN
                IF(updateRHSVector) THEN
                  IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                    IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                      & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                      mhs=0
                      DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                        MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                          & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                        DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                          mhs=mhs+1
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          !note mh value derivative
                          SUM=0.0_DP
                          X(1) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)
                          X(2) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
                          IF(DEPENDENT_BASIS1%NUMBER_OF_XI==3) THEN
                            X(3) = equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,1)
                          END IF
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1) THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(X(1)**3*RHO_PARAM+3.0_DP*MU_PARAM*10.0_DP**2- &
                                & 3.0_DP*RHO_PARAM*X(2)**2*X(1))/(10.0_DP**4))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2) &
                            & THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-4.0_DP*MU_PARAM/10.0_DP/10.0_DP*EXP((X(1)-X(2))/10.0_DP))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3) &
                           & THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(16.0_DP*MU_PARAM*PI**2/10.0_DP**2*COS(2.0_DP*PI*X(2)/10.0_DP)* &
                                & COS(2.0_DP*PI*X(1)/10.0_DP)- &
                                & 2.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)*RHO_PARAM*PI/10.0_DP)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4) &
                           & THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=PHIMS*(2.0_DP*SIN(X(1))*COS(X(2)))*MU_PARAM
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-2.0_DP*COS(X(1))*SIN(X(2)))*MU_PARAM
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5) &
                           & THEN
                            !do nothing
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(RHO_PARAM*X(1)**3+6.0_DP*RHO_PARAM*X(1)*X(3)*X(2)+ &
                                & 6.0_DP*MU_PARAM*10.0_DP**2- &
                                & 3.0_DP*RHO_PARAM*X(2)**2*X(1)-3.0_DP*RHO_PARAM*X(3)*X(1)**2-3.0_DP*RHO_PARAM*X(3)*X(2)**2)/ &
                                  & (10.0_DP**4))
                            ELSE IF(mh==3) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-2.0_DP/3.0_DP*(6.0_DP*RHO_PARAM*X(1)*X(3)*X(2)+RHO_PARAM*X(1)**3+ &
                                & 6.0_DP*MU_PARAM*10.0_DP**2- &
                                & 3.0_DP*RHO_PARAM*X(1)*X(3)**2-3.0_DP*RHO_PARAM*X(2)*X(1)**2-3.0_DP*RHO_PARAM*X(2)*X(3)**2)/ &
                                & (10.0_DP**4))
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2) THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*((-4.0_DP*MU_PARAM*EXP((X(1)-X(2))/10.0_DP)-2.0_DP*MU_PARAM*EXP((X(2)-X(3))/10.0_DP)+ &
                                & RHO_PARAM*EXP((X(3)-X(2))/10.0_DP)*10.0_DP)/10.0_DP**2)
                            ELSE IF(mh==3) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-(4.0_DP*MU_PARAM*EXP((X(3)-X(1))/10.0_DP)+2.0_DP*MU_PARAM*EXP((X(2)-X(3))/10.0_DP)+ &
                                & RHO_PARAM*EXP((X(3)-X(2))/10.0_DP)*10.0_DP)/10.0_DP** 2)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3) THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              SUM=0.0_DP
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              SUM=PHIMS*(2.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)*(18.0_DP*COS(2.0_DP*PI*X(1)/10.0_DP)* &
                                & MU_PARAM*PI*SIN(2.0_DP*PI*X(3)/10.0_DP)-3.0_DP*RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)**2* &
                                & SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP-2.0_DP*RHO_PARAM*SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP+ &
                                & 2.0_DP*RHO_PARAM*SIN(2.0_DP*PI*X(2)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(3)/10.0_DP)**2)*PI/ &
                                & 10.0_DP**2)
                            ELSE IF(mh==3) THEN
                              !Calculate SUM
                              SUM=PHIMS*(-2.0_DP*PI*COS(2.0_DP*PI*X(3)/10.0_DP)*RHO_PARAM*SIN(2.0_DP*PI*X(3)/10.0_DP)* &
                                & (-1.0_DP+COS(2.0_DP*PI*X(2)/10.0_DP)**2)/10.0_DP)
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4) THEN
                            IF(mh==1) THEN
                              !Calculate SUM
                              !SUM=PHIMS*(2.0_DP*SIN(X(1))*COS(X(2)))*MU_PARAM
                            ELSE IF(mh==2) THEN
                              !Calculate SUM
                              !SUM=PHIMS*(-2.0_DP*COS(X(1))*SIN(X(2)))*MU_PARAM
                            END IF
                          ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                            !do nothing
                          END IF
                          !Calculate RH VECTOR
                           rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)+SUM*JGW
                        END DO !ms
                      END DO !mh
                    ELSE
                      rhsVector%elementVector%vector(mhs)=0.0_DP
                    END IF
                  END IF
                END IF
              END IF

              !Calculate nonlinear vector
              IF(updateNonlinearResidual) THEN
                ! Get interpolated velocity and velocity gradient values for nonlinear term
                U_VALUE(1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
                U_VALUE(2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
                U_DERIV(1,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S1)
                U_DERIV(1,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S2)
                U_DERIV(2,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S1)
                U_DERIV(2,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S2)
                IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                  U_VALUE(3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
                  U_DERIV(3,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S1)
                  U_DERIV(3,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S2)
                  U_DERIV(3,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S3)
                  U_DERIV(1,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S3)
                  U_DERIV(2,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S3)
                ELSE
                  U_VALUE(3)=0.0_DP
                  U_DERIV(3,1)=0.0_DP
                  U_DERIV(3,2)=0.0_DP
                  U_DERIV(3,3)=0.0_DP
                  U_DERIV(1,3)=0.0_DP
                  U_DERIV(2,3)=0.0_DP
                END IF
                !Here W_VALUES must be ZERO if ALE part of linear matrix
                W_VALUE=0.0_DP
                mhs=0
                DO mh=1,(FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1)
                  MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                  JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                    & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                  DXI_DX=0.0_DP

                  DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                    DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      DXI_DX(mi,ni)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(mi,ni)
                    END DO
                  END DO

                  DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                    !note mh value derivative
                    SUM=0.0_DP
                    ! Convective form
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      SUM=SUM+RHO_PARAM*(PHIMS)*( &
                        & (U_VALUE(1))*(U_DERIV(mh,ni)*DXI_DX(ni,1))+ &
                        & (U_VALUE(2))*(U_DERIV(mh,ni)*DXI_DX(ni,2))+ &
                        & (U_VALUE(3))*(U_DERIV(mh,ni)*DXI_DX(ni,3)))
                    END DO !ni

                    nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)+SUM*JGW

                  END DO !ms
                END DO !mh
              END IF
            END IF

            !------------------------------------------------------------------
            ! R e s i d u a l - b a s e d    S t a b i l i s a t i o n
            !------------------------------------------------------------------
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
              CALL NavierStokes_ResidualBasedStabilisation(EQUATIONS_SET,ELEMENT_NUMBER,ng, &
               & MU_PARAM,RHO_PARAM,.FALSE.,err,error,*999)
            END IF

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!                                        !!!!!
            !!!!!         1 D  T R A N S I E N T         !!!!!
            !!!!!                                        !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Start with matrix calculations
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              Q_VALUE=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              Q_DERIV=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,FIRST_PART_DERIV)
              A_VALUE=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              A_DERIV=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,FIRST_PART_DERIV)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3, &
               & alpha,err,error,*999)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,8, &
               & G0_PARAM,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
              A0_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              A0_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(1,FIRST_PART_DERIV)
              E_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              E_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(2,FIRST_PART_DERIV)
              H_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
              H_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(3,FIRST_PART_DERIV)
              beta = (4.0_DP*(SQRT(PI))*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)  !(kg/m2/s2)
              kappa = 8.0_DP*PI*MU_PARAM/RHO_PARAM ! viscous resistance operator

              ! If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
              IF(A_VALUE < A0_PARAM*0.001_DP) THEN
                A_VALUE = A0_PARAM*0.001_DP
              END IF

              mhs=0
              !Loop Over Element Rows
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(mh)%DOMAIN%TOPOLOGY%ELEMENTS
                DXI_DX=0.0_DP
                !Calculate dxi_dx in 3D
                DO xiIdx=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                  DO coordIdx=1,equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE) &
                    & %ptr%NUMBER_OF_X_DIMENSIONS
                    DXI_DX(1,1)=DXI_DX(1,1)+(equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)% &
                      & PTR%DXI_DX(xiIdx,coordIdx))**2.0_DP
                  END DO !coordIdx
                END DO !xiIdx
                DXI_DX(1,1)=SQRT(DXI_DX(1,1))
                !Loop Over Element rows
                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  DPHIMS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,FIRST_PART_DERIV,ng)
                  mhs=mhs+1
                  nhs=0
                  IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                    !Loop Over Element Columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        DPHINS_DXI(1)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)
                        nhs=nhs+1

                        !!!-- D A M P I N G  M A T R I X --!!!
                        IF(updateDampingMatrix) THEN
                          !Momentum Equation, dQ/dt
                          IF(mh==1 .AND. nh==1) THEN
                            SUM=PHINS*PHIMS
                            dampingMatrix%elementMatrix%matrix(mhs,nhs)= &
                              & dampingMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                          !Mass Equation, dA/dt
                          IF(mh==2 .AND. nh==2) THEN
                            SUM=PHINS*PHIMS
                            dampingMatrix%elementMatrix%matrix(mhs,nhs)= &
                              & dampingMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                        END IF

                        !!!-- S T I F F N E S S  M A T R I X --!!!
                        IF(updateStiffnessMatrix) THEN
                          IF(mh==1 .AND. nh==2) THEN
                            !Momentum Equation, linearisable A0 terms
                            SUM=-PHINS*PHIMS*(beta*SQRT(A0_PARAM)/RHO_PARAM)*(H_DERIV/H_PARAM + E_DERIV/E_PARAM)*DXI_DX(1,1)
                            stiffnessMatrix%elementMatrix%matrix(mhs,nhs)= &
                              & stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                            !Momentum Equation, gravitational force
                            SUM=PHINS*PHIMS*G0_PARAM
                            stiffnessMatrix%elementMatrix%matrix(mhs,nhs)= &
                              & stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                          !Mass Equation, dQ/dX, flow derivative
                          IF(mh==2 .AND. nh==1) THEN
                            SUM=DPHINS_DXI(1)*DXI_DX(1,1)*PHIMS
                            stiffnessMatrix%elementMatrix%matrix(mhs,nhs)= &
                              & stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*JGW
                          END IF
                        END IF

                      END DO !ns
                    END DO !nh
                  END IF

                  !!!-- N O N L I N E A R  V E C T O R --!!!
                  IF(updateNonlinearResidual) THEN
                    !Momentum Equation
                    IF(mh==1) THEN
                      SUM=((2.0_DP*alpha*(Q_VALUE/A_VALUE)*Q_DERIV - &
                        & (alpha*((Q_VALUE/A_VALUE)**2.0_DP)*A_DERIV)+(beta/RHO_PARAM)* &           !Convective
                        & ((SQRT(A_VALUE)/2.0_DP)*A_DERIV+ &                                        !A  gradient
                        & (A_VALUE/(2.0_DP*SQRT(A0_PARAM))-(A_VALUE**1.5_DP)/A0_PARAM)*A0_DERIV+ &  !A0 gradient
                        & (A_VALUE*(SQRT(A_VALUE)))*(H_DERIV/H_PARAM) + &                           !H  gradient (nonlinear part)
                        & (A_VALUE*(SQRT(A_VALUE)))*(E_DERIV/E_PARAM)))* &                          !E  gradient (nonlinear part)
                        & DXI_DX(1,1)+kappa*(Q_VALUE/A_VALUE))*PHIMS                                !Viscosity
                      nonlinearMatrices%elementResidual%vector(mhs)= &
                        & nonlinearMatrices%elementResidual%vector(mhs)+SUM*JGW
                    END IF
                  END IF

                END DO !ms
              END DO !mh
            END IF
          END DO !ng

          IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            IF(updateNonlinearResidual) THEN
              ELEMENTS_TOPOLOGY=>dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%components(1)%domain%topology%elements
              numberOfElementNodes=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_NODES
              numberOfParameters=ELEMENTS_TOPOLOGY%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS
              firstNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(1)
              lastNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(numberOfElementNodes)
              !Get material constants
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
                & RHO_PARAM,err,error,*999)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4, &
                & pExternal,err,error,*999)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,5, &
                & Lref,err,error,*999)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,6, &
                & Tref,err,error,*999)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,7, &
                & Mref,err,error,*999)

              !!!-- P R E S S U R E    C A L C U L A T I O N --!!!
              !Loop over the element nodes and versions
              DO nodeIdx=1,numberOfElementNodes
                nodeNumber=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nodeIdx)
                derivativeIdx = 1
                versionIdx=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)% &
                 & elementVersions(derivativeIdx,nodeIdx)
                !Get current Area values
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                 & derivativeIdx,nodeNumber,2,area,err,error,*999)
                IF(area < A0_PARAM*0.001_DP) area = A0_PARAM*0.001_DP
                !Get material parameters
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                 & derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                 & derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
                  & derivativeIdx,nodeNumber,3,H_PARAM,err,error,*999)
                beta = (4.0_DP*(SQRT(PI))*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)  !(kg/m2/s2)
                !Pressure equation in mmHg
                pressure=(pExternal+beta*(SQRT(area)-SQRT(A0_PARAM)))!/(Mref/(Lref*Tref**2.0))!*0.0075_DP
                !Update the dependent field
                IF(ELEMENT_NUMBER<=dependentField%DECOMPOSITION%TOPOLOGY%ELEMENTS%NUMBER_OF_ELEMENTS) THEN
                  CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,1,nodeNumber,1,pressure,err,error,*999)
                END IF
              END DO

              !!!-- B R A N C H   F L U X   U P W I N D I N G --!!!
              !----------------------------------------------------
              ! In order to enforce conservation of mass and momentum across discontinuous
              ! branching topologies, flux is upwinded against the conservative branch values
              ! established by the characteristic solver.
              DO nodeIdx=1,numberOfElementNodes
                nodeNumber=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nodeIdx)
                numberOfVersions=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%NODES%NODES(nodeNumber)%DERIVATIVES(1)%numberOfVersions

                ! Find the branch node on this element
                IF(numberOfVersions>1) THEN
                  derivativeIdx = 1
                  elementVersionNumber=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)% &
                   & elementVersions(derivativeIdx,nodeIdx)

                  ! Find the wave direction - incoming or outgoing
                  DO componentIdx = 1,2
                    CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                     & elementVersionNumber,derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
                    IF(ABS(normalWave) > ZERO_TOLERANCE) THEN
                      normal = normalWave
                    END IF
                  END DO

                  ! Get materials parameters for node on this element
                  CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)
                  CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)
                  CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & elementVersionNumber,derivativeIdx,nodeNumber,3,H_PARAM,err,error,*999)
                  beta = (4.0_DP*(SQRT(PI))*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)

                  !Get current Q & A values for node on this element
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,1,Q_VALUE,err,error,*999)
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,2,A_VALUE,err,error,*999)

                  !Get upwind Q & A values based on the branch (characteristics) solver
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,1,QUpwind,err,error,*999)
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
                   & elementVersionNumber,derivativeIdx,nodeNumber,2,AUpwind,err,error,*999)

                  ! If A goes negative during nonlinear iteration, set to positive value to avoid segfault
                  IF(A_VALUE < A0_PARAM*0.001_DP) THEN
                    A_VALUE = A0_PARAM*0.001_DP
                  END IF

                  !Momentum Equation: F_upwind - F_Current
                  momentum = ((alpha*(QUpwind**2.0_DP)/AUpwind+(AUpwind**1.5_DP-A0_PARAM**1.5_DP)*(beta/(3.0_DP*RHO_PARAM))) &
                         & - (alpha*(Q_VALUE**2.0_DP)/A_VALUE+(A_VALUE**1.5_DP-A0_PARAM**1.5_DP)*(beta/(3.0_DP*RHO_PARAM))))*normal

                  !Continuity Equation
                  mass = (QUpwind-Q_VALUE)*normal

                  !Add momentum/mass contributions to first/last node accordingly
                  IF(nodeNumber==firstNode) THEN
                    nonlinearMatrices%elementResidual%vector(1)= &
                      & nonlinearMatrices%elementResidual%vector(1)+momentum
                    nonlinearMatrices%elementResidual%vector(numberOfParameters+1)= &
                      & nonlinearMatrices%elementResidual%vector(numberOfParameters+1)+mass
                  ELSE IF(nodeNumber==lastNode) THEN
                    nonlinearMatrices%elementResidual%vector(numberOfParameters)= &
                      & nonlinearMatrices%elementResidual%vector(numberOfParameters)+momentum
                    nonlinearMatrices%elementResidual%vector(numberOfParameters*2)= &
                      & nonlinearMatrices%elementResidual%vector(numberOfParameters*2)+mass
                  END IF
                END IF !version>1
              END DO !loop nodes
              ! Update any distributed pressure field values
              CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & err,error,*999)
              CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & err,error,*999)
            END IF
          END IF

          ! F a c e   I n t e g r a t i o n
          IF(updateNonlinearResidual) THEN
            !If specified, also perform boundary line (2D) or face (3D) integration for neumann boundary conditions
            CALL NavierStokes_FiniteElementBoundaryIntegrate(EQUATIONS_SET,ELEMENT_NUMBER,FIELD_VARIABLE,.FALSE.,err,error,*999)
          END IF

          !!!--   A S S E M B L E   M A T R I C E S  &  V E C T O R S   --!!!
          mhs_min=mhs
          mhs_max=nhs
          nhs_min=mhs
          nhs_max=nhs
          IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
            IF(stiffnessMatrix%firstAssembly) THEN
              IF(updateStiffnessMatrix) THEN
                DO mhs=mhs_min+1,mhs_max
                  DO nhs=1,nhs_min
                    !Transpose pressure type entries for mass equation
                    stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=-stiffnessMatrix%elementMatrix%matrix(nhs,mhs)
                  END DO
                END DO
              END IF
            END IF
          END IF

        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a Navier-Stokes equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementResidualEvaluate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS,INDEPENDENT_BASIS
    TYPE(DECOMPOSITION_TYPE), POINTER :: DECOMPOSITION
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: ELEMENTS_TOPOLOGY
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: ng,mh,mhs,mi,ms,nh,nhs,ni,ns,x,xiIdx,coordIdx
    INTEGER(INTG) :: derivativeIdx,elementVersionNumber,firstNode,lastNode,nodeIdx,nodeNumber
    INTEGER(INTG) :: numberOfElementNodes,numberOfParameters,numberOfVersions,componentIdx
    INTEGER(INTG) :: FIELD_VAR_TYPE,MESH_COMPONENT_NUMBER,MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: JGW,SUM,DXI_DX(3,3),DPHIMS_DXI(3),DPHINS_DXI(3),PHIMS,PHINS
    REAL(DP) :: U_VALUE(3),W_VALUE(3),U_DERIV(3,3),Q_VALUE,Q_DERIV,A_VALUE,A_DERIV,alpha,beta,normal,normalWave,kappa
    REAL(DP) :: MU_PARAM,RHO_PARAM,A0_PARAM,A0_DERIV,E_PARAM,E_DERIV,H_PARAM,H_DERIV,mass,momentum1,momentum2,muScale
    LOGICAL  :: updateJacobianMatrix

    ENTERS("NavierStokes_FiniteElementJacobianEvaluate",err,error,*999)

    DXI_DX=0.0_DP
    updateJacobianMatrix=.FALSE.

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
          & err,error,*999)
      END IF
      NULLIFY(EQUATIONS)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%specification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          !Set some general and case-specific pointers
          dependentField=>equations%interpolation%dependentField
          independentField=>equations%interpolation%independentField
          geometricField=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          vectorMatrices=>vectorEquations%vectorMatrices
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          vectorMapping=>vectorEquations%vectorMapping
          SELECT CASE(EQUATIONS_SET%specification(3))
          CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
            linearMatrices=>vectorMatrices%linearMatrices
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            stiffnessMatrix=>linearMatrices%MATRICES(1)%ptr
            linearMapping=>vectorMapping%linearMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
          CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
            linearMatrices=>vectorMatrices%linearMatrices
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            stiffnessMatrix=>linearMatrices%MATRICES(1)%ptr
            linearMapping=>vectorMapping%linearMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            !SOURCE_VECTOR=>vectorMatrices%SOURCE_VECTOR
            stiffnessMatrix%elementMatrix%matrix=0.0_DP
            nonlinearMatrices%elementResidual%vector=0.0_DP
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
          CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
            nonlinearMapping=>vectorMapping%nonlinearMapping
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            jacobianMatrix%elementJacobian%matrix=0.0_DP
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            linearMapping=>vectorMapping%linearMapping
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
          CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            nonlinearMapping=>vectorMapping%nonlinearMapping
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            jacobianMatrix%elementJacobian%matrix=0.0_DP
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            linearMapping=>vectorMapping%linearMapping
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
          CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
            nonlinearMapping=>vectorMapping%nonlinearMapping
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            jacobianMatrix%elementJacobian%matrix=0.0_DP
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            linearMapping=>vectorMapping%linearMapping
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
          CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            DECOMPOSITION => dependentField%DECOMPOSITION
            MESH_COMPONENT_NUMBER = DECOMPOSITION%MESH_COMPONENT_NUMBER
            nonlinearMapping=>vectorMapping%nonlinearMapping
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            jacobianMatrix%elementJacobian%matrix=0.0_DP
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            linearMapping=>vectorMapping%linearMapping
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
          CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE)
            independentField=>equations%interpolation%independentField
            INDEPENDENT_BASIS=>independentField%DECOMPOSITION%DOMAIN(independentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)% &
              & PTR%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            nonlinearMapping=>vectorMapping%nonlinearMapping
            nonlinearMatrices=>vectorMatrices%nonlinearMatrices
            jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
            jacobianMatrix%elementJacobian%matrix=0.0_DP
            dynamicMatrices=>vectorMatrices%dynamicMatrices
            dynamicMapping=>vectorMapping%dynamicMapping
            FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
            linearMapping=>vectorMapping%linearMapping
            IF(ASSOCIATED(jacobianMatrix)) updateJacobianMatrix=jacobianMatrix%updateJacobian
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_MESH_VELOCITY_SET_TYPE,ELEMENT_NUMBER,equations% &
              & interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%specification(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          !Loop over all Gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              W_VALUE(1)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              W_VALUE(2)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                W_VALUE(3)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
              END IF
            ELSE
              W_VALUE=0.0_DP
            END IF

            ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
              ! Note the constant from the U_VARIABLE is a scale factor
              muScale = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              ! Get the gauss point based value returned from the CellML solver
              CALL Field_ParameterSetGetLocalGaussPoint(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ng,ELEMENT_NUMBER,1,MU_PARAM,err,error,*999)
              MU_PARAM=MU_PARAM*muScale
            ELSE
              MU_PARAM = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
            END IF
            RHO_PARAM = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)

            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN

              U_VALUE(1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              U_VALUE(2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              U_DERIV(1,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S1)
              U_DERIV(1,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S2)
              U_DERIV(2,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S1)
              U_DERIV(2,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S2)
              IF(FIELD_VARIABLE%NUMBER_OF_COMPONENTS==4) THEN
                U_VALUE(3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
                U_DERIV(3,1)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S1)
                U_DERIV(3,2)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S2)
                U_DERIV(3,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(3,PART_DERIV_S3)
                U_DERIV(1,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,PART_DERIV_S3)
                U_DERIV(2,3)=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,PART_DERIV_S3)
              ELSE
                U_VALUE(3)=0.0_DP
                U_DERIV(3,1)=0.0_DP
                U_DERIV(3,2)=0.0_DP
                U_DERIV(3,3)=0.0_DP
                U_DERIV(1,3)=0.0_DP
                U_DERIV(2,3)=0.0_DP
              END IF
              !Start with calculation of partial matrices
              !Here W_VALUES must be ZERO if ALE part of linear matrix
              W_VALUE=0.0_DP
            END IF

            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE) THEN
              !Loop over field components
              mhs=0

              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)

                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF(updateJacobianMatrix) THEN
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS-1
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        &(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        !Calculate some general values needed below
                        DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                          DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DXI_DX(mi,ni)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & DXI_DX(mi,ni)
                          END DO
                          DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        END DO !ni
                        PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        SUM=0.0_DP
                        IF(updateJacobianMatrix) THEN
                          !Calculate J1 only
                          DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            SUM=SUM+(PHINS*U_DERIV(mh,ni)*DXI_DX(ni,nh)*PHIMS*RHO_PARAM)
                          END DO
                          !Calculate MATRIX
                          jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs) &
                             & +SUM*JGW
                          !Calculate J2 only
                          IF(nh==mh) THEN
                            SUM=0.0_DP
                            !Calculate SUM
                            DO x=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                              DO mi=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                                SUM=SUM+RHO_PARAM*(U_VALUE(x)-W_VALUE(x))*DPHINS_DXI(mi)*DXI_DX(mi,x)*PHIMS
                              END DO !mi
                            END DO !x
                            !Calculate MATRIX
                            jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs) &
                              & +SUM*JGW
                          END IF
                        END IF
                      END DO !ns
                    END DO !nh
                  END IF
                END DO !ms
              END DO !mh
              ! Stabilisation terms
              IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                & EQUATIONS_SET%specification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
                CALL NavierStokes_ResidualBasedStabilisation(EQUATIONS_SET,ELEMENT_NUMBER,ng,MU_PARAM,RHO_PARAM,.TRUE., &
                  & err,error,*999)
              END IF
            END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                        !!!!!
!!!!!         1 D  T R A N S I E N T         !!!!!
!!!!!                                        !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !Start with Matrix Calculations
            IF(EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
              & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
              Q_VALUE=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              Q_DERIV=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,FIRST_PART_DERIV)
              A_VALUE=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              A_DERIV=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(2,FIRST_PART_DERIV)
              CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3, &
                & alpha,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
              A0_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              A0_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(1,FIRST_PART_DERIV)
              E_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
              E_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(2,FIRST_PART_DERIV)
              H_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
              H_DERIV=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(3,FIRST_PART_DERIV)
              beta = (4.0_DP*SQRT(PI)*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)     !(kg/m2/s2)
              kappa = 8.0_DP*PI*MU_PARAM/RHO_PARAM ! viscous resistance operator

              ! If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
              IF(A_VALUE < A0_PARAM*0.001_DP) THEN
                A_VALUE = A0_PARAM*0.001_DP
              END IF

              mhs=0
              !Loop Over Element Rows
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                  & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
                ELEMENTS_TOPOLOGY=>FIELD_VARIABLE%COMPONENTS(mh)%DOMAIN%TOPOLOGY%ELEMENTS
                DXI_DX(1,1)=0.0_DP
                !Calculate dxi_dx in 3D
                DO xiIdx=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                  DO coordIdx=1,equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE) &
                    & %ptr%NUMBER_OF_X_DIMENSIONS
                    DXI_DX(1,1)=DXI_DX(1,1)+(equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)% &
                      & PTR%DXI_DX(xiIdx,coordIdx))**2.0_DP
                  END DO !coordIdx
                END DO !xiIdx
                DXI_DX(1,1)=SQRT(DXI_DX(1,1))
                !Loop Over Element rows
                DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                  PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  DPHIMS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,FIRST_PART_DERIV,ng)
                  mhs=mhs+1
                  nhs=0
                  !Loop Over Element Columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                    DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                      & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                    QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                      &(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                    DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                      PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                      DPHINS_DXI(1)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ns,FIRST_PART_DERIV,ng)
                      nhs=nhs+1
                      IF(updateJacobianMatrix) THEN

                        !Momentum Equation (dF/dQ)
                        IF(mh==1 .AND. nh==1) THEN
                          SUM=((alpha*2.0_DP*PHINS*Q_DERIV/A_VALUE +  &
                            & alpha*2.0_DP*Q_VALUE*DPHINS_DXI(1)/A_VALUE+ &
                            & (-2.0_DP)*alpha*Q_VALUE*PHINS*A_DERIV/(A_VALUE**2.0_DP))*DXI_DX(1,1)+ &   !Convective
                            & ((PHINS*kappa/A_VALUE)))*PHIMS                                                  !Viscosity
                          jacobianMatrix%elementJacobian%matrix(mhs,nhs)= &
                            & jacobianMatrix%elementJacobian%matrix(mhs,nhs)+SUM*JGW
                        END IF

                        !Momentum Equation (dF/dA)
                        IF(mh==1 .AND. nh==2) THEN
                          SUM=((((-2.0_DP*alpha*Q_VALUE*PHINS*Q_DERIV)/(A_VALUE**2.0_DP))+ &
                            & ((2.0_DP*alpha*PHINS*(Q_VALUE**2.0_DP)*A_DERIV)/(A_VALUE**3.0_DP))+ &
                            & (-alpha*((Q_VALUE/A_VALUE)**2.0_DP)*DPHINS_DXI(1))+ &                              !Convective
                            & ((0.5_DP*PHINS*(1.0_DP/SQRT(A_VALUE))*A_DERIV+SQRT(A_VALUE)*DPHINS_DXI(1))+ &      !Area Gradient
                            & ((1.0_DP/SQRT(A0_PARAM))-((3.0_DP/(A0_PARAM))*SQRT(A_VALUE)))*(A0_DERIV) + &       !Ref Area Gradient
                            & (2.0_DP*PHINS*1.5_DP*SQRT(A_VALUE))*H_DERIV/H_PARAM+ &                             !Thickness Gradient
                            & (2.0_DP*PHINS*1.5_DP*SQRT(A_VALUE))*E_DERIV/E_PARAM) &                             !Elasticity Gradient
                            & *beta/(2.0_DP*RHO_PARAM))*DXI_DX(1,1)+(-PHINS*kappa*Q_VALUE/A_VALUE**2.0_DP))*PHIMS!Viscosity
                          jacobianMatrix%elementJacobian%matrix(mhs,nhs)= &
                            & jacobianMatrix%elementJacobian%matrix(mhs,nhs)+SUM*JGW
                        END IF

                      END IF
                    END DO !ns
                  END DO !nh
                END DO !ms
              END DO !mh
            END IF
          END DO !ng

          ! B o u n d a r y   I n t e g r a t i o n
          SELECT CASE(EQUATIONS_SET%specification(3))
          CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
            & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
            ! Calculate the Jacobian of the nonlinear boundary stabilisation term if beta > 0
            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
              & FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,beta,err,error,*999)
            IF (beta > ZERO_TOLERANCE) THEN
              CALL NavierStokes_FiniteElementBoundaryIntegrate(EQUATIONS_SET,ELEMENT_NUMBER,FIELD_VARIABLE,.TRUE.,ERR,ERROR,*999)
            END IF
          CASE DEFAULT
            ! Do nothing for other equation set subtypes
          END SELECT

          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
            & EQUATIONS_SET%specification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            ELEMENTS_TOPOLOGY=>dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%components(1)%domain%topology%elements
            numberOfElementNodes=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_NODES
            numberOfParameters=ELEMENTS_TOPOLOGY%MAXIMUM_NUMBER_OF_ELEMENT_PARAMETERS
            firstNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(1)
            lastNode=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(numberOfElementNodes)
            !Get material constants
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
              & RHO_PARAM,err,error,*999)

!!!-- B R A N C H   F L U X   U P W I N D I N G --!!!
            !----------------------------------------------------
            ! In order to enforce conservation of mass and momentum across discontinuous
            ! branching topologies, flux is upwinded against the conservative branch values
            ! established by the characteristic solver.
            DO nodeIdx=1,numberOfElementNodes
              nodeNumber=ELEMENTS_TOPOLOGY%ELEMENTS(ELEMENT_NUMBER)%ELEMENT_NODES(nodeIdx)
              numberOfVersions=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%NODES%NODES(nodeNumber)%DERIVATIVES(1)%numberOfVersions

              ! Find the branch node on this element
              IF(numberOfVersions>1) THEN
                derivativeIdx = 1
                elementVersionNumber=ELEMENTS_TOPOLOGY%DOMAIN%TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)% &
                  & elementVersions(derivativeIdx,nodeIdx)

                ! Find the wave direction - incoming or outgoing
                DO componentIdx = 1,2
                  CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & elementVersionNumber,derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
                  IF(ABS(normalWave) > ZERO_TOLERANCE) THEN
                    normal = normalWave
                  END IF
                END DO

                ! Get materials parameters for node on this element
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementVersionNumber,derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementVersionNumber,derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementVersionNumber,derivativeIdx,nodeNumber,3,H_PARAM,err,error,*999)
                beta = (4.0_DP*(SQRT(PI))*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)

                !Get current Q & A values
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementVersionNumber,derivativeIdx,nodeNumber,1,Q_VALUE,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & elementVersionNumber,derivativeIdx,nodeNumber,2,A_VALUE,err,error,*999)

                !Momentum Equation, d/dQ
                momentum1 = (-alpha*2.0_DP*Q_VALUE/A_VALUE)*normal

                !Momentum Equation, d/dA
                momentum2 = (alpha*(Q_VALUE/A_VALUE)**2.0_DP-1.5_DP*(A_VALUE**0.5_DP)*(beta/(3.0_DP*RHO_PARAM)))*normal

                !Continuity Equation , d/dQ
                mass = -1.0_DP*normal

                !Add momentum/mass contributions to first/last node accordingly
                IF(nodeNumber==firstNode) THEN
                  jacobianMatrix%elementJacobian%matrix(1,1)= &
                    & jacobianMatrix%elementJacobian%matrix(1,1)+momentum1
                  jacobianMatrix%elementJacobian%matrix(1,numberOfParameters+1)= &
                    & jacobianMatrix%elementJacobian%matrix(1,numberOfParameters+1)+momentum2
                  jacobianMatrix%elementJacobian%matrix(numberOfParameters+1,1)= &
                    & jacobianMatrix%elementJacobian%matrix(numberOfParameters+1,1)+mass
                ELSE IF(nodeNumber==lastNode) THEN
                  jacobianMatrix%elementJacobian%matrix(numberOfParameters,numberOfParameters)= &
                    & jacobianMatrix%elementJacobian%matrix(numberOfParameters,numberOfParameters)+momentum1
                  jacobianMatrix%elementJacobian%matrix(numberOfParameters,2*numberOfParameters)= &
                    & jacobianMatrix%elementJacobian%matrix(numberOfParameters,2*numberOfParameters)+momentum2
                  jacobianMatrix%elementJacobian%matrix(2*numberOfParameters,numberOfParameters)= &
                    & jacobianMatrix%elementJacobian%matrix(2*numberOfParameters,numberOfParameters)+mass
                END IF
              END IF !version>1
            END DO !loop nodes

          END IF

        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementJacobianEvaluate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem post solve.
  SUBROUTINE NAVIER_STOKES_POST_SOLVE(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EquationsType), POINTER :: equations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: equationsSetIdx,iteration,timestep,outputIteration,equationsSetNumber
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement
    LOGICAL :: convergedFlag,fluidEquationsSetFound

    ENTERS("NAVIER_STOKES_POST_SOLVE",err,error,*999)
    NULLIFY(SOLVER2)
    NULLIFY(SOLVERS)
    NULLIFY(dependentField)
    NULLIFY(fieldVariable)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE,PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
              ENDIF
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
              SELECT CASE(SOLVER%SOLVE_TYPE)
              CASE(SOLVER_NONLINEAR_TYPE)
                ! Characteristic solver- copy branch Q,A values to new parameter set
                dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                END IF
                iteration = CONTROL_LOOP%WHILE_LOOP%ITERATION_NUMBER
                IF(iteration == 1) THEN
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                END IF
              CASE(SOLVER_DYNAMIC_TYPE)
                ! Navier-Stokes solver: do nothing
              CASE DEFAULT
                localError="The solver type of "//TRIM(NumberToVString(SOLVER%SOLVE_TYPE,"*",err,error))// &
                  & " is invalid for a 1D Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                  END IF
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)

                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,err,error,*999)
                  END IF
                CASE DEFAULT
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                IF(SOLVER%GLOBAL_NUMBER == 1) THEN
                  ! DAE solver- do nothing
                ELSE
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the CellML DAE simple loop of a 1D0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END IF
              ELSE
                localError="The control loop type for solver "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END IF
            CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                  END IF
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)

                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,err,error,*999)
                  END IF
                CASE DEFAULT
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                IF(SOLVER%GLOBAL_NUMBER == 1) THEN
                  ! DAE solver- do nothing
                ELSE
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the CellML DAE simple loop of a 1D0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END IF
              ELSE
                localError="The control loop type for solver "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END IF
            CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
              SELECT CASE(SOLVER%GLOBAL_NUMBER)
              CASE(1)
                ! Characteristic solver- copy branch Q,A values to new parameter set
                dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                END IF
                CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)
              CASE(2)
                ! check characteristic/ N-S convergence at branches
  !                CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,err,error,*999)
              CASE(3)
                ! Advection solver output data if necessary
                IF(CONTROL_LOOP%WHILE_LOOP%CONTINUE_LOOP .EQV. .FALSE.) THEN
                  ! 1D NSE solver output data if N-S/Chars converged
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              CASE DEFAULT
                localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                  & " is invalid for a 1D Navier-Stokes and Advection problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                  END IF
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,err,error,*999)
                  END IF
                CASE DEFAULT
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                ! DAE and advection solvers - output data if post advection solve
                IF(SOLVER%SOLVERS%CONTROL_LOOP%SUB_LOOP_INDEX == 3) THEN
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              ELSE
                localError="The control loop type for solver "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END IF
            CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_equations%SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr%DEPENDENT%DEPENDENT_FIELD
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%ptr)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
                  END IF
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,err,error,*999)
                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,err,error,*999)
                  END IF
                CASE DEFAULT
                  localError="The solver global number of "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                ! DAE and advection solvers - output data if post advection solve
                IF(SOLVER%SOLVERS%CONTROL_LOOP%SUB_LOOP_INDEX == 3) THEN
                  CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
                END IF
              ELSE
                localError="The control loop type for solver "//TRIM(NumberToVString(SOLVER%GLOBAL_NUMBER,"*",err,error))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END IF
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,startTime,stopTime,currentTime,timeIncrement, &
                  & timestep,outputIteration,err,error,*999)
                IF(ASSOCIATED(SOLVER%SOLVER_EQUATIONS)) THEN
                  convergedFlag = .FALSE.
                  CALL NavierStokes_CalculateBoundaryFlux3D0D(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                    & EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS%equationsSet,err,error,*999)
                ENDIF
                CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
              ENDIF
            CASE(PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                DO equationsSetNumber=1,SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  ! If this is a coupled constitutive (non-Newtonian) viscosity problem, update shear rate values
                  !  to be passed to the CellML solver at beginning of next timestep
                  IF(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(equationsSetNumber)%PTR% &
                    &  equations%equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                    CALL NavierStokes_ShearRateCalculate(SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING% &
                      & EQUATIONS_SETS(equationsSetNumber)%PTR%EQUATIONS%equationsSet,err,error,*999)
                  END IF
                END DO
              ENDIF
            CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%WHILE_LOOP)) THEN
                SELECT CASE(SOLVER%GLOBAL_NUMBER)
                CASE(1)
                  ! Characteristic solver- copy branch Q,A values to new parameter set
                  dependentField=>SOLVER%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR%DEPENDENT%DEPENDENT_FIELD
                  CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
                  IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_UPWIND_VALUES_SET_TYPE)%PTR)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
                     & FIELD_UPWIND_VALUES_SET_TYPE,ERR,ERROR,*999)
                  END IF
                  CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)

                CASE(2)
                  ! ! 1D Navier-Stokes solver
                  IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL==3) THEN
                    ! check characteristic/ N-S convergence at branches
  !                    CALL NavierStokes_CoupleCharacteristics(CONTROL_LOOP,SOLVER,ERR,ERROR,*999)
                  END IF
                CASE DEFAULT
                  localError="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,ERR,ERROR,*999)
                END SELECT
              ELSE IF(ASSOCIATED(SOLVER%SOLVERS%CONTROL_LOOP%SIMPLE_LOOP)) THEN
                IF(SOLVER%GLOBAL_NUMBER == 1) THEN
                  ! DAE solver- do nothing
                ELSE
                  localError="The solver global number of "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                    & " is invalid for the CellML DAE simple loop of a 1D0D coupled Navier-Stokes problem."
                  CALL FlagError(localError,ERR,ERROR,*999)
                END IF
              ELSE
                localError="The control loop type for solver "//TRIM(NUMBER_TO_VSTRING(SOLVER%GLOBAL_NUMBER,"*",ERR,ERROR))// &
                  & " is invalid for the a 1D0D coupled Navier-Stokes problem."
                CALL FlagError(localError,ERR,ERROR,*999)
              END IF
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Post solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER2,err,error,*999)
                IF(ASSOCIATED(SOLVER2%DYNAMIC_SOLVER)) THEN
                  SOLVER2%DYNAMIC_SOLVER%ALE=.TRUE.
                ELSE
                  CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
                END IF
              !Post solve for the dynamic solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                NULLIFY(solverEquations)
                CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
                NULLIFY(solverMapping)
                CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
                equationsSetIdx=1
                fluidEquationsSetFound=.FALSE.
                DO WHILE(.NOT.fluidEquationsSetFound.AND.equationsSetIdx<=solverMapping%NUMBER_OF_EQUATIONS_SETS)
                  equations=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(equationsSetIdx)%equations
                  NULLIFY(equationsSet)
                  CALL Equations_EquationsSetGet(equations,equationsSet,err,error,*999)
                  IF(equationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                    & .AND.equationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                    & .AND.(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
                    & .OR.equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
                    fluidEquationsSetFound=.TRUE.
                  ELSE
                    equationsSetIdx=equationsSetIdx+1
                  END IF
                END DO
                IF(.NOT.fluidEquationsSetFound) THEN
                  localError="Fluid equations set not found."
                  CALL FlagError(localError,Err,Error,*999)
                END IF
                !CALL NavierStokes_WallShearStressCalculate(equationsSet,err,error,*999)
                CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
              END IF
            CASE DEFAULT
              localError="The third problem specification of  "// &
                & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes fluid mechanics problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Control loop is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    END IF

    EXITS("NAVIER_STOKES_POST_SOLVE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_POST_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE NAVIER_STOKES_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Navier-Stokes flow pre solve
  SUBROUTINE NavierStokes_PreSolveUpdateBoundaryConditions(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET,SOLID_EQUATIONS_SET,FLUID_EQUATIONS_SET
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE), POINTER :: SOLID_DEPENDENT
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE), POINTER :: FLUID_GEOMETRIC
    TYPE(EquationsType), POINTER :: EQUATIONS,SOLID_EQUATIONS,FLUID_EQUATIONS
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: INTERPOLATED_POINT(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: INTERPOLATION_PARAMETERS(:)
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,geometricField,materialsField
    TYPE(FIELD_TYPE), POINTER :: independentField,SOLID_dependentField,FLUID_GEOMETRIC_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentFieldVariable,independentFieldVariable
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS,SOLID_SOLVER_EQUATIONS,FLUID_SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING,SOLID_SOLVER_MAPPING,FLUID_SOLVER_MAPPING
    TYPE(SOLVER_TYPE), POINTER :: Solver2
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeIdx,derivativeIdx,versionIdx,variableIdx,numberOfSourceTimesteps,timeIdx,componentIdx
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE,GLOBAL_DERIV_INDEX,node_idx,variable_type
    INTEGER(INTG) :: variable_idx,local_ny,ANALYTIC_FUNCTION_TYPE,component_idx,deriv_idx,dim_idx,version_idx
    INTEGER(INTG) :: element_idx,en_idx,I,J,K,number_of_nodes_xic(3),search_idx,localDof,globalDof,componentBC,previousNodeNumber
    INTEGER(INTG) :: componentNumberVelocity,numberOfDimensions,numberOfNodes,numberOfGlobalNodes
    INTEGER(INTG) :: dependentVariableType,independentVariableType,dependentDof,independentDof,userNodeNumber,localNodeNumber
    INTEGER(INTG) :: EquationsSetIndex,SolidNodeNumber,FluidNodeNumber,equationsSetIdx
    INTEGER(INTG) :: currentTimeLoopIteration,outputIterationNumber,numberOfFittedNodes,computationalNode
    INTEGER(INTG), ALLOCATABLE :: InletNodes(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,DISPLACEMENT_VALUE,VALUE,XI_COORDINATES(3),timeData,QP,QPP,componentValues(3)
    REAL(DP) :: T_COORDINATES(20,3),MU_PARAM,RHO_PARAM,X(3),FluidGFValue,SolidDFValue,NewLaplaceBoundaryValue,Lref,Tref,Mref
    REAL(DP) :: startTime,stopTime,currentTime,timeIncrement
    REAL(DP), POINTER :: MESH_VELOCITY_VALUES(:), GEOMETRIC_PARAMETERS(:), BOUNDARY_VALUES(:)
    REAL(DP), POINTER :: TANGENTS(:,:),NORMAL(:),TIME,ANALYTIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    REAL(DP), POINTER :: independentParameters(:),dependentParameters(:)
    REAL(DP), ALLOCATABLE :: nodeData(:,:),qSpline(:),qValues(:),tValues(:),BoundaryValues(:),fittedNodes(:)
    LOGICAL :: ghostNode,nodeExists,importDataFromFile,ALENavierStokesEquationsSetFound=.FALSE.
    LOGICAL :: SolidEquationsSetFound=.FALSE.,SolidNodeFound=.FALSE.,FluidEquationsSetFound=.FALSE.,parameterSetCreated
    CHARACTER(70) :: inputFile,tempString

    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)
    NULLIFY(EQUATIONS)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(ANALYTIC_FIELD)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(materialsField)
    NULLIFY(independentField)
    NULLIFY(ANALYTIC_VARIABLE)
    NULLIFY(FIELD_VARIABLE)
    NULLIFY(GEOMETRIC_VARIABLE)
    NULLIFY(MATERIALS_VARIABLE)
    NULLIFY(DOMAIN)
    NULLIFY(DOMAIN_NODES)
    NULLIFY(INTERPOLATED_POINT)
    NULLIFY(INTERPOLATION_PARAMETERS)
    NULLIFY(MESH_VELOCITY_VALUES)
    NULLIFY(GEOMETRIC_PARAMETERS)
    NULLIFY(BOUNDARY_VALUES)
    NULLIFY(TANGENTS)
    NULLIFY(NORMAL)
    NULLIFY(TIME)
    NULLIFY(ANALYTIC_PARAMETERS)
    NULLIFY(MATERIALS_PARAMETERS)
    NULLIFY(independentParameters)
    NULLIFY(dependentParameters)

    ENTERS("NavierStokes_PreSolveUpdateBoundaryConditions",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,startTime,stopTime,currentTime,timeIncrement, &
          & currentTimeLoopIteration,outputIterationNumber,ERR,ERROR,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet

                    ! Fitting boundary condition- get values from file
                    ! TODO: this should be generalised with input filenames specified from the example file when IO is improved
                    IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
                      !Read in field values to independent field
                      NULLIFY(independentFieldVariable)
                      NULLIFY(dependentFieldVariable)
                      independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                      dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      independentVariableType=independentField%VARIABLES(1)%VARIABLE_TYPE
                      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentFieldVariable,err,error,*999)
                      dependentVariableType=dependentField%VARIABLES(1)%VARIABLE_TYPE
                      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentFieldVariable,err,error,*999)
                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_equations%BOUNDARY_CONDITIONS, &
                        & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                      !Read in field data from file
                      !Loop over nodes and update independent field values - if a fixed fitted boundary, also update dependent
                      IF(ASSOCIATED(independentField)) THEN
                        componentNumberVelocity = 1
                        numberOfDimensions = dependentFieldVariable%NUMBER_OF_COMPONENTS - 1
                        ! Get the nodes on this computational domain
                        IF(independentFieldVariable%COMPONENTS(componentNumberVelocity)%INTERPOLATION_TYPE== &
                          & FIELD_NODE_BASED_INTERPOLATION) THEN
                          domain=>independentFieldVariable%COMPONENTS(componentNumberVelocity)%DOMAIN
                          IF(ASSOCIATED(domain)) THEN
                            IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                              DOMAIN_NODES=>domain%TOPOLOGY%NODES
                              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                numberOfNodes = DOMAIN_NODES%NUMBER_OF_NODES
                                numberOfGlobalNodes = DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES
                              ELSE
                                CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Domain topology is not associated.",err,error,*999)
                            END IF
                          ELSE
                            CALL FlagError("Domain is not associated.",err,error,*999)
                          END IF
                        ELSE
                          CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                        END IF

                        ! Construct the filename based on the computational node and time step
                        inputFile = './../interpolatedData/fitData' //TRIM(NUMBER_TO_VSTRING(currentTimeLoopIteration, &
                          & "*",ERR,ERROR)) // '.dat'

                        INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                        IF(importDataFromFile) THEN
                          !Read fitted data from input file (if exists)
                          IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                            & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Updating independent field and boundary nodes from "// &
                            & inputFile,err,error,*999)
                          OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                          READ(10,*) numberOfFittedNodes
                          ALLOCATE(fittedNodes(numberOfFittedNodes))
                          READ(10,*) fittedNodes
                          DO nodeIdx=1, numberOfFittedNodes
                            userNodeNumber=INT(fittedNodes(nodeIdx),INTG)
                            CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(domain%Topology,userNodeNumber,nodeExists, &
                              & localNodeNumber,ghostNode,err,error,*999)
                            IF(nodeExists .AND. .NOT. ghostNode) THEN
                              ! Node found on this computational node
                              READ(10,*) (componentValues(componentIdx), componentIdx=1,numberOfDimensions)
                              DO componentIdx=1,numberOfDimensions
                                dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(localNodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                                independentDof = independentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                  & NODE_PARAM2DOF_MAP%NODES(localNodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                                VALUE = componentValues(componentIdx)
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(INdependentField,independentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,independentDof,VALUE,ERR,ERROR,*999)
                                CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField,dependentVariableType,1,1, &
                                  & userNodeNumber,componentIdx,localDof,globalDof,err,error,*999)
                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(globalDof)
                                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,dependentVariableType, &
                                    & FIELD_VALUES_SET_TYPE,localDof,VALUE,ERR,ERROR,*999)
                                END IF
                              END DO !componentIdx
                            ELSE
                              ! Dummy read if this node not on this computational node
                              READ(10,*)
                            END IF
                          END DO
                          DEALLOCATE(fittedNodes)
                          CLOSE(UNIT=10)
                          ! Update any distributed field values
                          CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType,FIELD_VALUES_SET_TYPE, &
                            & err,error,*999)
                          CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType,FIELD_VALUES_SET_TYPE, &
                            & err,error,*999)
                          CALL Field_ParameterSetUpdateStart(INdependentField,independentVariableType,FIELD_VALUES_SET_TYPE, &
                            & err,error,*999)
                          CALL Field_ParameterSetUpdateFinish(INdependentField,independentVariableType,FIELD_VALUES_SET_TYPE, &
                            & err,error,*999)
                        END IF !check import file exists
                      ELSE
                        CALL FlagError("Equations set independent field is not associated.",ERR,ERROR,*999)
                      END IF
                    END IF !Equations set independent

                    ! Analytic equations
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      !Standard analytic functions
                      IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                        ! Update analytic time value with current time
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                        !Calculate analytic values
                        BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,ERR,ERROR,*999)
                        END IF
                      ELSE IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                        & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN.OR. &
                        & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                        & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE) THEN
                        IF(ASSOCIATED(EQUATIONS_SET)) THEN
                          IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                            dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            IF(ASSOCIATED(dependentField)) THEN
                              geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                              IF(ASSOCIATED(geometricField)) THEN
                                ! Geometric parameters
                                CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                                  & err,error,*999)
                                NULLIFY(GEOMETRIC_VARIABLE)
                                CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
                                NULLIFY(GEOMETRIC_PARAMETERS)
                                CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                  & GEOMETRIC_PARAMETERS,err,error,*999)
                                ! Analytic parameters
                                ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                                NULLIFY(ANALYTIC_VARIABLE)
                                NULLIFY(ANALYTIC_PARAMETERS)
                                IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                  CALL Field_VariableGet(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                    & ANALYTIC_PARAMETERS,err,error,*999)
                                END IF
                                ! Materials parameters
                                NULLIFY(materialsField)
                                NULLIFY(MATERIALS_VARIABLE)
                                NULLIFY(MATERIALS_PARAMETERS)
                                IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                                  materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                  CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,err,error,*999)
                                  CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                    & MATERIALS_PARAMETERS,err,error,*999)
                                END IF
                                DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                                  variable_type=dependentField%VARIABLES(variable_idx)%VARIABLE_TYPE
                                  FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                    CALL Field_ParameterSetEnsureCreated(dependentField,variable_type, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                    DO componentIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                      IF(FIELD_VARIABLE%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                        & FIELD_NODE_BASED_INTERPOLATION) THEN
                                        DOMAIN=>FIELD_VARIABLE%COMPONENTS(componentIdx)%DOMAIN
                                        IF(ASSOCIATED(DOMAIN)) THEN
                                          IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                            DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                            IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                              !Should be replaced by boundary node flag
                                              DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                                element_idx=DOMAIN%topology%nodes%nodes(node_idx)%surrounding_elements(1)
                                                CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,element_idx, &
                                                  & INTERPOLATION_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                en_idx=0
                                                XI_COORDINATES=0.0_DP
                                                number_of_nodes_xic(1)=DOMAIN%topology%elements%elements(element_idx)% &
                                                  & basis%number_of_nodes_xic(1)
                                                number_of_nodes_xic(2)=DOMAIN%topology%elements%elements(element_idx)% &
                                                  & basis%number_of_nodes_xic(2)
                                                IF(NUMBER_OF_DIMENSIONS==3) THEN
                                                  number_of_nodes_xic(3)=DOMAIN%topology%elements%elements(element_idx)%basis% &
                                                    & number_of_nodes_xic(3)
                                                ELSE
                                                  number_of_nodes_xic(3)=1
                                                END IF
                                                !\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                                                IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4 .OR. &
                                                  & DOMAIN%topology%elements%maximum_number_of_element_parameters==9 .OR. &
                                                  & DOMAIN%topology%elements%maximum_number_of_element_parameters==16 .OR. &
                                                  & DOMAIN%topology%elements%maximum_number_of_element_parameters==8 .OR. &
                                                  & DOMAIN%topology%elements%maximum_number_of_element_parameters==27 .OR. &
                                                  & DOMAIN%topology%elements%maximum_number_of_element_parameters==64) THEN
                                                  DO K=1,number_of_nodes_xic(3)
                                                    DO J=1,number_of_nodes_xic(2)
                                                      DO I=1,number_of_nodes_xic(1)
                                                        en_idx=en_idx+1
                                                        IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                          & element_nodes(en_idx)==node_idx) EXIT
                                                        XI_COORDINATES(1)=XI_COORDINATES(1)+(1.0_DP/(number_of_nodes_xic(1)-1))
                                                      END DO !I
                                                      IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                        & element_nodes(en_idx)==node_idx) EXIT
                                                      XI_COORDINATES(1)=0.0_DP
                                                      XI_COORDINATES(2)=XI_COORDINATES(2)+(1.0_DP/(number_of_nodes_xic(2)-1))
                                                    END DO !J
                                                    IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                      & element_nodes(en_idx)==node_idx) EXIT
                                                    XI_COORDINATES(1)=0.0_DP
                                                    XI_COORDINATES(2)=0.0_DP
                                                    IF(number_of_nodes_xic(3)/=1) THEN
                                                      XI_COORDINATES(3)=XI_COORDINATES(3)+(1.0_DP/(number_of_nodes_xic(3)-1))
                                                    END IF
                                                  END DO !K
                                                  CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,XI_COORDINATES, &
                                                    & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                ELSE
                                                  !\todo: Use boundary flag
                                                  IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==3) THEN
                                                    T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                    T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                  ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==6) THEN
                                                    T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                    T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                    T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                                                    T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                                                    T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                                                  ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                                    & NUMBER_OF_DIMENSIONS==2) THEN
                                                    T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                                                    T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                                                    T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                    T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                                                    T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                  ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==4) THEN
                                                    T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                    T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                    T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                  ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==10.AND. &
                                                    & NUMBER_OF_DIMENSIONS==3) THEN
                                                    T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                    T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                    T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                                                    T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                                                    T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                                                    T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                                                    T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                                                  ELSE IF(DOMAIN%topology%elements%maximum_number_of_element_parameters==20) THEN
                                                    T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                                                    T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                                                    T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                    T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                                                    T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                                                    T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                                                    T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                                                    T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                                                    T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                                                  END IF
                                                  DO K=1,DOMAIN%topology%elements%maximum_number_of_element_parameters
                                                    IF(DOMAIN%topology%elements%elements(element_idx)% &
                                                      & element_nodes(K)==node_idx) EXIT
                                                  END DO !K
                                                  IF(NUMBER_OF_DIMENSIONS==2) THEN
                                                    CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                                                      & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                  ELSE IF(NUMBER_OF_DIMENSIONS==3) THEN
                                                    CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                                                      & INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                                  END IF
                                                END IF
                                                X=0.0_DP
                                                DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                                  X(dim_idx)=INTERPOLATED_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dim_idx,1)
                                                END DO !dim_idx
                                                !Loop over the derivatives
                                                DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                                  ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                                  GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                                    & GLOBAL_DERIVATIVE_INDEX
                                                  materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                                                  !Define MU_PARAM, density=1
                                                  MU_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                    & parameters%cmiss%dataDP(1)
                                                  !Define RHO_PARAM, density=2
                                                  RHO_PARAM=materialsField%variables(1)%parameter_sets%parameter_sets(1)%ptr% &
                                                    & parameters%cmiss%dataDP(2)
                                                  CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(ANALYTIC_FUNCTION_TYPE,X, &
                                                    & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX,componentIdx, &
                                                    & NUMBER_OF_DIMENSIONS,FIELD_VARIABLE%NUMBER_OF_COMPONENTS, &
                                                    & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,err,error,*999)
                                                  !Default to version 1 of each node derivative
                                                  local_ny=FIELD_VARIABLE%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                                  CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,dependentField% &
                                                    & VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr,BOUNDARY_CONDITIONS_VARIABLE, &
                                                    & err,error,*999)
                                                  IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                                    BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                      & CONDITION_TYPES(local_ny)
                                                    IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                        & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                                        & VALUE,err,error,*999)
                                                    END IF
                                                  ELSE
                                                    CALL FlagError("Boundary conditions U variable is not associated.", &
                                                      & err,error,*999)
                                                  END IF
                                                END DO !deriv_idx
                                              END DO !node_idx
                                            ELSE
                                              CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                            END IF
                                          ELSE
                                            CALL FlagError("Domain topology is not associated.",err,error,*999)
                                          END IF
                                        ELSE
                                          CALL FlagError("Domain is not associated.",err,error,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                      END IF
                                    END DO !componentIdx
                                    CALL Field_ParameterSetUpdateStart(dependentField,variable_type, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                    CALL Field_ParameterSetUpdateFinish(dependentField,variable_type, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                    CALL Field_ParameterSetUpdateStart(dependentField,variable_type, &
                                      & FIELD_VALUES_SET_TYPE,err,error,*999)
                                    CALL Field_ParameterSetUpdateFinish(dependentField,variable_type, &
                                      & FIELD_VALUES_SET_TYPE,err,error,*999)
                                  ELSE
                                    CALL FlagError("Field variable is not associated.",err,error,*999)
                                  END IF
                                END DO !variable_idx
                                CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,&
                                  & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
                              ELSE
                                CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                            END IF
                          ELSE
                            CALL FlagError("Equations set analytic is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Equations set is not associated.",err,error,*999)
                        END IF
                      END IF !Standard/unit analytic subtypes

                    END IF ! Analytic boundary conditions

                    !TODO implement non-analytic time-varying boundary conditions (i.e. from file)
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
               ENDIF

            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                !If analytic flow waveform, calculate and update
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                      EQUATIONS_SET=>equations%equationsSet
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                        CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                          & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                          ! Calculate analytic values
                          CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*999)
                        CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                          ! Perform spline interpolation of values from a file
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                          DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
                            dependentVariableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
                            NULLIFY(dependentFieldVariable)
                            CALL Field_VariableGet(dependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
                              & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                              IF(ASSOCIATED(dependentFieldVariable)) THEN
                                DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                  IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
                                    IF(ASSOCIATED(domain)) THEN
                                      IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                          ! Create the analytic field values type on the dependent field if it does not exist
                                          CALL FIELD_PARAMETER_SET_CREATED(dependentField,dependentVariableType, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,parameterSetCreated,ERR,ERROR,*999)
                                          IF (.NOT. parameterSetCreated) THEN
                                            CALL FIELD_PARAMETER_SET_CREATE(dependentField,dependentVariableType, &
                                              & FIELD_ANALYTIC_VALUES_SET_TYPE,ERR,ERROR,*999)
                                          END IF
                                          !Loop over the local nodes excluding the ghosts.
                                          DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                            userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                            DO derivativeIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                              DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                  & numberOfVersions
                                                dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                  & VERSIONS(versionIdx)
                                                ! Update dependent field value if this is a splint BC
                                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                  & CONDITION_TYPES(dependentDof)
                                                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                                  !Update analytic field if file exists and dependent field if boundary condition set
                                                  inputFile = './input/interpolatedData/1D/'
                                                  IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) THEN
                                                    inputFile = TRIM(inputFile) // 'U/component'
                                                  END IF
                                                  WRITE(tempString,"(I1.1)") componentIdx
                                                  inputFile = TRIM(inputFile) // tempString(1:1) // '/derivative'
                                                  WRITE(tempString,"(I1.1)") derivativeIdx
                                                  inputFile = TRIM(inputFile) // tempString(1:1) // '/version'
                                                  WRITE(tempString,"(I1.1)") versionIdx
                                                  inputFile = TRIM(inputFile) // tempString(1:1) // '/'
                                                  WRITE(tempString,"(I4.4)") userNodeNumber
                                                  inputFile = TRIM(inputFile) // tempString(1:4) // '.dat'
                                                  inputFile = TRIM(inputFile)
                                                  INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                                                  IF(importDataFromFile) THEN
                                                    !Read fitted data from input file (if exists)
                                                    OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                                                    ! Header timeData = numberOfTimesteps
                                                    READ(10,*) timeData
                                                    numberOfSourceTimesteps = INT(timeData)
                                                    ALLOCATE(nodeData(numberOfSourceTimesteps,2))
                                                    ALLOCATE(qValues(numberOfSourceTimesteps))
                                                    ALLOCATE(tValues(numberOfSourceTimesteps))
                                                    ALLOCATE(qSpline(numberOfSourceTimesteps))
                                                    nodeData = 0.0_DP
                                                    ! Read in time and dependent value
                                                    DO timeIdx=1,numberOfSourceTimesteps
                                                      READ(10,*) (nodeData(timeIdx,component_idx), component_idx=1,2)
                                                    END DO
                                                    CLOSE(UNIT=10)
                                                    tValues = nodeData(:,1)
                                                    qValues = nodeData(:,2)
                                                    CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues, &
                                                      & 2,0.0_DP,2,0.0_DP,qSpline,err,error,*999)
                                                    CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues,qSpline, &
                                                      & CURRENT_TIME,VALUE,QP,QPP,err,error,*999)
                                                    DEALLOCATE(nodeData)
                                                    DEALLOCATE(qSpline)
                                                    DEALLOCATE(qValues)
                                                    DEALLOCATE(tValues)
                                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                      & dependentVariableType,FIELD_VALUES_SET_TYPE,dependentDof, &
                                                      & VALUE,ERR,ERROR,*999)
                                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                      & dependentVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE,dependentDof, &
                                                      & VALUE,ERR,ERROR,*999)
                                                  END IF
                                                END IF ! check if import data file exists
                                              END DO !versionIdx
                                            END DO !derivativeIdx
                                          END DO !nodeIdx
                                          ! Update distributed field values
                                          CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                            & FIELD_VALUES_SET_TYPE,err,error,*999)
                                          CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                            & FIELD_VALUES_SET_TYPE,err,error,*999)
                                          CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                          CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                        ELSE
                                          CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("Domain topology is not associated.",err,error,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                  END IF
                                END DO !componentIdx
                              ELSE
                                CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                              END IF
                            END IF
                          END DO !variableIdx
                        CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART)
                          ! Using heart lumped parameter model for input
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                          materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                          CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,5, &
                            & Lref,err,error,*999)
                          CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,6, &
                            & Tref,err,error,*999)
                          CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,7, &
                            & Mref,err,error,*999)
                          DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
                            dependentVariableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
                            NULLIFY(dependentFieldVariable)
                            CALL Field_VariableGet(dependentField,dependentVariableType,dependentFieldVariable,err,error,*999)
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,dependentFieldVariable, &
                              & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                            IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                              IF(ASSOCIATED(dependentFieldVariable)) THEN
                                DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                  IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                    & FIELD_NODE_BASED_INTERPOLATION) THEN
                                    domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
                                    IF(ASSOCIATED(domain)) THEN
                                      IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                        DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                        IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                          !Loop over the local nodes excluding the ghosts.
                                          DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                            userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                            DO derivativeIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                              DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                 & numberOfVersions
                                                dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                                  & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                                                  & VERSIONS(versionIdx)
                                                BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                  & CONDITION_TYPES(dependentDof)
                                                IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE, &
                                                    & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,1,VALUE, &
                                                    & err,error,*999)
                                                  ! Convert Q from ml/s to non-dimensionalised form.
                                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,dependentVariableType, &
                                                    & FIELD_VALUES_SET_TYPE,dependentDof,((Lref**3.0)/Tref)*VALUE,err,error,*999)
                                                END IF
                                              END DO !versionIdx
                                            END DO !derivativeIdx
                                          END DO !nodeIdx
                                          ! Update distributed field values
                                          CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                            & FIELD_VALUES_SET_TYPE,err,error,*999)
                                          CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                            & FIELD_VALUES_SET_TYPE,err,error,*999)
                                        ELSE
                                          CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                        END IF
                                      ELSE
                                        CALL FlagError("Domain topology is not associated.",err,error,*999)
                                      END IF
                                    ELSE
                                      CALL FlagError("Domain is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                                  END IF
                                END DO !componentIdx
                              ELSE
                                CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                              END IF
                            END IF
                          END DO !variableIdx
                        CASE DEFAULT
                          ! Do nothing (might have another use for analytic equations)
                        END SELECT
                      END IF ! Check for analytic equations
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations are not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                END IF ! solver equations associated
                ! Update any multiscale boundary values (coupled 0D or non-reflecting)
                CALL NavierStokes_UpdateMultiscaleBoundary(EQUATIONS_SET,BOUNDARY_CONDITIONS,TIME_INCREMENT,err,error,*999)
            CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
              & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
              ! TODO: this should be set up so it uses the same pre_solve steps as the individual 3D/1D equations sets
              SELECT CASE(SOLVER%SOLVE_TYPE)
                ! --- D y n a m i c    S o l v e r s ---
              CASE(SOLVER_DYNAMIC_TYPE)
                CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,startTime,stopTime,currentTime,timeIncrement,currentTimeLoopIteration, &
                  & outputIterationNumber,ERR,ERROR,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    DO equationsSetIdx = 1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equationsSetIdx)%PTR
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                          ! --- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
                        CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                          ! Fitting boundary condition- get values from file
                          ! TODO: this should be generalised when IO is improved
                          IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
                            !Read in field values to independent field
                            NULLIFY(independentFieldVariable)
                            NULLIFY(dependentFieldVariable)
                            independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                            dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                            independentVariableType=independentField%VARIABLES(1)%VARIABLE_TYPE
                            CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE, &
                              & independentFieldVariable,ERR,ERROR,*999)
                            dependentVariableType=dependentField%VARIABLES(1)%VARIABLE_TYPE
                            CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentFieldVariable,ERR,ERROR,*999)
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_EQUATIONS%BOUNDARY_CONDITIONS, &
                              & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                            !Read in field data from file
                            !Loop over nodes and update independent field values. If a fixed fitted boundary, also update dependent
                            IF(ASSOCIATED(independentField)) THEN
                              componentNumberVelocity = 1
                              numberOfDimensions = dependentFieldVariable%NUMBER_OF_COMPONENTS - 1
                              ! Get the nodes on this computational domain
                              IF(independentFieldVariable%COMPONENTS(componentNumberVelocity)%INTERPOLATION_TYPE== &
                                & FIELD_NODE_BASED_INTERPOLATION) THEN
                                domain=>independentFieldVariable%COMPONENTS(componentNumberVelocity)%DOMAIN
                                IF(ASSOCIATED(domain)) THEN
                                  IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      numberOfNodes = DOMAIN_NODES%NUMBER_OF_NODES
                                      numberOfGlobalNodes = DOMAIN_NODES%NUMBER_OF_GLOBAL_NODES
                                    ELSE
                                      CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                  END IF
                                ELSE
                                  CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                END IF
                              ELSE
                                CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                              END IF

                              ! Construct the filename based on the computational node and time step
                              inputFile = './../interpolatedData/fitData' //TRIM(NUMBER_TO_VSTRING(currentTimeLoopIteration, &
                                & "*",ERR,ERROR)) // '.dat'

                              INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                              IF(importDataFromFile) THEN
                                !Read fitted data from input file (if exists)
                                computationalNode = ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
                                IF(computationalNode==0) THEN
                                  CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Updating independent field and boundary nodes from " &
                                    & //inputFile,ERR,ERROR,*999)
                                END IF
                                OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')

                                READ(10,*) numberOfFittedNodes
                                ALLOCATE(fittedNodes(numberOfFittedNodes))
                                READ(10,*) fittedNodes
                                DO nodeIdx=1, numberOfFittedNodes
                                  userNodeNumber=INT(fittedNodes(nodeIdx),INTG)
                                  CALL DOMAIN_TOPOLOGY_NODE_CHECK_EXISTS(domain%Topology,userNodeNumber,nodeExists, &
                                    & localNodeNumber,ghostNode,err,error,*999)
                                  IF(nodeExists .AND. .NOT. ghostNode) THEN
                                    ! Node found on this computational node
                                    READ(10,*) (componentValues(componentIdx), componentIdx=1,numberOfDimensions)
                                    DO componentIdx=1,numberOfDimensions
                                      dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(localNodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                                      independentDof = independentFieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                        & NODE_PARAM2DOF_MAP%NODES(localNodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                                      VALUE = componentValues(componentIdx)
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(independentField,independentVariableType, &
                                        & FIELD_VALUES_SET_TYPE,independentDof,VALUE,ERR,ERROR,*999)
                                      CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField,dependentVariableType,1,1, &
                                        & userNodeNumber,componentIdx,localDof,globalDof,err,error,*999)
                                      BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(globalDof)
                                      IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,dependentVariableType, &
                                          & FIELD_VALUES_SET_TYPE,localDof,VALUE,ERR,ERROR,*999)
                                      END IF
                                    END DO !componentIdx
                                  ELSE
                                    ! Dummy read if this node not on this computational node
                                    READ(10,*)
                                  END IF
                                END DO
                                DEALLOCATE(fittedNodes)
                                CLOSE(UNIT=10)
                                ! Update any distributed field values
                                CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                                CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                                CALL Field_ParameterSetUpdateStart(independentField,independentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                                CALL Field_ParameterSetUpdateFinish(independentField,independentVariableType, &
                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                              END IF !check import file exists
                            ELSE
                              CALL FlagError("Equations set independent field is not associated.",ERR,ERROR,*999)
                            END IF
                          END IF !Equations set independent
                          ! Analytic equations
                          IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                            !Standard analytic functions
                            IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                              ! Update analytic time value with current time
                              EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                              !Calculate analytic values
                              BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS, &
                                  & ERR,ERROR,*999)
                              END IF
                            ELSE
                              localError="Analytic equations type "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC% &
                                & ANALYTIC_FUNCTION_TYPE,"*",err,error))//" is not yet implemented for a 3D Navier-Stokes"// &
                                & " equations set for a multiscale problem."
                              CALL FlagError(localError,err,error,*999)
                            END IF
                          END IF
                          ! --- 1 D    N a v i e r - S t o k e s   E q u a t i o n s ---
                        CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                          SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                            !If analytic flow waveform, calculate and update
                            SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                            IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                              EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                              IF(ASSOCIATED(EQUATIONS)) THEN
                                EQUATIONS_SET=>EQUATIONS%equationsSet
                                IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                                  SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                  CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                                    & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                                    EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                                    BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                                    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                                      ! Calculate analytic values
                                      CALL NavierStokes_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS, &
                                        & ERR,ERROR,*999)
                                    ELSE
                                      CALL FlagError("Boundary conditions are not associated.",ERR,ERROR,*999)
                                    END IF
                                  CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                                    ! Perform spline interpolation of values from a file
                                    EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=CURRENT_TIME
                                    BOUNDARY_CONDITIONS=>SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
                                    dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                                    ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                                    DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
                                      dependentVariableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
                                      NULLIFY(dependentFieldVariable)
                                      CALL Field_VariableGet(dependentField,dependentVariableType,dependentFieldVariable, &
                                        & ERR,ERROR,*999)
                                      CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS, &
                                        & dependentFieldVariable,BOUNDARY_CONDITIONS_VARIABLE,ERR,ERROR,*999)
                                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                        IF(ASSOCIATED(dependentFieldVariable)) THEN
                                          DO componentIdx=1,dependentFieldVariable%NUMBER_OF_COMPONENTS
                                            IF(dependentFieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE== &
                                              & FIELD_NODE_BASED_INTERPOLATION) THEN
                                              domain=>dependentFieldVariable%COMPONENTS(componentIdx)%DOMAIN
                                              IF(ASSOCIATED(domain)) THEN
                                                IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                                                  DOMAIN_NODES=>domain%TOPOLOGY%NODES
                                                  IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                                    !Loop over the local nodes excluding the ghosts.
                                                    DO nodeIdx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                                      userNodeNumber=DOMAIN_NODES%NODES(nodeIdx)%USER_NUMBER
                                                      DO derivativeIdx=1,DOMAIN_NODES%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                                                        DO versionIdx=1,DOMAIN_NODES%NODES(nodeIdx)% &
                                                          & DERIVATIVES(derivativeIdx)%numberOfVersions
                                                          ! Update analytic field if file exists and
                                                          ! the dependent field if boundary condition set
                                                          inputFile = './input/interpolatedData/1D/'
                                                          IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) THEN
                                                            inputFile = TRIM(inputFile) // 'U/component'
                                                          END IF
                                                          WRITE(tempString,"(I1.1)") componentIdx
                                                          inputFile = TRIM(inputFile) // tempString(1:1) // '/derivative'
                                                          WRITE(tempString,"(I1.1)") derivativeIdx
                                                          inputFile = TRIM(inputFile) // tempString(1:1) // '/version'
                                                          WRITE(tempString,"(I1.1)") versionIdx
                                                          inputFile = TRIM(inputFile) // tempString(1:1) // '/'
                                                          WRITE(tempString,"(I4.4)") userNodeNumber
                                                          inputFile = TRIM(inputFile) // tempString(1:4) // '.dat'
                                                          inputFile = TRIM(inputFile)
                                                          INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                                                          IF(importDataFromFile) THEN
                                                            ! Create the analytic field values type on the dependent field
                                                            ! if it does not exist
                                                            CALL Field_ParameterSetEnsureCreated(dependentField, &
                                                              & dependentVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                                                              & err,error,*999)
                                                            !Read fitted data from input file (if exists)
                                                            OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                                                            ! Header timeData = numberOfTimesteps
                                                            READ(10,*) timeData
                                                            numberOfSourceTimesteps = INT(timeData)
                                                            ALLOCATE(nodeData(numberOfSourceTimesteps,2))
                                                            ALLOCATE(qValues(numberOfSourceTimesteps))
                                                            ALLOCATE(tValues(numberOfSourceTimesteps))
                                                            ALLOCATE(qSpline(numberOfSourceTimesteps))
                                                            nodeData = 0.0_DP
                                                            ! Read in time and dependent value
                                                            DO timeIdx=1,numberOfSourceTimesteps
                                                              READ(10,*) (nodeData(timeIdx,component_idx), component_idx=1,2)
                                                            END DO
                                                            CLOSE(UNIT=10)
                                                            tValues = nodeData(:,1)
                                                            qValues = nodeData(:,2)
                                                            CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues, &
                                                              & 2,0.0_DP,2,0.0_DP,qSpline,err,error,*999)
                                                            CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues, &
                                                              & qSpline,CURRENT_TIME,VALUE,QP,QPP,err,error,*999)

                                                            DEALLOCATE(nodeData)
                                                            DEALLOCATE(qSpline)
                                                            DEALLOCATE(qValues)
                                                            DEALLOCATE(tValues)

                                                            dependentDof = dependentFieldVariable%COMPONENTS(componentIdx)% &
                                                              & PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(nodeIdx)% &
                                                              & DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                              & dependentVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                                                              & dependentDof,VALUE,ERR,ERROR,*999)
                                                            ! Update dependent field value if this is a splint BC
                                                            BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                                              & CONDITION_TYPES(dependentDof)
                                                            IF(BOUNDARY_CONDITION_CHECK_VARIABLE== &
                                                              & BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                                & dependentVariableType,FIELD_VALUES_SET_TYPE,dependentDof, &
                                                                & VALUE,ERR,ERROR,*999)
                                                            END IF
                                                          END IF ! check if import data file exists
                                                        END DO !versionIdx
                                                      END DO !derivativeIdx
                                                    END DO !nodeIdx
                                                    ! Update distributed field values
                                                    CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                                      & FIELD_VALUES_SET_TYPE,err,error,*999)
                                                    CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                                      & FIELD_VALUES_SET_TYPE,err,error,*999)
                                                    CALL Field_ParameterSetUpdateStart(dependentField,dependentVariableType, &
                                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                                    CALL Field_ParameterSetUpdateFinish(dependentField,dependentVariableType, &
                                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                                                  ELSE
                                                    CALL FlagError("Domain topology nodes is not associated.",ERR,ERROR,*999)
                                                  END IF
                                                ELSE
                                                  CALL FlagError("Domain topology is not associated.",ERR,ERROR,*999)
                                                END IF
                                              ELSE
                                                CALL FlagError("Domain is not associated.",ERR,ERROR,*999)
                                              END IF
                                            ELSE
                                              CALL FlagError("Only node based interpolation is implemented.",ERR,ERROR,*999)
                                            END IF
                                          END DO !componentIdx
                                        ELSE
                                          CALL FlagError("Dependent field variable is not associated.",ERR,ERROR,*999)
                                        END IF
                                      END IF
                                    END DO !variableIdx
                                  CASE DEFAULT
                                    ! Do nothing (might have another use for analytic equations)
                                  END SELECT
                                END IF ! Check for analytic equations
                              ELSE
                                CALL FlagError("Equations are not associated.",ERR,ERROR,*999)
                              END IF
                            ELSE
                              CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
                            END IF
                          END IF ! solver equations associated
                          ! Update any multiscale boundary values (coupled 0D or non-reflecting)
                          CALL NavierStokes_UpdateMultiscaleBoundary(EQUATIONS_SET,BOUNDARY_CONDITIONS, &
                            & TIME_INCREMENT,err,error,*999)
                        CASE DEFAULT
                          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*", &
                            & err,error))//" is not valid for a multiscale dynamic Navier-Stokes solver."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    END DO
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="The solve type of "//TRIM(NUMBER_TO_VSTRING(SOLVER%SOLVE_TYPE,"*",err,error))// &
                  & " is invalid for pre_solve_update_BC step of a multiscale Navier-Stokes problem type."
                CALL FlagError(localError,err,error,*999)
              END SELECT

            CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>equations%equationsSet
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                      FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                      IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                          & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                      ELSE
                        CALL FlagError("Field U variable is not associated",err,error,*999)
                      END IF
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                        NULLIFY(BOUNDARY_VALUES)
                        CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                        CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_NONLINEAR_TYPE,BOUNDARY_VALUES, &
                          & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                          & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                        DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                          variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                          FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                            DO componentIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                              DOMAIN=>FIELD_VARIABLE%COMPONENTS(componentIdx)%DOMAIN
                              IF(ASSOCIATED(DOMAIN)) THEN
                                IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                  DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                  IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                    !Loop over the local nodes excluding the ghosts.
                                    DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                      DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                        !Default to version 1 of each node derivative
                                        local_ny=FIELD_VARIABLE%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                          & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                        BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                          & CONDITION_TYPES(local_ny)
                                        IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                            & BOUNDARY_VALUES(local_ny),err,error,*999)
                                        END IF
                                      END DO !deriv_idx
                                    END DO !node_idx
                                  END IF
                                END IF
                              END IF
                            END DO !componentIdx
                          END IF
                        END DO !variable_idx

                        !\todo: This part should be read in out of a file eventually
                      ELSE
                        CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations are not associated.",err,error,*999)
              END IF
              CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,err,error,*999)
            CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the dynamic solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",err,error,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO componentIdx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(componentIdx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                          ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !componentIdx
                            END IF
                          END DO !variable_idx
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              END IF
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              !Pre solve for the linear solver
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",err,error,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_MOVED_WALL,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                          CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          !\todo: This part should be read in out of a file eventually
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                !Pre solve for the dynamic solver
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                  & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                  IF(ASSOCIATED(EQUATIONS)) THEN
                    EQUATIONS_SET=>equations%equationsSet
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                        FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                        IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                            & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                        ELSE
                          CALL FlagError("Field U variable is not associated",err,error,*999)
                        END IF
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(MESH_VELOCITY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          NULLIFY(BOUNDARY_VALUES)
                          CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                            & NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                            & CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER,CURRENT_TIME,1.0_DP,err,error,*999)
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          DISPLACEMENT_VALUE=0.0_DP
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & MESH_VELOCITY_VALUES(local_ny),err,error,*999)
                                          ELSE IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & BOUNDARY_VALUES(local_ny),err,error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                          CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                          CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                            & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                        ELSE
                          CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations are not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
                CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              END IF
              ! do nothing ???
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
                NULLIFY(Solver2)
                !Pre solve for the linear solver
                IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                    & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                    EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                    IF(ASSOCIATED(EQUATIONS)) THEN
                      EQUATIONS_SET=>equations%equationsSet
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                          IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                            CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                              & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                          ELSE
                            CALL FlagError("Field U variable is not associated",err,error,*999)
                          END IF
                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                            CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_DIMENSIONS,err,error,*999)
                            !Update moving wall nodes from solid/fluid gap (as we solve for displacements of the mesh
                            !in Laplacian smoothing step).
                            IF(CONTROL_LOOP%problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                              & CONTROL_LOOP%problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                              CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,Solver2,err,error,*999)
                            ELSE
                              CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,4,Solver2,err,error,*999)
                            ENDIF
                            IF(.NOT.ASSOCIATED(Solver2)) CALL FlagError("Dynamic solver is not associated.",Err,Error,*999)
                            !Find the FiniteElasticity equations set as there is a NavierStokes equations set too
                            SOLID_SOLVER_EQUATIONS=>Solver2%SOLVER_EQUATIONS
                            IF(ASSOCIATED(SOLID_SOLVER_EQUATIONS)) THEN
                              SOLID_SOLVER_MAPPING=>SOLID_SOLVER_equations%SOLVER_MAPPING
                              IF(ASSOCIATED(SOLID_SOLVER_MAPPING)) THEN
                                EquationsSetIndex=1
                                SolidEquationsSetFound=.FALSE.
                                DO WHILE (.NOT.SolidEquationsSetFound &
                                  & .AND.EquationsSetIndex<=SOLID_SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                                  SOLID_EQUATIONS=>SOLID_SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                                  IF(ASSOCIATED(SOLID_EQUATIONS)) THEN
                                    SOLID_EQUATIONS_SET=>SOLID_equations%equationsSet
                                    IF(ASSOCIATED(SOLID_EQUATIONS_SET)) THEN
                                      IF(SOLID_EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_ELASTICITY_CLASS &
                                        & .AND.SOLID_EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE &
                                        & .AND.((SOLID_EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE) &
                                        & .OR.(SOLID_EQUATIONS_SET%SPECIFICATION(3)== &
                                        & EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
                                        & .OR.(SOLID_EQUATIONS_SET%SPECIFICATION(3)== &
                                        & EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE))) THEN
                                        SolidEquationsSetFound=.TRUE.
                                      ELSE
                                        EquationsSetIndex=EquationsSetIndex+1
                                      END IF
                                    ELSE
                                      CALL FlagError("Solid equations set is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Solid equations not associated.",Err,Error,*999)
                                  END IF
                                END DO
                                IF(SolidEquationsSetFound.EQV..FALSE.) THEN
                                  localError="Solid equations set not found when trying to update boundary conditions."
                                  CALL FlagError(localError,Err,Error,*999)
                                END IF
                              ELSE
                                CALL FlagError("Solid solver mapping is not associated.",Err,Error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Solver equations for solid equations set not associated.",Err,Error,*999)
                            END IF
                            SOLID_DEPENDENT=>SOLID_EQUATIONS_SET%DEPENDENT
                            IF(.NOT.ASSOCIATED(SOLID_DEPENDENT%DEPENDENT_FIELD)) THEN
                              CALL FlagError("Solid equations set dependent field is not associated.",Err,Error,*999)
                            END IF
                            SOLID_dependentField=>SOLID_DEPENDENT%DEPENDENT_FIELD
                            !Find the NavierStokes equations set as there is a FiniteElasticity equations set too
                            FLUID_SOLVER_EQUATIONS=>Solver2%SOLVER_EQUATIONS
                            IF(ASSOCIATED(FLUID_SOLVER_EQUATIONS)) THEN
                              FLUID_SOLVER_MAPPING=>FLUID_SOLVER_equations%SOLVER_MAPPING
                              IF(ASSOCIATED(FLUID_SOLVER_MAPPING)) THEN
                                EquationsSetIndex=1
                                FluidEquationsSetFound=.FALSE.
                                DO WHILE (.NOT.FluidEquationsSetFound &
                                  & .AND.EquationsSetIndex<=FLUID_SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                                  FLUID_EQUATIONS=>FLUID_SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                                  IF(ASSOCIATED(SOLID_EQUATIONS)) THEN
                                    FLUID_EQUATIONS_SET=>FLUID_equations%equationsSet
                                    IF(ASSOCIATED(FLUID_EQUATIONS_SET)) THEN
                                      IF(FLUID_EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                                        & .AND.FLUID_EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                                        & .AND.(FLUID_EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
                                        & .OR.FLUID_EQUATIONS_SET%SPECIFICATION(3)== &
                                        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
                                        FluidEquationsSetFound=.TRUE.
                                      ELSE
                                        EquationsSetIndex=EquationsSetIndex+1
                                      END IF
                                    ELSE
                                      CALL FlagError("Fluid equations set is not associated.",err,error,*999)
                                    END IF
                                  ELSE
                                    CALL FlagError("Fluid equations not associated.",Err,Error,*999)
                                  END IF
                                END DO
                                IF(FluidEquationsSetFound.EQV..FALSE.) THEN
                                  localError="Fluid equations set not found when trying to update boundary conditions."
                                  CALL FlagError(localError,Err,Error,*999)
                                END IF
                              ELSE
                                CALL FlagError("Fluid solver mapping is not associated.",Err,Error,*999)
                              END IF
                            ELSE
                              CALL FlagError("Fluid equations for fluid equations set not associated.",Err,Error,*999)
                            END IF
                            FLUID_GEOMETRIC=>FLUID_EQUATIONS_SET%GEOMETRY
                            IF(.NOT.ASSOCIATED(FLUID_GEOMETRIC%GEOMETRIC_FIELD)) THEN
                              CALL FlagError("Fluid equations set geometric field is not associated",Err,Error,*999)
                            END IF
                            FLUID_GEOMETRIC_FIELD=>FLUID_GEOMETRIC%GEOMETRIC_FIELD
                            !DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_idx=1
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                            & CONDITION_TYPES(local_ny)
                                          !Update moved wall nodes only
                                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_MOVED_WALL) THEN
                                            !NOTE: assuming same mesh and mesh nodes for fluid domain and moving mesh domain
                                            FluidNodeNumber=node_idx
                                            DO search_idx=1,SIZE(Solver2%SOLVER_equations%SOLVER_MAPPING% &
                                              & INTERFACE_CONDITIONS(1)%ptr%INTERFACE% &
                                              & NODES%COUPLED_NODES(2,:))
                                              IF(Solver2%SOLVER_equations%SOLVER_MAPPING% &
                                                & INTERFACE_CONDITIONS(1)%ptr%INTERFACE% &
                                                & NODES%COUPLED_NODES(2,search_idx)==node_idx) THEN
                                                SolidNodeNumber=Solver2%SOLVER_equations%SOLVER_MAPPING% &
                                                  & INTERFACE_CONDITIONS(1)%ptr%INTERFACE% &
                                                  & NODES%COUPLED_NODES(1,search_idx)!might wanna put a break here
                                                SolidNodeFound=.TRUE.
                                              END IF
                                            END DO
                                            IF(.NOT.SolidNodeFound &
                                              & .OR.FluidNodeNumber==0) CALL FlagError("Solid interface node not found.", &
                                              & Err,Error,*999)
                                            !Default to version number 1
                                            IF(variable_idx==1) THEN
                                              CALL FIELD_PARAMETER_SET_GET_NODE(FLUID_GEOMETRIC_FIELD,variable_type, &
                                                & FIELD_VALUES_SET_TYPE,1,deriv_idx, &
                                                & FluidNodeNumber,component_idx,FluidGFValue,Err,Error,*999)
                                            ELSE
                                              FluidGFValue=0.0_DP
                                            END IF
                                            CALL FIELD_PARAMETER_SET_GET_NODE(SOLID_dependentField,variable_type, &
                                              & FIELD_VALUES_SET_TYPE,1,deriv_idx, &
                                              & SolidNodeNumber,component_idx,SolidDFValue,Err,Error,*999)
                                            NewLaplaceBoundaryValue=SolidDFValue-FluidGFValue
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & NewLaplaceBoundaryValue,Err,Error,*999)
                                          END IF
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                            !END DO !variable_idx
                          ELSE
                            CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                          END IF
                        ELSE
                          CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations are not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations are not associated.",err,error,*999)
                  END IF
                  CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,err,error,*999)
                  !Pre solve for the dynamic solver
                ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  !   CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Velocity field change boundary conditions... ",err,error,*999)
                  !   SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  !   IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  !     SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  !     !Find the NavierStokes equations set as there is a finite elasticity equations set too
                  !     EquationsSetIndex=1
                  !     ALENavierStokesEquationsSetFound=.FALSE.
                  !     DO WHILE (.NOT.ALENavierStokesEquationsSetFound &
                  !       & .AND.EquationsSetIndex<=SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS)
                  !       EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(EquationsSetIndex)%EQUATIONS
                  !       IF(ASSOCIATED(EQUATIONS)) THEN
                  !         EQUATIONS_SET=>equations%equationsSet
                  !         IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  !           IF(EQUATIONS_SET%SPECIFICATION(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                  !             & .AND.EQUATIONS_SET%SPECIFICATION(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                  !             & .AND.EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
                  !             & .AND.EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  !             ALENavierStokesEquationsSetFound=.TRUE.
                  !           ELSE
                  !             EquationsSetIndex=EquationsSetIndex+1
                  !           END IF
                  !         ELSE
                  !           CALL FlagError("ALE Navier-Stokes equations set is not associated.",err,error,*999)
                  !         END IF
                  !       ELSE
                  !         CALL FlagError("ALE equations not associated.",Err,Error,*999)
                  !       END IF
                  !     END DO
                  !     IF(ALENavierStokesEquationsSetFound.EQV..FALSE.) THEN
                  !       localError="ALE NavierStokes equations set not found when trying to update boundary conditions."
                  !       CALL FlagError(localError,Err,Error,*999)
                  !     END IF
                  !     !Get boundary conditions
                  !     BOUNDARY_CONDITIONS=>SOLVER_equations%BOUNDARY_CONDITIONS
                  !     IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                  !       FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
                  !       IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  !         CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,FIELD_VARIABLE, &
                  !           & BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                  !       ELSE
                  !         CALL FlagError("Field U variable is not associated",err,error,*999)
                  !       END IF
                  !       IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                  !         CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  !           & NUMBER_OF_DIMENSIONS,err,error,*999)
                  !         NULLIFY(MESH_VELOCITY_VALUES)
                  !         CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  !           & FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                  !         NULLIFY(BOUNDARY_VALUES)
                  !         CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  !           & FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  !         !Get update for time-dependent boundary conditions
                  !         IF(CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER==1) THEN
                  !           componentBC=1
                  !           CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                  !             & SOLVER%SOLVE_TYPE,InletNodes, &
                  !             & BoundaryValues,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                  !             & CURRENT_TIME,CONTROL_LOOP%TIME_LOOP%STOP_TIME,err,error,*999)
                  !           DO node_idx=1,SIZE(InletNodes)
                  !             CALL FIELD_PARAMETER_SET_UPDATE_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  !               & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,InletNodes(node_idx),componentBC, &
                  !               & BoundaryValues(node_idx),err,error,*999)
                  !           END DO
                  !         ELSE
                  !           !Figure out which component we're applying BC at
                  !           IF(CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER==2) THEN
                  !             componentBC=1
                  !           ELSE
                  !             componentBC=2
                  !           END IF
                  !           !Get inlet nodes and the corresponding velocities
                  !           CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                  !             & SOLVER%SOLVE_TYPE,InletNodes, &
                  !             & BoundaryValues,BOUNDARY_CONDITION_FIXED_INLET,CONTROL_LOOP%TIME_LOOP%INPUT_NUMBER, &
                  !             & CURRENT_TIME,CONTROL_LOOP%TIME_LOOP%STOP_TIME,err,error,*999)
                  !           DO node_idx=1,SIZE(InletNodes)
                  !             CALL FIELD_PARAMETER_SET_UPDATE_NODE(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  !               & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,InletNodes(node_idx),componentBC, &
                  !               & BoundaryValues(node_idx),err,error,*999)
                  !           END DO
                  !         END IF
                  !         CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  !           & FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,MESH_VELOCITY_VALUES,err,error,*999)
                  !         CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                  !           & FIELD_U_VARIABLE_TYPE,FIELD_BOUNDARY_SET_TYPE,BOUNDARY_VALUES,err,error,*999)
                  !       ELSE
                  !         CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                  !       END IF
                  !     ELSE
                  !       CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                  !     END IF
                  !   ELSE
                  !     CALL FlagError("Solver equations are not associated.",err,error,*999)
                  !   END IF
                  !   CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  !     & FIELD_VALUES_SET_TYPE,err,error,*999)
                  !   CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  !     & FIELD_VALUES_SET_TYPE,err,error,*999)
                END IF
                ! do nothing ???
              CASE DEFAULT
                localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is not valid for a FiniteElasticity-NavierStokes problem type of a multi physics problem class."
                CALL FlagError(localError,Err,Error,*999)
              END SELECT
            CASE DEFAULT
              localError="Problem type "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",err,error))// &
                & " is not valid for NAVIER_STOKES_PRE_SOLVE of a multi physics problem class."
              CALL FlagError(localError,Err,Error,*999)
            END SELECT
          CASE DEFAULT
            localError="The first problem specification of "// &
              & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(1),"*",err,error))// &
              & " is not valid for NavierStokes_PreSolveUpdateBoundaryConditions."
            CALL FlagError(localError,Err,Error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("NavierStokes_PreSolveUpdateBoundaryConditions",err,error)
    EXITS("NavierStokes_PreSolveUpdateBoundaryConditions")
    RETURN 1

  END SUBROUTINE NavierStokes_PreSolveUpdateBoundaryConditions

  !
  !================================================================================================================================
  !

  !>Update mesh velocity and move mesh for ALE Navier-Stokes problem
  SUBROUTINE NavierStokes_PreSolveAleUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EQUATIONS_SET_TYPE), POINTER :: laplaceEquationsSet,fluidEquationsSet,solidEquationsSet
    TYPE(FIELD_TYPE), POINTER :: laplaceDependentField,laplaceGeometricField,fluidIndependentField,fluidGeometricField, &
      & solidDependentField,interfaceGeometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fluidGeometricVariable,interfaceGeometricVariable
    TYPE(INTERFACE_TYPE), POINTER :: fsiInterface
    TYPE(INTERFACE_CONDITION_TYPE), POINTER :: fsiInterfaceCondition
    TYPE(NODES_TYPE), POINTER :: interfaceNodes
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: laplaceSolverEquations,fluidSolverEquations,fsiSolverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: laplaceSolverMapping,fsiSolverMapping,fluidSolverMapping,solidSolverMapping
    TYPE(SOLVER_TYPE), POINTER :: dynamicSolver,laplaceSolver
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: componentIdx,currentIteration,derivativeIdx,equationsSetIdx,fluidNode,fluidNumberOfDimensions, &
      & geometricMeshComponent,interfaceNumberOfDimensions,inputType,inputOption,laplaceNumberOfDimensions,localDOF, &
      & nodeIdx,outputIteration,solidNode,variableIdx,variableType,versionIdx
    REAL(DP) :: alpha,currentTime,previousSolidNodePosition,solidDelta,solidNodePosition,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: meshDisplacementValues(:)
    LOGICAL :: fluidEquationsSetFound=.FALSE.
    LOGICAL :: solidEquationsSetFound=.FALSE.

    ENTERS("NavierStokes_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    NULLIFY(dynamicSolver)
    NULLIFY(fluidEquationsSet)
    NULLIFY(fluidGeometricField)
    NULLIFY(fluidGeometricVariable)
    NULLIFY(fluidIndependentField)
    NULLIFY(fluidSolverEquations)
    NULLIFY(fluidSolverMapping)
    NULLIFY(fsiInterfaceCondition)
    NULLIFY(fsiInterface)
    NULLIFY(fsiSolverEquations)
    NULLIFY(fsiSolverMapping)
    NULLIFY(interfaceGeometricField)
    NULLIFY(interfaceGeometricVariable)
    NULLIFY(interfaceNodes)
    NULLIFY(laplaceEquationsSet)
    NULLIFY(laplaceDependentField)
    NULLIFY(laplaceGeometricField)
    NULLIFY(laplaceSolver)
    NULLIFY(laplaceSolverEquations)
    NULLIFY(laplaceSolverMapping)
    NULLIFY(meshDisplacementValues)
    NULLIFY(problem)
    NULLIFY(solidDependentField)
    NULLIFY(solidEquationsSet)
    NULLIFY(solidSolverMapping)
    NULLIFY(solvers)

    IF(.NOT.ASSOCIATED(SOLVER)) CALL FlagError("Solver is not associated.",err,error,*999)

    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification array is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have >= 3 entries for a Navier-Stokes problem.",err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime, &
      & currentIteration,outputIteration,err,error,*999)

    SELECT CASE(problem%specification(1))
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      SELECT CASE(problem%specification(3))
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_PGM_NAVIER_STOKES_SUBTYPE)
        !Update mesh within the dynamic solver
        IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          !Get the independent field for the ALE Navier-Stokes problem
          CALL Solvers_SolverGet(solvers,1,dynamicSolver,err,error,*999)
          CALL Solver_SolverEquationsGet(dynamicSolver,fluidSolverEquations,err,error,*999)
          CALL SolverEquations_SolverMappingGet(fluidSolverEquations,fluidSolverMapping,err,error,*999)
          CALL SolverMapping_EquationsSetGet(fluidSolverMapping,1,fluidEquationsSet,err,error,*999)
          CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
          CALL EquationsSet_IndependentFieldGet(fluidEquationsSet,fluidIndependentField,err,error,*999)
          !Get the data
          CALL Field_NumberOfComponentsGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidNumberOfDimensions,err,error,*999)
          !\todo: Introduce user calls instead of hard-coding 42/1
          !Copy input to Navier-Stokes' independent field
          inputType=42
          inputOption=1
          CALL Field_ParameterSetDataGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & meshDisplacementValues,err,error,*999)
          CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,meshDisplacementValues,fluidNumberOfDimensions,inputType, &
            & inputOption,currentIteration,1.0_DP,err,error,*999)
          CALL Field_ParameterSetUpdateStart(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateFinish(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & err,error,*999)
          !Use calculated values to update mesh
          CALL Field_ComponentMeshComponentGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          !CALL Field_ParameterSetDataGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          !  & meshDisplacementValues,err,error,*999)
          DO variableIdx=1,fluidGeometricField%NUMBER_OF_VARIABLES
            variableType=fluidGeometricField%VARIABLES(variableIdx)%VARIABLE_TYPE
            NULLIFY(fluidGeometricVariable)
            CALL Field_VariableGet(fluidGeometricField,variableType,fluidGeometricVariable,err,error,*999)
            DO componentIdx=1,fluidGeometricVariable%NUMBER_OF_COMPONENTS
              NULLIFY(domain)
              CALL FieldVariable_DomainGet(fluidGeometricVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                  DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                    localDOF=fluidGeometricVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                      & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(derivativeIdx)%versions(versionIdx)
                    CALL Field_ParameterSetAddLocalDOF(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOF, &
                      & meshDisplacementValues(localDOF),err,error,*999)
                  ENDDO !versionIdx
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
          ENDDO !variableIdx
          CALL Field_ParameterSetDataRestore(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & meshDisplacementValues,err,error,*999)
          CALL Field_ParameterSetUpdateStart(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Now use displacement values to calculate velocity values
          alpha=1.0_DP/timeIncrement
          CALL Field_ParameterSetsCopy(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
        ELSE
          CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
        END IF
      CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        !Update mesh within the dynamic solver
        IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
          IF(solver%DYNAMIC_SOLVER%ALE) THEN
            !Get the dependent field for the three component Laplace problem
            CALL Solvers_SolverGet(solvers,1,laplaceSolver,err,error,*999)
            CALL Solver_SolverEquationsGet(laplaceSolver,laplaceSolverEquations,err,error,*999)
            CALL SolverEquations_SolverMappingGet(laplaceSolverEquations,laplaceSolverMapping,err,error,*999)
            NULLIFY(laplaceEquationsSet)
            CALL SolverMapping_EquationsSetGet(laplaceSolverMapping,1,laplaceEquationsSet,err,error,*999)
            CALL EquationsSet_DependentFieldGet(laplaceEquationsSet,laplaceDependentField,err,error,*999)
            CALL EquationsSet_GeometricFieldGet(laplaceEquationsSet,laplaceGeometricField,err,error,*999)
            CALL Field_NumberOfComponentsGet(laplaceGeometricField,FIELD_U_VARIABLE_TYPE,laplaceNumberOfDimensions, &
              & err,error,*999)
            !Get the independent field for the ALE Navier-Stokes problem
            CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
            CALL Solver_SolverEquationsGet(dynamicSolver,fluidSolverEquations,err,error,*999)
            CALL SolverEquations_SolverMappingGet(fluidSolverEquations,fluidSolverMapping,err,error,*999)
            CALL SolverMapping_EquationsSetGet(fluidSolverMapping,1,fluidEquationsSet,err,error,*999)
            CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
            CALL EquationsSet_IndependentFieldGet(fluidEquationsSet,fluidIndependentField,err,error,*999)
            CALL Field_NumberOfComponentsGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidNumberOfDimensions, &
              & err,error,*999)
            !Copy result from Laplace mesh movement to Navier-Stokes' independent field
            IF(fluidNumberOfDimensions/=laplaceNumberOfDimensions) &
              & CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",err,error,*999)
            DO componentIdx=1,fluidNumberOfDimensions
              CALL Field_ParametersToFieldParametersCopy(laplaceDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,componentIdx, &
                & err,error,*999)
            ENDDO !componentIdx
            !Use calculated values to update mesh
            CALL Field_ComponentMeshComponentGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ParameterSetDataGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
              & meshDisplacementValues,err,error,*999)
            NULLIFY(fluidGeometricVariable)
            CALL Field_VariableGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidGeometricVariable,err,error,*999)
            DO componentIdx=1,fluidGeometricVariable%NUMBER_OF_COMPONENTS
              NULLIFY(domain)
              CALL FieldVariable_DomainGet(fluidGeometricVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                  !Default to version 1 of each node derivative
                  localDOF=fluidGeometricVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                    & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(derivativeIdx)%versions(1)
                  CALL Field_ParameterSetAddLocalDOF(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOF, &
                      & meshDisplacementValues(localDOF),err,error,*999)
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
            CALL Field_ParameterSetDataRestore(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
              & meshDisplacementValues,err,error,*999)
            CALL Field_ParameterSetUpdateStart(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            !Now use displacement values to calculate velocity values
            alpha=1.0_DP/timeIncrement
            CALL Field_ParameterSetsCopy(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
              & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
          ELSE
            CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
        END IF
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      SELECT CASE(problem%specification(2))
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        SELECT CASE(problem%specification(3))
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
          !Update mesh within the dynamic solver
          IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
            IF(solver%DYNAMIC_SOLVER%ALE) THEN
              !Get the dependent field for the Laplace problem
              IF(problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                & problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                CALL Solvers_SolverGet(solvers,3,laplaceSolver,err,error,*999)
              ELSE IF(problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                & problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                CALL Solvers_SolverGet(solvers,5,laplaceSolver,err,error,*999)
              ELSE
                CALL Solvers_SolverGet(solvers,4,laplaceSolver,err,error,*999)
              ENDIF
              CALL Solver_SolverEquationsGet(laplaceSolver,laplaceSolverEquations,err,error,*999)
              CALL SolverEquations_SolverMappingGet(laplaceSolverEquations,laplaceSolverMapping,err,error,*999)
              CALL SolverMapping_EquationsSetGet(laplaceSolverMapping,1,laplaceEquationsSet,err,error,*999)
              CALL EquationsSet_GeometricFieldGet(laplaceEquationsSet,laplaceGeometricField,err,error,*999)
              CALL EquationsSet_DependentFieldGet(laplaceEquationsSet,laplaceDependentField,err,error,*999)
              CALL Field_NumberOfComponentsGet(laplaceGeometricField,FIELD_U_VARIABLE_TYPE,laplaceNumberOfDimensions, &
                & err,error,*999)
              !Get the independent field for the ALE Navier-Stokes problem
              !Get the dynamic solver
              IF(problem%specification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                & problem%specification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
              ELSE IF(problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
                & problem%specification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
                CALL Solvers_SolverGet(solvers,4,dynamicSolver,err,error,*999)
              ELSE
                CALL Solvers_SolverGet(solvers,3,dynamicSolver,err,error,*999)
              ENDIF
              CALL Solver_SolverEquationsGet(dynamicSolver,fsiSolverEquations,err,error,*999)
              CALL SolverEquations_SolverMappingGet(fsiSolverEquations,fsiSolverMapping,err,error,*999)
              !Find the Navier Stokes equations set
              equationsSetIdx=1
              fluidEquationsSetFound=.FALSE.
              DO WHILE (.NOT.fluidEquationsSetFound.AND.equationsSetIdx<=fsiSolverMapping%NUMBER_OF_EQUATIONS_SETS)
                NULLIFY(fluidEquationsSet)
                CALL SolverMapping_EquationsSetGet(fsiSolverMapping,equationsSetIdx,fluidEquationsSet,err,error,*999)
                IF(fluidEquationsSet%specification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                  & .AND.fluidEquationsSet%specification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                  & .AND.(fluidEquationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
                  & .OR.fluidEquationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
                  fluidEquationsSetFound=.TRUE.
                ELSE
                  equationsSetIdx=equationsSetIdx+1
                ENDIF
              ENDDO
              IF(.NOT.fluidEquationsSetFound) &
                & CALL FlagError("ALE NavierStokes equations set not found when trying to update ALE mesh.",err,error,*999)
              CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
              CALL EquationsSet_IndependentFieldGet(fluidEquationsSet,fluidIndependentField,err,error,*999)
              CALL Field_NumberOfComponentsGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidNumberOfDimensions, &
                & err,error,*999)
              !Copy result from Laplace mesh movement to Navier-Stokes' independent field
              IF(fluidNumberOfDimensions/=laplaceNumberOfDimensions) &
                & CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",err,error,*999)
              DO componentIdx=1,fluidNumberOfDimensions
                CALL Field_ParametersToFieldParametersCopy(laplaceDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & componentIdx,fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,componentIdx, &
                  & err,error,*999)
              ENDDO !componentIdx
              !Find the solid mechanics equations set
              equationsSetIdx=1
              solidEquationsSetFound=.FALSE.
              DO WHILE (.NOT.solidEquationsSetFound.AND.equationsSetIdx<=fsiSolverMapping%NUMBER_OF_EQUATIONS_SETS)
                NULLIFY(solidEquationsSet)
                CALL SolverMapping_EquationsSetGet(fsiSolverMapping,equationsSetIdx,solidEquationsSet,err,error,*999)
                IF(solidEquationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
                  & solidEquationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
                  & (solidEquationsSet%specification(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE.OR. &
                  & solidEquationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE.OR. &
                  & solidEquationsSet%specification(3)==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE)) THEN
                  solidEquationsSetFound=.TRUE.
                ELSE
                  equationsSetIdx=equationsSetIdx+1
                ENDIF
              ENDDO
              IF(.NOT.fluidEquationsSetFound) &
                & CALL FlagError("Solid equations set not found when trying to update ALE mesh.",err,error,*999)
              CALL EquationsSet_DependentFieldGet(solidEquationsSet,solidDependentField,err,error,*999)
              !Loop over the interface nodes and see if there has been any change in the solid position since the
              !beginning of the time loop. Add this change to the mesh displacements.
              CALL SolverMapping_InterfaceConditionGet(fsiSolverMapping,1,fsiInterfaceCondition,err,error,*999)
              CALL InterfaceCondition_InterfaceGet(fsiInterfaceCondition,fsiInterface,err,error,*999)
              CALL Interface_NodesGet(fsiInterface,interfaceNodes,err,error,*999)
              CALL InterfaceCondition_GeometricFieldGet(fsiInterfaceCondition,interfaceGeometricField,err,error,*999)
              CALL Field_NumberOfComponentsGet(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,interfaceNumberOfDimensions, &
                & err,error,*999)
              CALL Field_VariableGet(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,interfaceGeometricVariable,err,error,*999)
              DO componentIdx=1,interfaceGeometricVariable%NUMBER_OF_COMPONENTS
                SELECT CASE(interfaceGeometricVariable%components(componentIdx)%INTERPOLATION_TYPE)
                CASE(FIELD_NODE_BASED_INTERPOLATION)
                  NULLIFY(domain)
                  CALL FieldVariable_DomainGet(interfaceGeometricVariable,componentIdx,domain,err,error,*999)
                  NULLIFY(domainTopology)
                  CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
                  NULLIFY(domainNodes)
                  CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
                  DO nodeIdx=1,domainNodes%TOTAL_NUMBER_OF_NODES
                    solidNode=interfaceNodes%COUPLED_NODES(1,nodeIdx)
                    fluidNode=interfaceNodes%COUPLED_NODES(2,nodeIdx)
                    DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                      DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                        CALL Field_ParameterSetGetNode(solidDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                          & versionIdx,derivativeIdx,solidNode,componentIdx,solidNodePosition,err,error,*999)
                        CALL Field_ParameterSetGetNode(solidDependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                          & versionIdx,derivativeIdx,solidNode,componentIdx,previousSolidNodePosition,err,error,*999)
                        solidDelta=solidNodePosition-previousSolidNodePosition
                        IF(ABS(solidDelta)>ZERO_TOLERANCE) THEN
                          CALL Field_ParameterSetAddNode(fluidIndependentField,FIELD_U_VARIABLE_TYPE, &
                            & FIELD_MESH_DISPLACEMENT_SET_TYPE,versionIdx,derivativeIdx,fluidNode,componentIdx, &
                            & solidDelta,err,error,*999)
                        ENDIF
                      ENDDO !versionIdx
                    ENDDO !derivativeIdx
                  ENDDO !nodeIdx
                CASE DEFAULT
                  CALL FlagError("Interface geometric component does not have node based interpolation.",err,error,*999)
                END SELECT
              ENDDO !componentIdx
              CALL Field_ParameterSetDataGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
                & meshDisplacementValues,err,error,*999)
              CALL Field_VariableGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidGeometricVariable,err,error,*999)
              DO componentIdx=1,fluidGeometricVariable%NUMBER_OF_COMPONENTS
                NULLIFY(domain)
                CALL FieldVariable_DomainGet(fluidGeometricVariable,componentIdx,domain,err,error,*999)
                NULLIFY(domainTopology)
                CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
                NULLIFY(domainNodes)
                CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
                !Loop over the local nodes excluding the ghosts.
                DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                  DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                    DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                      localDOF=fluidGeometricVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                        & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(derivativeIdx)%versions(versionIdx)
                      CALL Field_ParameterSetAddLocalDOF(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                        & localDOF,meshDisplacementValues(localDOF),err,error,*999)
                    ENDDO !versionIdx
                  ENDDO !derivativeIdx
                ENDDO !nodeIdx
              ENDDO !componentIdx
              CALL Field_ParameterSetDataRestore(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
                & meshDisplacementValues,err,error,*999)
              CALL Field_ParameterSetUpdateStart(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetUpdateFinish(fluidGeometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
              !Now use displacement values to calculate velocity values
              alpha=1.0_DP/timeIncrement
              CALL Field_ParameterSetsCopy(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
                & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
            ELSE
              CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
            & " is not valid for a Finite Elasticity-NavierS tokes type of a multi physics problem class."
        END SELECT
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(problem%specification(2),"*",err,error))// &
          & " is not valid for a multi physics problem class."
      END SELECT
    CASE DEFAULT
      localError="Problem class "//TRIM(NumberToVString(problem%specification(1),"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_PreSolveALEUpdateMesh")
    RETURN
999 ERRORSEXITS("NavierStokes_PreSolveALEUpdateMesh",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PreSolveALEUpdateMesh

  !
  !================================================================================================================================
  !
  !>Update mesh parameters for Laplace problem
  SUBROUTINE NavierStokes_PreSolveALEUpdateParameters(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(FIELD_TYPE), POINTER :: independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: component_idx,deriv_idx,local_ny,node_idx,variable_idx,variable_type
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP), POINTER :: MESH_STIFF_VALUES(:)

    ENTERS("NavierStokes_PreSolveALEUpdateParameters",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%specification(1))
          CASE(PROBLEM_FLUID_MECHANICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
            CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE,PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
              IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                !Get the independent field for the ALE Navier-Stokes problem
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    NULLIFY(MESH_STIFF_VALUES)
                    CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                      IF(ASSOCIATED(EQUATIONS)) THEN
                        independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                        IF(ASSOCIATED(independentField)) THEN
                          DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                            variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          !Calculation of K values dependent on current mesh topology
                                          MESH_STIFF_VALUES(local_ny)=1.0_DP
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD, &
                                            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                            & MESH_STIFF_VALUES(local_ny),err,error,*999)
                                        END DO !deriv_idx
                                      END DO !node_idx
                                    END IF
                                  END IF
                                END IF
                              END DO !component_idx
                            END IF
                          END DO !variable_idx
                        ELSE
                          CALL FlagError("Independent field is not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations are not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    END IF
                    CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Solver equations are not associated.",err,error,*999)
                END IF
              ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(3),"*",err,error))// &
                & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(PROBLEM_MULTI_PHYSICS_CLASS)
            SELECT CASE(CONTROL_LOOP%PROBLEM%specification(2))
            CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
              SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
              CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
                & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
                IF(SOLVER%SOLVE_TYPE==SOLVER_LINEAR_TYPE) THEN
                  !Get the independent field for the ALE Navier-Stokes problem
                  SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                  IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                    SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                    IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                      EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                      NULLIFY(MESH_STIFF_VALUES)
                      CALL FIELD_PARAMETER_SET_DATA_GET(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                      IF(ASSOCIATED(EQUATIONS_SET)) THEN
                        EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
                        IF(ASSOCIATED(EQUATIONS)) THEN
                          independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
                          IF(ASSOCIATED(independentField)) THEN
                            DO variable_idx=1,EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                              variable_type=EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                              FIELD_VARIABLE=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
                              IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                                DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                  DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                  IF(ASSOCIATED(DOMAIN)) THEN
                                    IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                      DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                      IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                        !Loop over the local nodes excluding the ghosts.
                                        DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                          DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                            !Default to version 1 of each node derivative
                                            local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                              & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                            !Calculation of K values dependent on current mesh topology
                                            MESH_STIFF_VALUES(local_ny)=1.0_DP
                                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(EQUATIONS_SET%INDEPENDENT% &
                                              & INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,local_ny, &
                                              & MESH_STIFF_VALUES(local_ny),err,error,*999)
                                          END DO !deriv_idx
                                        END DO !node_idx
                                      END IF
                                    END IF
                                  END IF
                                END DO !component_idx
                              END IF
                            END DO !variable_idx
                          ELSE
                            CALL FlagError("Independent field is not associated.",err,error,*999)
                          END IF
                        ELSE
                          CALL FlagError("Equations are not associated.",err,error,*999)
                        END IF
                      ELSE
                        CALL FlagError("Equations set is not associated.",err,error,*999)
                      END IF
                      CALL Field_ParameterSetDataRestore(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,MESH_STIFF_VALUES,err,error,*999)
                    ELSE
                      CALL FlagError("Solver mapping is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Solver equations are not associated.",err,error,*999)
                  END IF
                ELSE IF(SOLVER%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                  CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
                END IF
              CASE DEFAULT
                localError="The third problem specification of "// &
                  & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(3),"*",err,error))// &
                  & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The second problem specification of "// &
                & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(2),"*",err,error))// &
                & " is not valid for NavierStokes_PreSolveALEUpdateParameters of a multi physics problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The first problem specification of "// &
              & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(1),"*",err,error))// &
              & " is not valid for NavierStokes_PreSolveALEUpdateParameters."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("NavierStokes_PreSolveALEUpdateParameters")
    RETURN
999 ERRORSEXITS("NavierStokes_PreSolveALEUpdateParameters",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PreSolveALEUpdateParameters

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    INTEGER(INTG), INTENT(OUT) :: ERR
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELDS_TYPE), POINTER :: Fields
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError,METHOD,VFileName,FILENAME
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,FileNameLength
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,START_TIME,STOP_TIME
    CHARACTER(20) :: FILE,OUTPUT_FILE

    NULLIFY(Fields)

    ENTERS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Navier-Stokes problem.",err,error,*999)
          END IF
          CALL SYSTEM('mkdir -p ./output')
          SELECT CASE(CONTROL_LOOP%PROBLEM%specification(3))
          CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                  FILENAME="./output/"//"STATIC_SOLUTION"
                  METHOD="FORTRAN"
                  IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                  ENDIF
                  Fields=>EQUATIONS_SET%REGION%FIELDS
                  CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                  NULLIFY(Fields)
                ENDDO
              END IF
            END IF

          CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)

            IF(SOLVER%GLOBAL_NUMBER==2) THEN
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
                        WRITE(OUTPUT_FILE,'("TimeStep_",I0)') CURRENT_LOOP_ITERATION
                        FILE=OUTPUT_FILE
                        METHOD="FORTRAN"
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                          !Use standard field IO routines (also only export nodes after first step as not a moving mesh case)
                          FileNameLength = LEN_TRIM(OUTPUT_FILE)
                          VFileName = OUTPUT_FILE(1:FileNameLength)
                          IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                            & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          Fields=>EQUATIONS_SET%REGION%FIELDS
                          CALL FIELD_IO_NODES_EXPORT(Fields,VFileName,METHOD,err,error,*999)
                          IF(CURRENT_LOOP_ITERATION==0) THEN
                            IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                              & CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export elements... ",err,error,*999)
                            CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,METHOD,err,error,*999)
                          END IF
                          NULLIFY(Fields)
                          IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          ENDIF
                        END IF
                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                            CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,err,error,*999)
                          END IF
                        END IF
                      END IF
                    END IF
                  END DO
                END IF
              END IF
            ENDIF

          CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
               & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,ERR,ERROR,*999)
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%PTR
                  SELECT CASE(EQUATIONS_SET%specification(3))
                  CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CURRENT_TIME<=STOP_TIME) THEN
                        WRITE(OUTPUT_FILE,'("TimeStep3D_",I0)') CURRENT_LOOP_ITERATION
                        FILE=OUTPUT_FILE
                        METHOD="FORTRAN"
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                          !Use standard field IO routines (also only export nodes after first step as not a moving mesh case)
                          FileNameLength = LEN_TRIM(OUTPUT_FILE)
                          VFileName = OUTPUT_FILE(1:FileNameLength)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          Fields=>EQUATIONS_SET%REGION%FIELDS
                          CALL FIELD_IO_NODES_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                          IF(CURRENT_LOOP_ITERATION==0) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export elements... ",ERR,ERROR,*999)
                            CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,METHOD,ERR,ERROR,*999)
                          END IF
                          NULLIFY(Fields)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                        END IF
                      END IF
                    END IF

                  CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                     & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CURRENT_TIME<=STOP_TIME) THEN
                        DEPENDENT_REGION=>EQUATIONS_SET%REGION
                        FILE=OUTPUT_FILE
                        FILENAME="TimeStep1D_"//TRIM(NUMBER_TO_VSTRING(CURRENT_LOOP_ITERATION,"*",ERR,ERROR))
                        METHOD="FORTRAN"
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",ERR,ERROR,*999)
                          CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                          ! Only export elements on first iteration (non-moving mesh case)
                          IF(CURRENT_LOOP_ITERATION==0) THEN
                            CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                          END IF
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,ERR,ERROR,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",ERR,ERROR,*999)
                          CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,ERR,ERROR,*999)
                        END IF
                      END IF
                    END IF
                  CASE DEFAULT
                    localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%specification(3),"*",err,error))// &
                      & " is not valid for a dynamic solver output for a multiscale Navier-Stokes problem."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                END DO
              END IF
            END IF

          CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
             & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)

            CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
              & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,err,error,*999)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
              IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                !Make sure the equations sets are up to date
                DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                  EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                  IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                    IF(CURRENT_TIME<=STOP_TIME) THEN
                      IF(CURRENT_LOOP_ITERATION<10) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                      ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                        WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                      END IF
                      DEPENDENT_REGION=>EQUATIONS_SET%REGION
                      FILE=OUTPUT_FILE
                      FILENAME="./output/"//"MainTime_"//TRIM(NumberToVString(CURRENT_LOOP_ITERATION,"*",err,error))
                      METHOD="FORTRAN"
                      IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                        IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                        ENDIF
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,err,error,*999)
                        CALL FIELD_IO_ELEMENTS_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,err,error,*999)
                        IF(CONTROL_LOOP%outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,err,error,*999)
                          CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                        ENDIF
                        CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_DIMENSIONS,err,error,*999)
                      END IF
                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                          & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,err,error,*999)
                        END IF
                      END IF
                    END IF
                  END IF
                END DO
              END IF
            END IF
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%specification(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Control loop is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1

  END SUBROUTINE NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Sets up analytic parameters and calls NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE to evaluate solutions to analytic problems
  SUBROUTINE NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: interpolatedPoint(:)
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: interpolationParameters(:)
    TYPE(FIELD_TYPE), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable,analyticVariable,materialsVariable
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: componentIdx,derivativeIdx,dimensionIdx,local_ny,nodeIdx,numberOfDimensions,variableIdx,variableType,I,J,K
    INTEGER(INTG) :: numberOfNodesXiCoord(3),elementIdx,en_idx,boundaryCount,analyticFunctionType,globalDerivativeIndex,versionIdx
    INTEGER(INTG) :: boundaryConditionsCheckVariable,numberOfXi,nodeNumber,userNodeNumber,localDof,globalDof
    INTEGER(INTG) :: parameterIdx,numberOfParameters
    REAL(DP) :: TIME,VALUE,X(3),xiCoordinates(3),initialValue,T_COORDINATES(20,3),nodeAnalyticParameters(10)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)

    ENTERS("NavierStokes_BoundaryConditionsAnalyticCalculate",err,error,*999)

    boundaryCount=0
    xiCoordinates(3)=0.0_DP

    IF(ASSOCIATED(equationsSet)) THEN
      IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
        dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN
            ! Geometric parameters
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            NULLIFY(geometricVariable)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
            NULLIFY(geometricParameters)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,geometricParameters, &
              & err,error,*999)
            ! Analytic parameters
            analyticFunctionType=equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            analyticField=>equationsSet%ANALYTIC%ANALYTIC_FIELD
            NULLIFY(analyticVariable)
            NULLIFY(analyticParameters)
            IF(ASSOCIATED(analyticField)) THEN
              CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & analyticParameters,err,error,*999)
            END IF
            ! Materials parameters
            NULLIFY(materialsField)
            NULLIFY(materialsVariable)
            NULLIFY(materialsParameters)
            IF(ASSOCIATED(equationsSet%MATERIALS)) THEN
              materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
              CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & materialsParameters,err,error,*999)
            END IF
            TIME=equationsSet%ANALYTIC%ANALYTIC_TIME
            ! Interpolation parameters
            NULLIFY(interpolationParameters)
            CALL FIELD_INTERPOLATION_PARAMETERS_INITIALISE(geometricField,interpolationParameters,err,error,*999)
            NULLIFY(interpolatedPoint)
            CALL FIELD_INTERPOLATED_POINTS_INITIALISE(interpolationParameters,interpolatedPoint,err,error,*999)
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    IF(ASSOCIATED(boundaryConditions)) THEN
      DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
        variableType=dependentField%VARIABLES(variableIdx)%VARIABLE_TYPE
        fieldVariable=>dependentField%VARIABLE_TYPE_MAP(variableType)%ptr
        IF(ASSOCIATED(fieldVariable)) THEN
          IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_ANALYTIC_VALUES_SET_TYPE)%ptr)) &
            & CALL FIELD_PARAMETER_SET_CREATE(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          DO componentIdx=1,fieldVariable%NUMBER_OF_COMPONENTS
            boundaryCount=0
            IF(fieldVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              domain=>fieldVariable%COMPONENTS(componentIdx)%DOMAIN
              IF(ASSOCIATED(domain)) THEN
                IF(ASSOCIATED(domain%TOPOLOGY)) THEN
                  domainNodes=>domain%TOPOLOGY%NODES
                  IF(ASSOCIATED(domainNodes)) THEN
                    !Loop over the local nodes excluding the ghosts.
                    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
                      nodeNumber = domainNodes%NODES(nodeIdx)%local_number
                      userNodeNumber = domainNodes%NODES(nodeIdx)%user_number
                      elementIdx=domain%topology%nodes%nodes(nodeNumber)%surrounding_elements(1)
                      CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementIdx, &
                        & interpolationParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                      en_idx=0
                      xiCoordinates=0.0_DP
                      numberOfXi=domain%topology%elements%elements(elementIdx)%basis%number_of_xi
                      numberOfNodesXiCoord(1)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(1)
                      IF(numberOfXi>1) THEN
                        numberOfNodesXiCoord(2)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(2)
                      ELSE
                        numberOfNodesXiCoord(2)=1
                      END IF
                      IF(numberOfXi>2) THEN
                        numberOfNodesXiCoord(3)=domain%topology%elements%elements(elementIdx)%basis%number_of_nodes_xic(3)
                      ELSE
                        numberOfNodesXiCoord(3)=1
                      END IF

                      SELECT CASE(analyticFunctionType)
                      ! --- Calculate analytic profile for validation ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                        IF(variableIdx < 3) THEN
                          ! Get geometric position info for this node
                          DO dimensionIdx=1,numberOfDimensions
                            local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                              & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                            X(dimensionIdx)=geometricParameters(local_ny)
                          END DO !dimensionIdx
                          DO derivativeIdx=1,domainNodes%NODES(nodeNumber)%NUMBER_OF_DERIVATIVES
                            globalDerivativeIndex=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)% &
                              & GLOBAL_DERIVATIVE_INDEX
                            DO versionIdx=1,domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions
                              CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                & globalDerivativeIndex,componentIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                & analyticParameters,materialsParameters,VALUE,err,error,*999)
                              local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                            END DO !versionIdx
                          END DO !derivativeIdx
                        END IF ! variableIdx < 3

                      ! --- Set velocity boundary conditions with analytic value ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                        ! Get geometric position info for this node
                        DO dimensionIdx=1,numberOfDimensions
                          local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeNumber)%DERIVATIVES(1)%VERSIONS(1)
                          X(dimensionIdx)=geometricParameters(local_ny)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivativeIdx=1,domainNodes%NODES(nodeNumber)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          IF(componentIdx<=numberOfXi .OR. &
                            &  analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                            DO versionIdx=1,domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions
                              ! Get global and local dof indices
                              CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField,variableType,versionIdx,derivativeIdx, &
                               & userNodeNumber,componentIdx,localDof,globalDof,err,error,*999)
                              IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                                CALL FIELD_NUMBER_OF_COMPONENTS_GET(analyticField,FIELD_U_VARIABLE_TYPE, &
                                 & numberOfParameters,err,error,*999)
                                DO parameterIdx=1,numberOfParameters
                                  ! populate nodeAnalyticParameters
                                  CALL Field_ParameterSetGetLocalNode(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                   & versionIdx,derivativeIdx,nodeNumber,parameterIdx,nodeAnalyticParameters(parameterIdx), &
                                   & err,error,*999)
                                END DO
                                CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                  & globalDerivativeIndex,componentIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                  & nodeAnalyticParameters,materialsParameters,VALUE,err,error,*999)
                              ELSE
                                CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                  & globalDerivativeIndex,componentIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                                  & analyticParameters,materialsParameters,VALUE,err,error,*999)
                              END IF
                              ! update analytic field values
                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,localDof,VALUE,err,error,*999)
                              IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                                IF(domainNodes%NODES(nodeNumber)%BOUNDARY_NODE) THEN
                                  CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable, &
                                   & boundaryConditionsVariable,err,error,*999)
                                  IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                                    boundaryConditionsCheckVariable=boundaryConditionsVariable% &
                                     & CONDITION_TYPES(globalDof)
                                    ! update dependent field values if fixed inlet or pressure BC
                                    IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET .OR. &
                                     & boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_PRESSURE) THEN
                                      ! Set DOF values
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                       & FIELD_VALUES_SET_TYPE,localDof,VALUE,err,error,*999)
                                    ELSE IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_PRESSURE) THEN
                                      ! ! Set neumann boundary pressure value on pressure nodes
                                      ! CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                                      !  & FIELD_PRESSURE_VALUES_SET_TYPE,1,1,nodeNumber,componentIdx,VALUE,err,error,*999)
                                    END IF
                                  END IF
                                END IF
                              END IF
                            END DO !versionIdx
                          END IF
                        END DO !derivativeIdx

                      ! --- Set Flow rate boundary conditions with analytic value ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                         & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                        ! Get geometric position info for this node
                        DO dimensionIdx=1,numberOfDimensions
                          local_ny=geometricVariable%COMPONENTS(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                            & NODES(nodeIdx)%DERIVATIVES(1)%VERSIONS(1)
                          X(dimensionIdx)=geometricParameters(local_ny)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivativeIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          IF(componentIdx==1 .AND. variableType==FIELD_U_VARIABLE_TYPE) THEN
                            DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                              local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                                & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                              IF(domainNodes%NODES(nodeIdx)%BOUNDARY_NODE) THEN
                                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable, &
                                 & boundaryConditionsVariable,err,error,*999)
                                IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                                  boundaryConditionsCheckVariable=boundaryConditionsVariable%CONDITION_TYPES(local_ny)
                                  IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                                    CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                                      & globalDerivativeIndex,componentIdx,numberOfXi,fieldVariable%NUMBER_OF_COMPONENTS, &
                                      & analyticParameters,materialsParameters,VALUE,err,error,*999)
                                    !If we are a boundary node then set the analytic value on the boundary
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                      & FIELD_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                  ELSE
                                    CALL Field_ParameterSetGetLocalNode(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
                                      & versionIdx,derivativeIdx,nodeIdx,componentIdx,VALUE,err,error,*999)
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                  END IF
                                END IF
                              END IF
                            END DO !versionIdx
                          END IF
                        END DO !derivativeIdx

                      ! --- Legacy unit shape testing types ---
                      CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4,  &
                       & EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
                        !Quad/Hex
                        !\todo: Use boundary flag
                        IF(domain%topology%elements%maximum_number_of_element_parameters==4.AND.numberOfDimensions==2.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==9.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==16.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==8.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==27.OR. &
                          & domain%topology%elements%maximum_number_of_element_parameters==64) THEN
                          DO K=1,numberOfNodesXiCoord(3)
                            DO J=1,numberOfNodesXiCoord(2)
                              DO I=1,numberOfNodesXiCoord(1)
                                en_idx=en_idx+1
                                IF(domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                                xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXiCoord(1)-1))
                              END DO
                                IF(domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                                xiCoordinates(1)=0.0_DP
                                xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXiCoord(2)-1))
                            END DO
                            IF(domain%topology%elements%elements(elementIdx)%element_nodes(en_idx)==nodeIdx) EXIT
                            xiCoordinates(1)=0.0_DP
                            xiCoordinates(2)=0.0_DP
                            IF(numberOfNodesXiCoord(3)/=1) THEN
                              xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXiCoord(3)-1))
                            END IF
                          END DO
                          CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,xiCoordinates, &
                            & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                        !Tri/Tet
                        !\todo: Use boundary flag
                        ELSE
                          IF(domain%topology%elements%maximum_number_of_element_parameters==3) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                          ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==6) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                            T_COORDINATES(4,1:2)=[0.5_DP,0.5_DP]
                            T_COORDINATES(5,1:2)=[1.0_DP,0.5_DP]
                            T_COORDINATES(6,1:2)=[0.5_DP,1.0_DP]
                          ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==10.AND. &
                            & numberOfDimensions==2) THEN
                            T_COORDINATES(1,1:2)=[0.0_DP,1.0_DP]
                            T_COORDINATES(2,1:2)=[1.0_DP,0.0_DP]
                            T_COORDINATES(3,1:2)=[1.0_DP,1.0_DP]
                            T_COORDINATES(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                          ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==4) THEN
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                          ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==10.AND. &
                            & numberOfDimensions==3) THEN
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                            T_COORDINATES(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                            T_COORDINATES(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                            T_COORDINATES(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                            T_COORDINATES(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                          ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==20) THEN
                            T_COORDINATES(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                            T_COORDINATES(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                            T_COORDINATES(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                            T_COORDINATES(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                            T_COORDINATES(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                            T_COORDINATES(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                            T_COORDINATES(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                          END IF
                          DO K=1,domain%topology%elements%maximum_number_of_element_parameters
                            IF(domain%topology%elements%elements(elementIdx)%element_nodes(K)==nodeIdx) EXIT
                          END DO
                          IF(numberOfDimensions==2) THEN
                            CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:2), &
                              & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                          ELSE IF(numberOfDimensions==3) THEN
                            CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,T_COORDINATES(K,1:3), &
                              & interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                          END IF
                        END IF
                        X=0.0_DP
                        DO dimensionIdx=1,numberOfDimensions
                          X(dimensionIdx)=interpolatedPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(dimensionIdx,1)
                        END DO !dimensionIdx
                        !Loop over the derivatives
                        DO derivativeIdx=1,domainNodes%NODES(nodeIdx)%NUMBER_OF_DERIVATIVES
                          globalDerivativeIndex=domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)% &
                            & GLOBAL_DERIVATIVE_INDEX
                          CALL NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(analyticFunctionType,X,TIME,variableType, &
                            & globalDerivativeIndex,componentIdx,numberOfDimensions,fieldVariable%NUMBER_OF_COMPONENTS, &
                            & analyticParameters,materialsParameters,VALUE,err,error,*999)
                          DO versionIdx=1,domainNodes%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions
                            local_ny=fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
                              & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
                            CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                            IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                              IF(domainNodes%NODES(nodeIdx)%BOUNDARY_NODE) THEN
                                !If we are a boundary node then set the analytic value on the boundary
                                CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType, &
                                  & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                ! \todo: This is just a workaround for linear pressure fields in simplex element components
                                IF(componentIdx>numberOfDimensions) THEN
                                  IF(domain%topology%elements%maximum_number_of_element_parameters==3) THEN
                                    IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_1.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_2.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_3.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_4.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_TWO_DIM_5) THEN
                                      IF(-0.001_DP<X(1).AND.X(1)<0.001_DP.AND.-0.001_DP<X(2).AND.X(2)<0.001_DP.OR. &
                                        &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.-0.001_DP<X(2).AND. &
                                        & X(2)<0.001_DP.OR. &
                                        &  10.0_DP-0.001_DP<X(1).AND.X(1)<10.0_DP+0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<10.0_DP+0.001_DP.OR. &
                                        &  -0.001_DP<X(1).AND.X(1)<0.001_DP.AND.10.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<10.0_DP+0.001_DP) THEN
                                          CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField, &
                                            & variableType,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                          boundaryCount=boundaryCount+1
                                      END IF
                                    END IF
                                  ELSE IF(domain%topology%elements%maximum_number_of_element_parameters==4.AND. &
                                    & numberOfDimensions==3) THEN
                                    IF(analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_1.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_2.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_3.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_4.OR. &
                                      & analyticFunctionType==EQUATIONS_SET_STOKES_EQUATION_THREE_DIM_5) THEN
                                      IF(-5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(3).AND.X(3)<-5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & -5.0_DP-0.001_DP<X(1).AND.X(1)<-5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<5.0_DP+0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP.OR. &
                                        & 5.0_DP-0.001_DP<X(1).AND.X(1)<5.0_DP+0.001_DP.AND.-5.0_DP-0.001_DP<X(2).AND. &
                                        & X(2)<-5.0_DP+ 0.001_DP.AND.5.0_DP-0.001_DP<X(3).AND.X(3)<5.0_DP+0.001_DP) THEN
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField, &
                                          & variableType,local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                        boundaryCount=boundaryCount+1
                                      END IF
                                    END IF
                                    ! \todo: This is how it should be if adjacent elements would be working
                                  ELSE IF(boundaryCount==0) THEN
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType,&
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                    boundaryCount=boundaryCount+1
                                  END IF
                                END IF
                              ELSE
                                !Set the initial condition.
                                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variableType, &
                                  & FIELD_VALUES_SET_TYPE,local_ny,initialValue,err,error,*999)
                              END IF
                            END IF
                          END DO !versionIdx
                        END DO !derivativeIdx

                      CASE DEFAULT
                        localError="Analytic Function Type "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                          & " is not yet implemented for a Navier-Stokes problem."
                        CALL FlagError(localError,err,error,*999)
                      END SELECT

                    END DO !nodeIdx
                  ELSE
                    CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Domain topology is not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Domain is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
            END IF
          END DO !componentIdx
          ! Update ghost and boundary values from local
          IF(ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_PRESSURE_VALUES_SET_TYPE)%PTR)) THEN
            CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
              & err,error,*999)
            CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE, &
              & err,error,*999)
          END IF
          CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
        ELSE
          CALL FlagError("Field variable is not associated.",err,error,*999)
        END IF
      END DO !variableIdx
      CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & geometricParameters,err,error,*999)
      CALL FIELD_INTERPOLATED_POINTS_FINALISE(interpolatedPoint,err,error,*999)
      CALL FIELD_INTERPOLATION_PARAMETERS_FINALISE(interpolationParameters,err,error,*999)
    ELSE
      CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("NavierStokes_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("NavierStokes_BoundaryConditionsAnalyticCalculate")
    RETURN 1

  END SUBROUTINE NavierStokes_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !
  !>Calculates the various analytic values for NSE examples with exact solutions
  SUBROUTINE NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE(ANALYTIC_FUNCTION_TYPE,X,TIME,VARIABLE_TYPE,GLOBAL_DERIV_INDEX, &
    & componentNumber,NUMBER_OF_DIMENSIONS,NUMBER_OF_COMPONENTS,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: X(:) !<X(dimension_idx). The geometric position to evaluate at (includes Y,Z for higher dim problems)
    REAL(DP), INTENT(IN) :: TIME !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: GLOBAL_DERIV_INDEX !<The global derivative to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_DIMENSIONS !<The number of geometric dimensions
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_COMPONENTS !<The number of components for the dependent field
    REAL(DP), INTENT(IN) :: ANALYTIC_PARAMETERS(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: MATERIALS_PARAMETERS(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,n,m
    REAL(DP) :: L_PARAM,H_PARAM,U_PARAM,P_PARAM,MU_PARAM,NU_PARAM,RHO_PARAM,INTERNAL_TIME,CURRENT_TIME,K_PARAM
    REAL(DP) :: amplitude,yOffset,period,phaseShift,frequency,s,startTime,stopTime,tt,tmax,Qo
    REAL(DP) :: componentCoeff(4),delta(300),t(300),q(300)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE",err,error,*999)

    !\todo: Introduce user-defined or default values instead for density and viscosity
    INTERNAL_TIME=TIME
    CURRENT_TIME=TIME

     SELECT CASE(ANALYTIC_FUNCTION_TYPE)

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
       !For fully developed 2D laminar flow through a channel, NSE should yield a parabolic profile,
       !U = Umax(1-y^2/H^2), Umax = (-dP/dx)*(H^2/(2*MU)), Umax = (3/2)*Umean
       !Note: assumes a flat inlet profile (U_PARAM = Umean).
       !Nonlinear terms from NSE will effectively be 0 for Poiseuille flow
       IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
         MU_PARAM = MATERIALS_PARAMETERS(1)
         RHO_PARAM = MATERIALS_PARAMETERS(2)
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           L_PARAM = ANALYTIC_PARAMETERS(1) ! channel length in x-direction
           H_PARAM = ANALYTIC_PARAMETERS(2) ! channel height in y-direction
           U_PARAM = ANALYTIC_PARAMETERS(3) ! mean (inlet) velocity
           P_PARAM = ANALYTIC_PARAMETERS(4) ! pressure value at outlet
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(componentNumber==1) THEN
               !calculate u
               VALUE=(3.0_DP/2.0_DP)*U_PARAM*(1.0_DP-((X(2)-H_PARAM)**2)/(H_PARAM**2))
             ELSE IF(componentNumber==2) THEN
               !calculate v
               VALUE=0.0_DP
             ELSE IF(componentNumber==3) THEN
               !calculate p
               VALUE = (3.0_DP*MU_PARAM*U_PARAM*(X(1)-L_PARAM))/(H_PARAM**2)+P_PARAM
             ELSE
               CALL FlagError("Not implemented.",ERR,ERROR,*999)
             END IF
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
               & " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE( NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
               & " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       ELSE
         localError="The number of components does not correspond to the number of dimensions."
         CALL FlagError(localError,ERR,ERROR,*999)
       END IF

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
       !Exact solution to 2D laminar, dynamic, nonlinear Taylor-Green vortex decay
       IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
         MU_PARAM = MATERIALS_PARAMETERS(1)
         RHO_PARAM = MATERIALS_PARAMETERS(2)
         NU_PARAM = MU_PARAM/RHO_PARAM ! kinematic viscosity
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           U_PARAM = ANALYTIC_PARAMETERS(1) ! characteristic velocity (initial amplitude)
           L_PARAM = ANALYTIC_PARAMETERS(2) ! length scale for square
           K_PARAM = 2.0_DP*PI/L_PARAM   ! scale factor for equations
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(componentNumber==1) THEN
               !calculate u
               VALUE=-1.0_DP*U_PARAM*COS(K_PARAM*X(1))*SIN(K_PARAM*X(2))*EXP(-2.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME)
             ELSE IF(componentNumber==2) THEN
               !calculate v
               VALUE=U_PARAM*SIN(K_PARAM*X(1))*COS(K_PARAM*X(2))*EXP(-2.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME)
             ELSE IF(componentNumber==3) THEN
               !calculate p
               VALUE =-1.0_DP*(U_PARAM**2)*(RHO_PARAM/4.0_DP)*(COS(2.0_DP*K_PARAM*X(1))+ &
                 & COS(2.0_DP*K_PARAM*X(2)))*(EXP(-4.0_DP*(K_PARAM**2)*NU_PARAM*CURRENT_TIME))
             ELSE
               CALL FlagError("Not implemented.",ERR,ERROR,*999)
             END IF
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
               & " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE( NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
               & " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       ELSE
         localError="The number of components does not correspond to the number of dimensions."
         CALL FlagError(localError,ERR,ERROR,*999)
       END IF

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA)
       SELECT CASE(NUMBER_OF_DIMENSIONS)
       CASE(1)
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(componentNumber==1) THEN
               !Input function
               period = 800
               tt=MOD(TIME,period)
               tmax=150.0_DP
               Qo=100000.0_DP
               VALUE=(Qo*tt/(tmax**2.0_DP))*EXP(-(tt**2.0_DP)/(2.0_DP*(tmax**2.0_DP)))
             ELSE
               CALL FlagError("Incorrect component specification for Aorta flow rate waveform ",ERR,ERROR,*999)
             END IF
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
           ! Do nothing
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       CASE DEFAULT
         localError="Aorta flowrate waveform for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
           & " dimension problem has not yet been implemented."
         CALL FlagError(localError,ERR,ERROR,*999)
       END SELECT

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
       SELECT CASE(NUMBER_OF_DIMENSIONS)
       CASE(1)
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(componentNumber==1) THEN
               !Olufsen Aorta
               t(1)= 0.0011660 ; q(1)= 17.39051
               t(2)= 0.0215840 ; q(2)= 10.41978
               t(3)= 0.0340860 ; q(3)= 18.75892
               t(4)= 0.0731370 ; q(4)= 266.3842
               t(5)= 0.0857710 ; q(5)= 346.3755
               t(6)= 0.1029220 ; q(6)= 413.8419
               t(7)= 0.1154270 ; q(7)= 424.2680
               t(8)= 0.1483530 ; q(8)= 429.1147
               t(9)= 0.1698860 ; q(9)= 411.0127
               t(10)= 0.220794 ; q(10)= 319.151
               t(11)= 0.264856 ; q(11)= 207.816
               t(12)= 0.295415 ; q(12)= 160.490
               t(13)= 0.325895 ; q(13)= 70.0342
               t(14)= 0.346215 ; q(14)= 10.1939
               t(15)= 0.363213 ; q(15)= -5.1222
               t(16)= 0.383666 ; q(16)= 6.68963
               t(17)= 0.405265 ; q(17)= 24.0659
               t(18)= 0.427988 ; q(18)= 35.8762
               t(19)= 0.455272 ; q(19)= 58.8137
               t(20)= 0.477990 ; q(20)= 67.8414
               t(21)= 0.502943 ; q(21)= 57.3893
               t(22)= 0.535816 ; q(22)= 33.7142
               t(23)= 0.577789 ; q(23)= 20.4676
               t(24)= 0.602753 ; q(24)= 16.2763
               t(25)= 0.639087 ; q(25)= 22.5119
               t(26)= 0.727616 ; q(26)= 18.9721
               t(27)= 0.783235 ; q(27)= 18.9334
               t(28)= 0.800000 ; q(28)= 16.1121

               !Initialize variables
               period = 800
               m=1
               n=28
               !Compute derivation
               DO i=1,n-1
                 delta(i)=(q(i+1)-q(i))/(t(i+1)-t(i))
               END DO
               delta(n)=delta(n-1)+(delta(n-1)-delta(n-2))/(t(n-1)-t(n-2))*(t(n)-t(n-1))
               !Find subinterval
               DO j=1,n-1
                 IF(t(j) <= (TIME/period)) THEN
                   m=j
                 END IF
               END DO
               !Evaluate interpolant
               s=(TIME/period)-t(m)
               VALUE=(q(m)+s*delta(m))
             ELSE
               CALL FlagError("Incorrect component specification for Olufsen flow rate waveform ",ERR,ERROR,*999)
             END IF
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
           ! Do nothing
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       CASE DEFAULT
         localError="Olufsen flowrate waveform for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
           & " dimension problem has not yet been implemented."
         CALL FlagError(localError,ERR,ERROR,*999)
       END SELECT

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
       ! Returns a sinusoidal value for boundary nodes
       SELECT CASE(NUMBER_OF_DIMENSIONS)
       CASE(2,3)
         componentCoeff(1) = ANALYTIC_PARAMETERS(1)
         componentCoeff(2) = ANALYTIC_PARAMETERS(2)
         componentCoeff(3) = ANALYTIC_PARAMETERS(3)
         componentCoeff(4) = ANALYTIC_PARAMETERS(4)
         amplitude = ANALYTIC_PARAMETERS(5)
         yOffset = ANALYTIC_PARAMETERS(6)
         frequency = ANALYTIC_PARAMETERS(7)
         phaseShift = ANALYTIC_PARAMETERS(8)
         startTime = ANALYTIC_PARAMETERS(9)
         stopTime = ANALYTIC_PARAMETERS(10)
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(CURRENT_TIME > startTime - ZERO_TOLERANCE .AND. &
               &  CURRENT_TIME < stopTime + ZERO_TOLERANCE) THEN
               VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*CURRENT_TIME+phaseShift))
             ELSE
               VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*stopTime+phaseShift))
             END IF
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       CASE DEFAULT
         localError="Sinusoidal analytic types for "//TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",ERR,ERROR))// &
           & " dimensional problems have not yet been implemented."
         CALL FlagError(localError,ERR,ERROR,*999)
       END SELECT

     CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
       IF(NUMBER_OF_DIMENSIONS==1.AND.NUMBER_OF_COMPONENTS==3) THEN
         !Polynomial function
         SELECT CASE(VARIABLE_TYPE)
         CASE(FIELD_U_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             IF(componentNumber==1) THEN
               !calculate Q
               VALUE=X(1)**2/10.0_DP**2
             ELSE IF(componentNumber==2) THEN
               !calculate A
               VALUE=X(1)**2/10.0_DP**2
             ELSE IF(componentNumber==3) THEN
               !calculate P
               VALUE=X(1)**2/10.0_DP**2
             ELSE
               CALL FlagError("Not implemented.",ERR,ERROR,*999)
             END IF
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE(FIELD_DELUDELN_VARIABLE_TYPE)
           SELECT CASE(GLOBAL_DERIV_INDEX)
           CASE(NO_GLOBAL_DERIV)
             VALUE= 0.0_DP
           CASE(GLOBAL_DERIV_S1)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE(GLOBAL_DERIV_S1_S2)
             CALL FlagError("Not implemented.",ERR,ERROR,*999)
           CASE DEFAULT
             localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
               & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
               & " is invalid."
             CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         CASE DEFAULT
           localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
             & " is invalid."
           CALL FlagError(localError,ERR,ERROR,*999)
         END SELECT
       ELSE
         localError="The number of components does not correspond to the number of dimensions."
         CALL FlagError(localError,ERR,ERROR,*999)
       END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Polynomial function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(VARIABLE_TYPE)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=2.0_DP/3.0_DP*X(1)*(3.0_DP*MU_PARAM*10.0_DP**2-RHO_PARAM*X(1)**2*X(2))/(10.0_DP ** 4)
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   VALUE= 0.0_DP
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Exponential function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(VARIABLE_TYPE)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE= EXP((X(1)-X(2))/10.0_DP)
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE= 2.0_DP*MU_PARAM/10.0_DP*EXP((X(1)-X(2))/10.0_DP)
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE= 0.0_DP
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE= 0.0_DP
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE= 0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Sine and cosine function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(VARIABLE_TYPE)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=SIN(2.0_DP*PI*X(1)/10.0_DP)*SIN(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=4.0_DP*MU_PARAM*PI/10.0_DP*SIN(2.0_DP*PI*X(2)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)+ &
                       & 0.5_DP*RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*COS(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_index,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=16.0_DP*MU_PARAM*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/ 10.0_DP)*cos(2.0_DP*PI*X(1)/10.0_DP)
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
         IF(NUMBER_OF_DIMENSIONS==2.AND.NUMBER_OF_COMPONENTS==3) THEN
           !Taylor-Green vortex solution
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(VARIABLE_TYPE)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=SIN(X(1)/10.0_DP*2.0_DP*PI)*COS(X(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
                     VALUE=SIN(X(1)/10.0_DP*PI)*COS(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
!                      VALUE=SIN(X(1))*COS(X(2))
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=-COS(X(1)/10.0_DP*2.0_DP*PI)*SIN(X(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
                     VALUE=-COS(X(1)/10.0_DP*PI)*SIN(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
!                      VALUE=-COS(X(1))*SIN(X(2))
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*2.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*2.0_DP*PI))* &
                       & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
                     VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*PI))* &
                       & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
!                      VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1))+COS(2.0_DP*X(2)))
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Polynomial function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=X(2)**2/10.0_DP**2+X(3)**2/10.0_DP**2
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=X(1)**2/10.0_DP**2+X(3)**2/10.0_DP** 2
                   ELSE IF(componentNumber==3) THEN
                     !calculate w
                     VALUE=X(1)**2/10.0_DP**2+X(2)**2/10.0_DP** 2
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=2.0_DP/3.0_DP*X(1)*(6.0_DP*MU_PARAM*10.0_DP**2-RHO_PARAM*X(2)*X(1)**2-3.0_DP* &
                       & RHO_PARAM*X(2)* &
                       & X(3)**2-RHO_PARAM*X(3)*X(1)**2-3.0_DP*RHO_PARAM*X(3)*X(2)**2)/(10.0_DP**4)
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   VALUE=0.0_DP
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Exponential function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(3)-X(1))/10.0_DP)
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=EXP((X(1)-X(2))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(componentNumber==3) THEN
                     !calculate w
                     VALUE=EXP((X(3)-X(1))/10.0_DP)+EXP((X(2)-X(3))/10.0_DP)
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=1.0_DP/10.0_DP*(2.0_DP*MU_PARAM*EXP((X(1)-X(2))/10.0_DP)- &
                       & 2.0_DP*MU_PARAM*EXP((X(3)-X(1))/10.0_DP)+RHO_PARAM*10.0_DP*EXP((X(1)-X(3))/10.0_DP)+ &
                       & RHO_PARAM*10.0_DP*EXP((X(2)-X(1))/10.0_DP))
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(1)-X(2))+EXP(X(2)-X(3)))
                   ELSE IF(componentNumber==3) THEN
                     !calculate w
                     VALUE=-2.0_DP*MU_PARAM*(2.0_DP*EXP(X(3)-X(1))+EXP(X(2)-X(3)))
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Sine/cosine function
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=sin(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=2.0_DP*cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)*cos(2.0_DP*PI*X(2)/10.0_DP)
                   ELSE IF(componentNumber==3) THEN
                     !calculate w
                     VALUE=-cos(2.0_DP*PI*X(1)/10.0_DP)*sin(2.0_DP*PI*X(2)/10.0_DP)*cos(2.0_DP*PI*X(3)/10.0_DP)
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=-COS(2.0_DP*PI*X(1)/10.0_DP)*(-12.0_DP*MU_PARAM*PI*SIN(2.0_DP*PI*X(2)/10.0_DP)* &
                       & SIN(2.0_DP*PI*X(3)/10.0_DP)-RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP+ &
                       & 2.0_DP*RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(3)/10.0_DP)**2- &
                       & RHO_PARAM*COS(2.0_DP*PI*X(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*X(2)/10.0_DP)**2)/10.0_DP/2.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(GLOBAL_DERIV_INDEX)
                  CASE(NO_GLOBAL_DERIV)
                    IF(componentNumber==1) THEN
                      !calculate u
                      VALUE=0.0_DP
                    ELSE IF(componentNumber==2) THEN
                      !calculate v
                      VALUE=36*MU_PARAM*PI**2/10.0_DP**2*cos(2.0_DP*PI*X(2)/10.0_DP)*sin(2.0_DP*PI*X(3)/10.0_DP)* &
                        & cos(2.0_DP*PI*X(1)/10.0_DP)
                    ELSE IF(componentNumber==3) THEN
                      !calculate w
                      VALUE=0.0_DP
                    ELSE IF(componentNumber==4) THEN
                      !calculate p
                      VALUE=0.0_DP
                    ELSE
                      CALL FlagError("Not implemented.",ERR,ERROR,*999)
                    END IF
                  CASE(GLOBAL_DERIV_S1)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(GLOBAL_DERIV_S2)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE(GLOBAL_DERIV_S1_S2)
                    CALL FlagError("Not implemented.",ERR,ERROR,*999)
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                      & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                      & " is invalid."
                    CALL FlagError(localError,ERR,ERROR,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                  & " is invalid."
                CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          ELSE
            localError="The number of components does not correspond to the number of dimensions."
            CALL FlagError(localError,ERR,ERROR,*999)
          END IF

       CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
         IF(NUMBER_OF_DIMENSIONS==3.AND.NUMBER_OF_COMPONENTS==4) THEN
           !Taylor-Green vortex solution
           MU_PARAM = MATERIALS_PARAMETERS(1)
           RHO_PARAM = MATERIALS_PARAMETERS(2)
           SELECT CASE(variable_type)
             CASE(FIELD_U_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=SIN(X(1)/10.0_DP*PI)*COS(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=-COS(X(1)/10.0_DP*PI)*SIN(X(2)/10.0_DP*PI)*EXP(-2.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
                   ELSE IF(componentNumber==3) THEN
                     !calculate v
                     VALUE=0.0_DP
!                      VALUE=-COS(X(1))*SIN(X(2))
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1)/10.0_DP*PI)+COS(2.0_DP*X(2)/10.0_DP*PI))* &
                       & EXP(-4.0_DP*MU_PARAM/RHO_PARAM*CURRENT_TIME)
!                      VALUE=RHO_PARAM/4.0_DP*(COS(2.0_DP*X(1))+COS(2.0_DP*X(2)))
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE(FIELD_DELUDELN_VARIABLE_TYPE)
               SELECT CASE(GLOBAL_DERIV_INDEX)
                 CASE(NO_GLOBAL_DERIV)
                   IF(componentNumber==1) THEN
                     !calculate u
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==2) THEN
                     !calculate v
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==3) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE IF(componentNumber==4) THEN
                     !calculate p
                     VALUE=0.0_DP
                   ELSE
                     CALL FlagError("Not implemented.",ERR,ERROR,*999)
                   END IF
                 CASE(GLOBAL_DERIV_S1)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE(GLOBAL_DERIV_S1_S2)
                   CALL FlagError("Not implemented.",ERR,ERROR,*999)
                 CASE DEFAULT
                   localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                     & GLOBAL_DERIV_INDEX,"*",ERR,ERROR))// &
                     & " is invalid."
                   CALL FlagError(localError,ERR,ERROR,*999)
               END SELECT
             CASE DEFAULT
               localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",ERR,ERROR))// &
                 & " is invalid."
               CALL FlagError(localError,ERR,ERROR,*999)
           END SELECT
         ELSE
           localError="The number of components does not correspond to the number of dimensions."
           CALL FlagError(localError,ERR,ERROR,*999)
         END IF
        CASE DEFAULT
          localError="The analytic function type of "// &
            & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(localError,ERR,ERROR,*999)
      END SELECT
    EXITS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 ERRORSEXITS("NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE",err,error)
    RETURN 1

  END SUBROUTINE NAVIER_STOKES_ANALYTIC_FUNCTIONS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Update SUPG parameters for Navier-Stokes equation
  SUBROUTINE NavierStokes_ResidualBasedStabilisation(equationsSet,elementNumber,gaussNumber,mu,rho,jacobianFlag,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number
    INTEGER(INTG), INTENT(IN) :: gaussNumber !<The gauss point number
    REAL(DP), INTENT(IN) :: mu !<The dynamic viscosity
    REAL(DP), INTENT(IN) :: rho !<The fluid density
    LOGICAL, INTENT(IN) ::  jacobianFlag !<Flag indicating whether this was called from the jacobian or residual evaluation routine
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: basisVelocity,basisPressure
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,equationsSetField,geometricField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureVelocity,quadraturePressure
    INTEGER(INTG) :: fieldVariableType,meshComponent1,meshComponent2
    INTEGER(INTG) :: numberOfDimensions
    INTEGER(INTG) :: i,j,k,l,mhs,nhs,ms,ns,nh,mh,nj,ni,pressureIndex
    INTEGER(INTG) :: numberOfElementParameters(4),stabilisationType
    REAL(DP) :: PHIMS,PHINS
    REAL(DP) :: dPhi_dX_Velocity(27,3),dPhi_dX_Pressure(27,3),DPHINS2_DXI(3,3)
    REAL(DP) :: jacobianMomentum(3),jacobianContinuity
    REAL(DP) :: DXI_DX(3,3)
    REAL(DP) :: meshVelocity(3),velocity(3),velocityPrevious(3),velocityDeriv(3,3), &
      & velocity2Deriv(3,3,3),pressure,pressureDeriv(3)
    REAL(DP) :: JGW,SUM,SUM2,SUPG,PSPG,LSIC,crossStress,reynoldsStress,momentumTerm
    REAL(DP) :: uDotGu,doubleDotG,tauSUPS,traceG,nuLSIC,timeIncrement,elementInverse,C1,stabilisationValueDP
    REAL(DP) :: tauC,tauMp,tauMu
    REAL(DP) :: residualMomentum(3),residualContinuity
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: linearElement

    ENTERS("NavierStokes_ResidualBasedStabilisation",err,error,*999)

    ! Nullify all local pointers
    NULLIFY(basisVelocity)
    NULLIFY(basisPressure)
    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(equationsSetField)
    NULLIFY(quadratureVelocity)
    NULLIFY(quadraturePressure)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(geometricField)
    NULLIFY(fieldVariable)
    NULLIFY(vectorMatrices)
    NULLIFY(nonlinearMatrices)
    NULLIFY(jacobianMatrix)
    NULLIFY(vectorEquations)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a Navier-Stokes type equations set.", &
      & err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) &
        & CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
      CALL EquationsMatricesNonlinear_JacobianMatrixGet(nonlinearMatrices,1,jacobianMatrix,err,error,*999)

      !Set general and specific pointers
      fieldVariable=>nonlinearMapping%residualVariables(1)%ptr
      fieldVariableType=fieldVariable%VARIABLE_TYPE
      numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
      meshComponent1=fieldVariable%components(1)%MESH_COMPONENT_NUMBER
      meshComponent2=fieldVariable%components(fieldVariable%NUMBER_OF_COMPONENTS)%MESH_COMPONENT_NUMBER
      basisVelocity=>dependentField%decomposition%domain(meshComponent1)%ptr%topology%elements%elements(elementNumber)%basis
      basisPressure=>dependentField%decomposition%domain(meshComponent2)%ptr%topology%elements%elements(elementNumber)%basis

      IF(basisVelocity%INTERPOLATION_ORDER(1).LE.1) THEN
        linearElement = .TRUE.
      ELSE
        ! higher order element type- can calculate 2nd order terms
        linearElement = .FALSE.
      END IF

      quadratureVelocity=>basisVelocity%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
      quadraturePressure=>basisPressure%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

      ! Stabilisation type (default 1 for RBS, 2 for RBVM, 0 for none)
      CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,stabilisationValueDP, &
        & err,error,*999)
      stabilisationType=NINT(stabilisationValueDP)
      ! Skip if type 0
      IF(stabilisationType > 0) THEN
        ! Get time step size and calc time derivative
        timeIncrement=equationsSet%deltaTime
        !CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        ! & 3,timeIncrement,err,error,*999)
        ! TODO: put this somewhere more sensible. This is a workaround since we don't have access to the dynamic solver values
        !       at this level in the element loop
        IF(equationsSet%specification(3)/=EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.AND.timeIncrement < ZERO_TOLERANCE) THEN
          CALL FlagError("Please set the equations set field time increment to a value > 0.",err,error,*999)
        END IF
        ! Stabilisation type (default 1 for RBS)
        CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,stabilisationValueDP, &
          & err,error,*999)
        stabilisationType=NINT(stabilisationValueDP)
        ! User specified or previously calculated C1
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elementNumber,10, &
          & elementInverse,err,error,*999)

        ! Get previous timestep values
        velocityPrevious=0.0_DP
        IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber,equations% &
            & interpolation%prevDependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation% &
            & prevDependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          velocityPrevious(1:numberOfDimensions)=equations%interpolation%prevDependentInterpPoint(fieldVariableType)%ptr% &
            & values(1:numberOfDimensions,NO_PART_DERIV)
        END IF

        ! Interpolate current solution velocity/pressure field values
        CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
          & dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        IF(linearElement) THEN
          ! Get 1st order derivatives for current timestep value
          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation% &
            & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        ELSE
          ! Get 2nd order derivatives for current timestep value
          CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation%&
            & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        END IF
        velocity=0.0_DP
        velocityDeriv=0.0_DP
        velocity2Deriv=0.0_DP
        pressure=0.0_DP
        pressureDeriv=0.0_DP
        DO i=1,numberOfDimensions
          velocity(i)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr%VALUES(i,NO_PART_DERIV)
          velocityDeriv(i,1)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr%VALUES(i,PART_DERIV_S1)
          velocityDeriv(i,2)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr%VALUES(i,PART_DERIV_S2)
          IF(.NOT. linearElement) THEN
            velocity2Deriv(i,1,1)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
              & VALUES(i,PART_DERIV_S1_S1)
            velocity2Deriv(i,1,2)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
              & VALUES(i,PART_DERIV_S1_S2)
            velocity2Deriv(i,2,1)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
              & VALUES(i,PART_DERIV_S1_S2)
            velocity2Deriv(i,2,2)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
              & VALUES(i,PART_DERIV_S2_S2)
          END IF
          IF(numberOfDimensions > 2) THEN
            velocityDeriv(i,3)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr%VALUES(i,PART_DERIV_S3)
            IF(.NOT. linearElement) THEN
              velocity2Deriv(i,1,3)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
                & VALUES(i,PART_DERIV_S1_S3)
              velocity2Deriv(i,2,3)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
                & VALUES(i,PART_DERIV_S2_S3)
              velocity2Deriv(i,3,1)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
                & VALUES(i,PART_DERIV_S1_S3)
              velocity2Deriv(i,3,2)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
                & VALUES(i,PART_DERIV_S2_S3)
              velocity2Deriv(i,3,3)=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr% &
                & VALUES(i,PART_DERIV_S3_S3)
            END IF
          END IF
        END DO
        pressureIndex = numberOfDimensions + 1
        pressure=equations%interpolation%dependentInterpPoint(fieldVariableType)%ptr%VALUES(pressureIndex,NO_PART_DERIV)
        pressureDeriv(1)=equations%interpolation%dependentInterpPoint(fieldVariableType)% &
          & PTR%VALUES(pressureIndex,PART_DERIV_S1)
        pressureDeriv(2)=equations%interpolation%dependentInterpPoint(fieldVariableType)% &
          & PTR%VALUES(pressureIndex,PART_DERIV_S2)
        IF(numberOfDimensions > 2) THEN
          pressureDeriv(3)=equations%interpolation%dependentInterpPoint(fieldVariableType)% &
            & PTR%VALUES(pressureIndex,PART_DERIV_S3)
        END IF
        DXI_DX=0.0_DP
        DO i=1,numberOfDimensions
          DO j=1,numberOfDimensions
            DXI_DX(j,i)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
              & DXI_DX(j,i)
          END DO
        END DO

        meshVelocity=0.0_DP
        IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,elementNumber,equations% &
            & interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation% &
            & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          meshVelocity(1:numberOfDimensions)=equations%interpolation%independentInterpPoint(fieldVariableType)%ptr% &
            & values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        ! Get number of element parameters for each dependent component
        numberOfElementParameters=basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS
        numberOfElementParameters(numberOfDimensions+1)=basisPressure%NUMBER_OF_ELEMENT_PARAMETERS
        ! Calculate dPhi/dX
        dPhi_dX_Velocity=0.0_DP
        dPhi_dX_Pressure=0.0_DP
        DO ms=1,numberOfElementParameters(1)
          DO nj=1,numberOfDimensions
            dPhi_dX_Velocity(ms,nj)=0.0_DP
            DO ni=1,numberOfDimensions
              dPhi_dX_Velocity(ms,nj)=dPhi_dX_Velocity(ms,nj) + &
                & quadratureVelocity%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),gaussNumber)* &
                & DXI_DX(ni,nj)
            END DO
          END DO
        END DO
        DO ms=1,numberOfElementParameters(numberOfDimensions+1)
          DO nj=1,numberOfDimensions
            dPhi_dX_Pressure(ms,nj)=0.0_DP
            DO ni=1,numberOfDimensions
              dPhi_dX_Pressure(ms,nj)=dPhi_dX_Pressure(ms,nj) + &
                & quadraturePressure%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),gaussNumber)* &
                & DXI_DX(ni,nj)
            END DO
          END DO
        END DO
        JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
          & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)

        !----------------------------------------------------------------------------------
        ! C a l c u l a t e   d i s c r e t e   r e s i d u a l s
        !----------------------------------------------------------------------------------
        SUM = 0.0_DP
        residualMomentum = 0.0_DP
        residualContinuity = 0.0_DP
        ! Calculate momentum residual
        DO i=1,numberOfDimensions
          SUM = 0.0_DP
          ! velocity time derivative
          IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
!!TODO: Should interpolate previous mesh velocity so that the delta velocity should be
!!      ((velocity - meshVelocity) - (previousVelocity - previousMeshVelocity)) however
!!      mesh velocity doesn't really change that much and so previousMeshVelocity ~ meshVelocity
!!      and so it will cancel out.
            SUM = rho*(velocity(i)-velocityPrevious(i))/timeIncrement
          END IF
          DO j=1,numberOfDimensions
            ! pressure gradient
            SUM = SUM + pressureDeriv(j)*DXI_DX(j,i)
            DO k=1,numberOfDimensions
              !Convective term
              SUM = SUM +rho*((velocity(j)-meshVelocity(j))*(velocityDeriv(i,k)*DXI_DX(k,j)))
              IF(.NOT. linearElement) THEN
                DO l=1,numberOfDimensions
                  ! viscous stress: only if quadratic or higher basis defined for laplacian
                  SUM = SUM - mu*(velocity2Deriv(i,k,l)*DXI_DX(k,j)*DXI_DX(l,j))
                END DO
              END IF
            END DO
          END DO
          residualMomentum(i) = SUM
        END DO
        ! Calculate continuity residual
        SUM = 0.0_DP
        DO i=1,numberOfDimensions
          DO j=1,numberOfDimensions
            SUM= SUM + velocityDeriv(i,j)*DXI_DX(j,i)
          END DO
        END DO
        residualContinuity = SUM

        ! Constant of element inverse inequality
        IF(elementInverse > -ZERO_TOLERANCE) THEN
          ! Use user-defined value if specified (default -1)
          C1 = elementInverse
        ELSE IF(linearElement) THEN
          C1=3.0_DP
        ELSE
          IF(numberOfDimensions==2 .AND. basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS==9 &
            & .AND. basisVelocity%INTERPOLATION_ORDER(1)==2) THEN
            C1=24.0_DP
          ELSE IF(numberOfDimensions==3 .AND. basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS==27 &
            & .AND. basisVelocity%INTERPOLATION_ORDER(1)==2) THEN
            C1=12.0_DP
            !TODO: Expand C1 for more element types
          ELSE
            CALL FlagError("Element inverse estimate undefined on element " &
              & //TRIM(NumberToVString(elementNumber,"*",err,error)),err,error,*999)
          END IF
        END IF
        ! Update element inverse value if calculated
        IF(ABS(C1-elementInverse) > ZERO_TOLERANCE) THEN
          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & elementNumber,10,C1,err,error,*999)
          ! Should probably move this field update it only happens one time for each element, when C1 undefined
!!TODO: CHANGE THIS. IT WILL CURRENTLY UPDATE FOR EVERY GAUSS POINT FOR EVERY ELEMENT.
          CALL Field_ParameterSetUpdateStart(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
        END IF

        !----------------------------------------------------------
        ! S t a b i l i z a t i o n    C o n s t a n t s    (Taus)
        !----------------------------------------------------------
        IF(stabilisationType == 1 .OR. stabilisationType == 2) THEN
          ! Bazilevs method for calculating tau
          pointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          uDotGu = 0.0_DP
          DO i=1,numberOfDimensions
            DO j=1,numberOfDimensions
              uDotGu = uDotGu + (velocity(i)-meshVelocity(i))*pointMetrics%GU(i,j)*(velocity(j)-meshVelocity(j))
            END DO
          END DO
          doubleDotG = 0.0_DP
          DO i=1,numberOfDimensions
            DO j=1,numberOfDimensions
              doubleDotG = doubleDotG + pointMetrics%GU(i,j)*pointMetrics%GU(i,j)
            END DO
          END DO
          ! Calculate tauSUPS (used for both PSPG and SUPG weights)
          IF(equationsSet%specification(3) == EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
            tauSUPS = (uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
          ELSE
            tauSUPS = ((4.0_DP/(timeIncrement**2.0_DP)) + uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
          END IF

          ! Calculate nu_LSIC (Least-squares incompressibility constraint)
          CALL Trace(pointMetrics%GU(1:numberOfDimensions,1:numberOfDimensions),traceG,err,error,*999)
          nuLSIC = 1.0_DP/(tauSUPS*traceG)

          tauMp = tauSUPS
          tauMu = tauSUPS
          tauC = nuLSIC

        ELSE
          CALL FlagError("A tau factor has not been defined for the stabilisation type of " &
            & //TRIM(NumberToVString(stabilisationType,"*",err,error)),err,error,*999)
        END IF

        !-------------------------------------------------------------------------------------------------
        ! A d d   s t a b i l i z a t i o n   f a c t o r s   t o   e l e m e n t   m a t r i c e s
        !-------------------------------------------------------------------------------------------------
        jacobianMomentum = 0.0_DP
        jacobianContinuity = 0.0_DP
        mhs = 0
        DO mh=1,numberOfDimensions+1
          DO ms=1,numberOfElementParameters(mh)
            mhs = mhs + 1
            IF(mh <= numberOfDimensions) THEN
              PHIMS=quadratureVelocity%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)
              JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)
            ELSE
              PHIMS=quadraturePressure%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)
              JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & quadraturePressure%GAUSS_WEIGHTS(gaussNumber)
            END IF
            !------------------
            ! J A C O B I A N
            !------------------
            IF(jacobianFlag) THEN
              nhs = 0
              DO nh=1,numberOfDimensions+1
                DO ns=1,numberOfElementParameters(nh)
                  nhs=nhs+1
                  ! Note that we still need to assemble the vector momentum jacobian for PSPG in the continuity row
                  IF(nh <= numberOfDimensions) THEN
                    PHINS=quadratureVelocity%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,gaussNumber)
                  ELSE
                    PHINS=quadraturePressure%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,gaussNumber)
                  END IF

                  ! Calculate jacobians of the discrete residual terms
                  jacobianMomentum = 0.0_DP
                  IF(nh == numberOfDimensions+1) THEN
                    ! d(Momentum(mh))/d(Pressure)
                    DO i=1,numberOfDimensions
                      jacobianMomentum(i) = dPhi_dX_Pressure(ns,i)
                    END DO
                    jacobianContinuity=0.0_DP
                  ELSE
                    DPHINS2_DXI=0.0_DP
                    IF(.NOT. linearElement) THEN
                      DPHINS2_DXI(1,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,gaussNumber)
                      DPHINS2_DXI(1,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,gaussNumber)
                      DPHINS2_DXI(2,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,gaussNumber)
                      DPHINS2_DXI(2,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,gaussNumber)
                      IF(numberOfDimensions > 2) THEN
                        DPHINS2_DXI(1,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,gaussNumber)
                        DPHINS2_DXI(2,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,gaussNumber)
                        DPHINS2_DXI(3,1)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,gaussNumber)
                        DPHINS2_DXI(3,2)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,gaussNumber)
                        DPHINS2_DXI(3,3)=quadratureVelocity%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,gaussNumber)
                      END IF
                    END IF
                    ! d(Momentum)/d(Velocity(nh))
                    jacobianMomentum = 0.0_DP
                    DO i=1,numberOfDimensions
                      SUM = 0.0_DP
                      !Note: Convective term split using product rule
                      !Convective term 1: applies to all velocity components
                      DO j=1,numberOfDimensions
                        SUM = SUM + rho*PHINS*velocityDeriv(i,j)*DXI_DX(j,nh)
                      END DO
                      !Diagonal terms
                      IF(i==nh) THEN
                        !Transient
                        IF(equationsSet%specification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                          SUM = SUM + rho*PHINS/timeIncrement
                        END IF
                        !Convective 2: nh component only
                        DO j=1,numberOfDimensions
                          SUM = SUM + rho*(velocity(j)-meshVelocity(j))*dPhi_dX_Velocity(ns,j)
                        END DO
                        IF(.NOT. linearElement) THEN
                          !Viscous laplacian term
                          DO j=1,numberOfDimensions
                            DO k=1,numberOfDimensions
                              DO l=1,numberOfDimensions
                                SUM=SUM-mu*DPHINS2_DXI(k,l)*DXI_DX(k,j)*DXI_DX(l,j)
                              END DO
                            END DO
                          END DO
                        END IF
                      END IF
                      jacobianMomentum(i)=SUM
                    END DO
                    ! Continuity/velocity
                    jacobianContinuity = dPhi_dX_Velocity(ns,nh)
                  END IF
                  ! Calculate jacobian of discrete residual * RBS factors (apply product rule if neccesary)

                  ! PSPG: Pressure stabilising Petrov-Galerkin
                  IF(mh == numberOfDimensions+1) THEN
                    PSPG = 0.0_DP
                    SUM = 0.0_DP
                    DO i=1,numberOfDimensions
                      SUM = SUM + dPhi_dX_Pressure(ms,i)*jacobianMomentum(i)
                    END DO
                    PSPG = tauMp*SUM/rho*JGW

                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+PSPG

                    ! SUPG: Streamline upwind/Petrov-Galerkin
                    ! LSIC: Least-squares incompressibility constraint
                  ELSE
                    SUPG=0.0_DP
                    LSIC=0.0_DP

                    SUM=0.0_DP
                    IF(nh <= numberOfDimensions) THEN
                      SUPG= SUPG + PHINS*dPhi_dX_Velocity(ms,nh)*residualMomentum(mh)
                    END IF
                    DO i=1,numberOfDimensions
                      SUM = SUM + (velocity(i)-meshVelocity(i))*dPhi_dX_Velocity(ms,i)
                    END DO
                    SUPG = tauMu*(SUPG + SUM*jacobianMomentum(mh))

                    SUM=0.0_DP
                    DO i=1,numberOfDimensions
                      SUM = SUM + dPhi_dX_Velocity(ms,i)
                    END DO
                    LSIC = tauC*rho*dPhi_dX_Velocity(ms,mh)*jacobianContinuity

                    momentumTerm = (SUPG + LSIC)*JGW

                    IF(stabilisationType == 2) THEN
                      ! Additional terms for RBVM
                      crossStress=0.0_DP
                      reynoldsStress=0.0_DP
                      crossStress = 0.0_DP
                      IF(nh <= numberOfDimensions) THEN
                        IF(mh == nh) THEN
                          DO i=1,numberOfDimensions
                            crossStress= crossStress + dPhi_dX_Velocity(ns,i)*residualMomentum(i)
                          END DO
                        END IF
                      END IF
                      SUM2=0.0_DP
                      DO i=1,numberOfDimensions
                        SUM=0.0_DP
                        ! dU_mh/dX_i
                        DO j=1,numberOfDimensions
                          SUM= SUM + velocityDeriv(mh,j)*DXI_DX(j,i)
                        END DO
                        ! Jm_i*dU_mh/dX_i
                        SUM2 = SUM2 + jacobianMomentum(i)*SUM
                      END DO
                      crossStress = -tauMu*(crossStress + SUM2)

                      reynoldsStress = 0.0_DP
                      SUM = 0.0_DP
                      !Rm_mh.Rm_i.dPhi/dX_i
                      DO i=1,numberOfDimensions
                        SUM = SUM + jacobianMomentum(mh)*residualMomentum(i)*dPhi_DX_Velocity(ms,i)
                        SUM = SUM + jacobianMomentum(i)*residualMomentum(mh)*dPhi_DX_Velocity(ms,i)
                      END DO
                      reynoldsStress = -tauMu*tauMu*SUM

                      momentumTerm = momentumTerm + (crossStress + reynoldsStress)*JGW
                    END IF

                    ! Add stabilisation to element jacobian
                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)= &
                      & jacobianMatrix%elementJacobian%matrix(mhs,nhs)+momentumTerm

                  END IF
                END DO
              END DO

              !-----------------
              ! R E S I D U A L
              !-----------------
            ELSE
              ! PSPG: Pressure stabilising Petrov-Galerkin
              IF(mh == numberOfDimensions+1) THEN
                SUM = 0.0_DP
                DO i=1,numberOfDimensions
                  SUM = SUM + dPhi_dX_Pressure(ms,i)*residualMomentum(i)
                END DO
                PSPG = SUM*(tauMp/rho)*JGW
                nonlinearMatrices%elementResidual%vector(mhs)= &
                  & nonlinearMatrices%elementResidual%vector(mhs) + PSPG

                ! SUPG: Streamline upwind/Petrov-Galerkin
                ! LSIC: Least-squares incompressibility constraint
              ELSE
                SUPG=0.0_DP
                LSIC=0.0_DP

                ! u_i*Rm_mh*dv_mh/dx_i
                SUM=0.0_DP
                DO i=1,numberOfDimensions
                  SUM = SUM + (velocity(i)-meshVelocity(i))*dPhi_dX_Velocity(ms,i)
                END DO
                SUPG = tauMu*SUM*residualMomentum(mh)

                LSIC = tauC*rho*dPhi_dX_Velocity(ms,mh)*residualContinuity
                momentumTerm = (SUPG + LSIC)*JGW

                IF(stabilisationType ==2) THEN
                  ! Additional terms for RBVM
                  crossStress=0.0_DP
                  reynoldsStress=0.0_DP
                  SUM2 = 0.0_DP
                  DO i=1,numberOfDimensions
                    SUM = 0.0_DP
                    ! dU_mh/dX_i
                    DO j=1,numberOfDimensions
                      SUM= SUM + velocityDeriv(mh,j)*DXI_DX(j,i)
                    END DO
                    ! Rm_i.dU_mh/dX_i
                    SUM2= SUM2 + residualMomentum(i)*SUM
                  END DO
                  crossStress= -tauMu*PHIMS*SUM2

                  reynoldsStress = 0.0_DP
                  SUM = 0.0_DP
                  !Rm_mh.Rm_i.dPhi/dX_i
                  DO i=1,numberOfDimensions
                    SUM = SUM + dPhi_dX_Velocity(ms,i)*residualMomentum(i)*residualMomentum(mh)
                  END DO
                  reynoldsStress = -SUM*(tauMu*tauMu)/rho
                  momentumTerm = momentumTerm + (crossStress + reynoldsStress)*JGW
                END IF

                ! Add stabilisation to element residual
                nonlinearMatrices%elementResidual%vector(mhs)= &
                  & nonlinearMatrices%elementResidual%vector(mhs) + momentumTerm
              END IF
            END IF ! jacobian/residual
          END DO !ms
        END DO !mh

      END IF ! check stabilisation type
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
        & " is not a valid subtype to use SUPG weighting functions."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_ResidualBasedStabilisation")
    RETURN
999 ERRORSEXITS("NavierStokes_ResidualBasedStabilisation",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_ResidualBasedStabilisation

  !
  !================================================================================================================================
  !

  !>Calculate element-level scale factors: CFL, cell Reynolds number
  SUBROUTINE NavierStokes_CalculateElementMetrics(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: basisVelocity
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: equationsSetField
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureVelocity
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable

    INTEGER(INTG) :: fieldVariableType,meshComponent1
    INTEGER(INTG) :: numberOfDimensions,mh
    INTEGER(INTG) :: gaussNumber
    INTEGER(INTG) :: i,j,ms
    INTEGER(INTG) :: numberOfElementParameters
    INTEGER(INTG) :: LWORK,INFO
    REAL(DP) :: cellReynoldsNumber,cellCourantNumber,timeIncrement
    REAL(DP) :: dPhi_dX_Velocity(27,3)
    REAL(DP) :: DXI_DX(3,3)
    REAL(DP) :: velocity(3),meshVelocity(3),avgVelocity(3),velocityNorm,velocityPrevious(3),velocityDeriv(3,3)
    REAL(DP) :: PHIMS,JGW,SUM,SUM2,mu,rho,normCMatrix,normKMatrix,normMMatrix,muScale
    REAL(DP) :: CMatrix(27,3),KMatrix(27,3),MMatrix(27,3)
    REAL(DP) :: svd(3),U(27,27),VT(3,3)
    REAL(DP), ALLOCATABLE :: WORK(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_CalculateElementMetrics",err,error,*999)

    ! Nullify all local pointers
    NULLIFY(basisVelocity)
    NULLIFY(equations)
    NULLIFY(vectorEquations)
    NULLIFY(vectorMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(equationsSetField)
    NULLIFY(quadratureVelocity)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(independentField)
    NULLIFY(fieldVariable)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      ENDIF
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      !Set general and specific pointers
      fieldVariable=>nonlinearMapping%residualVariables(1)%ptr
      fieldVariableType=fieldVariable%VARIABLE_TYPE
      numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
      meshComponent1=fieldVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
      basisVelocity=>dependentField%DECOMPOSITION%DOMAIN(meshComponent1)%ptr% &
        & TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)%BASIS
      quadratureVelocity=>basisVelocity%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
      numberOfElementParameters=basisVelocity%NUMBER_OF_ELEMENT_PARAMETERS

      ! Get time step size
      timeIncrement=equationsSet%deltaTime
      !CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      ! & 3,timeIncrement,err,error,*999)

      ! Loop over gauss points
      CMatrix = 0.0_DP
      MMatrix = 0.0_DP
      KMatrix = 0.0_DP
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%INTERPOLATION% &
        & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      END IF
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)

      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
        & equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
        & equations%interpolation%prevDependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
          & equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
      ENDIF

      avgVelocity = 0.0_DP
      DO gaussNumber = 1,quadratureVelocity%NUMBER_OF_GAUSS

        ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
        IF(equationsSet%specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
          ! Note the constant from the U_VARIABLE is a scale factor
          muScale = mu
          ! Get the gauss point based value returned from the CellML solver
          CALL Field_ParameterSetGetLocalGaussPoint(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,gaussNumber,elementNumber,1,mu,err,error,*999)
          mu=mu*muScale
        END IF

        ! Get previous timestep values
        velocityPrevious=0.0_DP
        CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation%&
          & prevDependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        velocityPrevious(1:numberOfDimensions)=equations%interpolation%prevDependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
          & values(1:numberOfDimensions,NO_PART_DERIV)

        ! Interpolate current solution velocity and first deriv field values
        ! Get 1st order derivatives for current timestep value
        CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber, &
          & equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        velocity=0.0_DP
        velocityDeriv=0.0_DP
        DO i=1,numberOfDimensions
          velocity(i)=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(i,NO_PART_DERIV)
          velocityDeriv(i,1)=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)% &
            & PTR%VALUES(i,PART_DERIV_S1)
          velocityDeriv(i,2)=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)% &
            & PTR%VALUES(i,PART_DERIV_S2)
          IF(numberOfDimensions > 2) THEN
            velocityDeriv(i,3)=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)% &
              & PTR%VALUES(i,PART_DERIV_S3)
          END IF
        END DO

        meshVelocity=0.0_DP
        IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber, &
            & equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          meshVelocity(1:numberOfDimensions)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
            & values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        ! get dXi/dX deriv
        CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,equations%interpolation%&
          & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        DXI_DX=0.0_DP
        DXI_DX(1:numberOfDimensions,1:numberOfDimensions)= &
          & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
          & DXI_DX(1:numberOfDimensions,1:numberOfDimensions)

        ! Calculate dPhi/dX
        dPhi_dX_Velocity=0.0_DP
        DO ms=1,numberOfElementParameters
          DO i=1,numberOfDimensions
            dPhi_dX_Velocity(ms,i)=0.0_DP
            DO j=1,numberOfDimensions
              dPhi_dX_Velocity(ms,i)=dPhi_dX_Velocity(ms,i) + &
                & quadratureVelocity%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(j),gaussNumber)* &
                & DXI_DX(j,i)
            END DO
          END DO
        END DO

        JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
          & quadratureVelocity%GAUSS_WEIGHTS(gaussNumber)
        DO mh=1,numberOfDimensions
          DO ms=1,numberOfElementParameters
            PHIMS=quadratureVelocity%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,gaussNumber)

            ! c_(a,i)
            SUM=0.0_DP
            DO i=1,numberOfDimensions
              DO j=1,numberOfDimensions
                SUM = SUM + (velocity(i)-meshVelocity(i))*velocityDeriv(mh,j)*DXI_DX(j,i)
              END DO
            END DO
            CMatrix(ms,mh)=CMatrix(ms,mh) + rho*PHIMS*SUM*JGW

            ! ~k_(a,i)
            SUM=0.0_DP
            DO i=1,numberOfDimensions
              SUM = SUM + (velocity(i)-meshVelocity(i))*dPhi_dX_Velocity(ms,i)
            END DO
            SUM2=0.0_DP
            DO i=1,numberOfDimensions
              DO j=1,numberOfDimensions
                SUM2 = SUM2 + (velocity(i)-meshVelocity(i))*velocityDeriv(mh,j)*DXI_DX(j,i)
              END DO
            END DO
            KMatrix(ms,mh)=KMatrix(ms,mh)+rho*SUM*SUM2*JGW

            ! m_(a,i)
!!TODO: Should interpolate previous mesh velocity so that the delta velocity should be
!!      ((velocity - meshVelocity) - (previousVelocity - previousMeshVelocity)) however
!!      mesh velocity doesn't really change that much and so previousMeshVelocity ~ meshVelocity
!!      and so it will cancel out.
            MMatrix(ms,mh)=MMatrix(ms,mh)+rho*PHIMS*(velocity(mh)-velocityPrevious(mh))/timeIncrement*JGW

          END DO !ms
        END DO !mh

        avgVelocity= avgVelocity + (velocity-meshVelocity)/quadratureVelocity%NUMBER_OF_GAUSS
      END DO ! gauss loop

      LWORK=MAX(1,3*MIN(numberOfElementParameters,numberOfDimensions)+ &
        & MAX(numberOfElementParameters,numberOfDimensions),5*MIN(numberOfElementParameters,numberOfDimensions))
      ALLOCATE(WORK(LWORK))

      ! compute the singular value decomposition (SVD) using LAPACK
      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,CMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
      normCMatrix=svd(1)
      IF(INFO /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END IF

      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,KMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
      normKMatrix=svd(1)
      IF(INFO /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END IF

      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,MMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,WORK,LWORK,INFO)
      normMMatrix=svd(1)
      IF(INFO /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END IF
      DEALLOCATE(WORK)

      CALL L2Norm(avgVelocity,velocityNorm,err,error,*999)
      cellReynoldsNumber = 0.0_DP
      cellCourantNumber = 0.0_DP
      IF(velocityNorm > ZERO_TOLERANCE) THEN
        IF(normKMatrix > ZERO_TOLERANCE) THEN
          cellReynoldsNumber = velocityNorm**2.0_DP/(mu/rho)*normCMatrix/normKMatrix
        END IF
        IF(normMMatrix > ZERO_TOLERANCE) THEN
          cellCourantNumber = timeIncrement/2.0_DP*normCMatrix/normMMatrix
        END IF
      END IF
      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & elementNumber,2,velocityNorm,err,error,*999)
      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & elementNumber,3,cellCourantNumber,err,error,*999)
      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_ELEMENT(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & elementNumber,4,cellReynoldsNumber,err,error,*999)

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
        & " is not a valid subtype to use SUPG weighting functions."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_CalculateElementMetrics")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateElementMetrics",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_CalculateElementMetrics

  !
  !================================================================================================================================
  !

  !>Calculates the boundary integration term of the finite element formulation for Navier-Stokes equation,
  !>required for pressure and multidomain boundary conditions. Optionally also includes a boundary stabilisation term.
  SUBROUTINE NavierStokes_FiniteElementBoundaryIntegrate(equationsSet,elementNumber,dependentVariable,jacobianFlag,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<The equations set to calculate the RHS term for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate the RHS term for
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    LOGICAL, INTENT(IN) ::  jacobianFlag !<Flag indicating whether this was called from the jacobian or residual evaluation routine
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(BASIS_TYPE), POINTER :: basis1,basis2,dependentBasis1,dependentBasis2
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(DECOMPOSITION_LINE_TYPE), POINTER :: line
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: geometricField,equationsSetField,dependentField,independentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters,geometricInterpolationParameters, &
      & pressureInterpolationParameters,independentInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint,geometricInterpolatedPoint, &
      & pressureInterpolatedPoint,independentInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme1,quadratureScheme2
    INTEGER(INTG) :: boundaryIdx, boundaryNumber, xiDirection(4), orientation
    INTEGER(INTG) :: componentIdx, componentIdx2, gaussIdx
    INTEGER(INTG) :: elementBaseDofIdx, nodeIdx, elementNodeIdx
    INTEGER(INTG) :: nodeDerivativeIdx,meshComponentNumber1,globalNodeDerivativeIdx,elementParameterIdx
    INTEGER(INTG) :: parameterIdx,elementDof
    INTEGER(INTG) :: parameterIdx2,elementDof2,elementBaseDofIdx2,nodeIdx2,elementNodeIdx2,nodeDerivativeIdx2
    INTEGER(INTG) :: meshComponentNumber2,globalNodeDerivativeIdx2,elementParameterIdx2
    INTEGER(INTG) :: numberOfDimensions,numberOfElementBoundaries,boundaryType,ni
    REAL(DP) :: pressure,density,jacobianGaussWeights,beta,normalFlow
    REAL(DP) :: meshVelocity(3),velocity(3),normalProjection(3),unitNormal(3),stabilisationTerm
    REAL(DP) :: boundaryNormal(3)
    REAL(DP) :: boundaryValue,normalDifference,normalTolerance
    REAL(DP) :: phim,phin
    REAL(DP) :: dPhinDXi(3)
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: integratedBoundary

    ENTERS("NavierStokes_FiniteElementBoundaryIntegrate",err,error,*999)

    NULLIFY(decomposition)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis1,dependentBasis2)
    NULLIFY(equations)
    NULLIFY(equationsSetField)
    NULLIFY(vectorMatrices)
    NULLIFY(face,line)
    NULLIFY(basis1,basis2)
    NULLIFY(quadratureScheme1,quadratureScheme2)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(pressureInterpolatedPoint)
    NULLIFY(pressureInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(nonlinearMatrices)
    NULLIFY(dependentField)
    NULLIFY(geometricField)
    NULLIFY(independentField)
    NULLIFY(vectorEquations)
    NULLIFY(vectorMapping)
    NULLIFY(nonlinearMapping)
    NULLIFY(vectorMatrices)

    ! Get pointers and perform sanity checks
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)

      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
      IF(jacobianFlag) THEN
        jacobianMatrix=>nonlinearMatrices%JACOBIANS(1)%PTR
      ENDIF
      IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
        & equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      ENDIF

      ! Check whether this element contains an integrated boundary type
      CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        & elementNumber,9,boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      integratedBoundary = .FALSE.
      IF(boundaryType == BOUNDARY_CONDITION_PRESSURE) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_FIXED_PRESSURE) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_COUPLING_STRESS) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_FIXED_CELLML) integratedBoundary = .TRUE.

      !Get the mesh decomposition and basis
      numberOfDimensions = dependentVariable%NUMBER_OF_COMPONENTS - 1
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      meshComponentNumber1=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
      meshComponentNumber2=dependentVariable%COMPONENTS(numberOfDimensions+1)%MESH_COMPONENT_NUMBER
      dependentBasis1=>decomposition%DOMAIN(meshComponentNumber1)%PTR%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(elementNumber)%BASIS
      dependentBasis2=>decomposition%DOMAIN(meshComponentNumber2)%PTR%TOPOLOGY%ELEMENTS% &
        & ELEMENTS(elementNumber)%BASIS
      decompElement=>DECOMPOSITION%TOPOLOGY%ELEMENTS%ELEMENTS(elementNumber)

      IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
        & equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        independentInterpolationParameters=>equations%interpolation%independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        independentInterpolatedPoint=>equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
      ENDIF

      !Determine if this is a 2D or 3D problem with line/face parameters calculated
      IF(integratedBoundary) THEN
        IF(numberOfDimensions /= 2 .AND. numberOfDimensions /=3) THEN
          localError="Invalid number of dimensions ("//TRIM(NUMBER_TO_VSTRING(numberOfDimensions,"*",ERR,ERROR))// &
            & ") for a 2D or 3D Navier-Stokes problem."
          CALL FlagError(localError,ERR,ERROR,*999)
        END IF
        IF(numberOfDimensions == 3 .AND. decomposition%CALCULATE_FACES) THEN
          numberOfElementBoundaries = dependentBasis1%NUMBER_OF_LOCAL_FACES
        ELSE IF(numberOfDimensions == 2 .AND. decomposition%CALCULATE_LINES) THEN
          numberOfElementBoundaries = dependentBasis1%NUMBER_OF_LOCAL_LINES
        ELSE
          integratedBoundary = .FALSE.
        END IF
      END IF

      !Only apply to required element boundaries
      IF(integratedBoundary) THEN
        !Get the density
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,&
         & FIELD_VALUES_SET_TYPE,2,density,err,error,*999)
        ! Get the boundary element parameters
        CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 1,beta,err,error,*999)
        boundaryNormal = 0.0_DP
        DO componentIdx=1,3
          CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & elementNumber,componentIdx+4,boundaryNormal(componentIdx),err,error,*999)
        END DO

        ! Loop over the boundaries (lines or faces) for this element
        DO boundaryIdx=1,numberOfElementBoundaries
          geometricInterpolationParameters=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
          dependentInterpolationParameters=>equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr

          ! Get 3D face specific parameters
          IF(numberOfDimensions == 3) THEN
            IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
              boundaryNumber=decompElement%ELEMENT_FACES(boundaryIdx)
            ELSE
              CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
            END IF
            face=>decomposition%TOPOLOGY%FACES%FACES(boundaryNumber)
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(face%BOUNDARY_FACE)) CYCLE
            SELECT CASE(dependentBasis1%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              xiDirection(3)=ABS(face%XI_DIRECTION)
              xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
              xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
              orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*face%XI_DIRECTION)
            CASE(BASIS_SIMPLEX_TYPE)
              orientation=1
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis1%TYPE,"*",ERR,ERROR))// &
                & " is not yet implemented for Navier-Stokes boundary integration."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
            basis1=>decomposition%DOMAIN(meshComponentNumber1)%PTR%TOPOLOGY%FACES%FACES(boundaryNumber)%BASIS
            basis2=>decomposition%DOMAIN(meshComponentNumber2)%PTR%TOPOLOGY%FACES%FACES(boundaryNumber)%BASIS
            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
              & geometricInterpolationParameters,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
              & dependentInterpolationParameters,err,error,*999)
            IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
                & independentInterpolationParameters,err,error,*999)
            ENDIF

          ! Get 2D line specific parameters
          ELSE IF(numberOfDimensions == 2) THEN
            IF(ALLOCATED(decompElement%ELEMENT_LINES)) THEN
              boundaryNumber=decompElement%ELEMENT_LINES(boundaryIdx)
            ELSE
              CALL FlagError("Decomposition element lines is not allocated.",err,error,*999)
            END IF
            line=>decomposition%TOPOLOGY%LINES%LINES(boundaryNumber)
            !This speeds things up but is also important, as non-boundary lines have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(line%BOUNDARY_LINE)) CYCLE
            SELECT CASE(dependentBasis1%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              xiDirection(2)=ABS(line%XI_DIRECTION)
              xiDirection(1)=OTHER_XI_DIRECTIONS2(xiDirection(2))
              orientation=SIGN(1,OTHER_XI_ORIENTATIONS2(xiDirection(1))*line%XI_DIRECTION)
            CASE(BASIS_SIMPLEX_TYPE)
              orientation=1
            CASE DEFAULT
              localError="Line integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis1%TYPE,"*",ERR,ERROR))// &
                & " is not yet implemented for Navier-Stokes boundary integration."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
            basis1=>decomposition%DOMAIN(meshComponentNumber1)%PTR%TOPOLOGY%LINES%LINES(boundaryNumber)%BASIS
            basis2=>decomposition%DOMAIN(meshComponentNumber2)%PTR%TOPOLOGY%LINES%LINES(boundaryNumber)%BASIS
            CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
              & geometricInterpolationParameters,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
              & dependentInterpolationParameters,err,error,*999)
            IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
                & independentInterpolationParameters,err,error,*999)
            ENDIF
          END IF

          quadratureScheme1=>basis1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          quadratureScheme2=>basis2%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(dependentVariable%VARIABLE_TYPE)%ptr
          pointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          ! Loop over gauss points
          DO gaussIdx=1,quadratureScheme1%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & geometricInterpolatedPoint,err,error,*999)
            normalTolerance=0.1_DP
            unitNormal = 0.0_DP
            velocity = 0.0_DP
            meshVelocity = 0.0_DP
            pressure = 0.0_DP
            IF(numberOfDimensions == 3) THEN
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE,pointMetrics,err,error,*999)
              ! Make sure this is the boundary that corresponds with the provided normal (could be a wall rather than inlet/outlet)
              CALL CrossProduct(pointMetrics%DX_DXI(:,1),pointMetrics%DX_DXI(:,2),normalProjection,err,error,*999)
              normalProjection = normalProjection*orientation
              CALL Normalise(normalProjection,unitNormal,err,error,*999)
              CALL L2Norm(boundaryNormal-unitNormal,normalDifference,err,error,*999)
              IF(normalDifference>normalTolerance) EXIT
              CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
                & dependentInterpolationParameters,err,error,*999)
            ELSE
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_LINE_TYPE,pointMetrics,err,error,*999)
              ! Make sure this is the boundary that corresponds with the provided normal (could be a wall rather than inlet/outlet)
              normalProjection = [pointMetrics%DX_DXI(2,1),pointMetrics%DX_DXI(1,1),0.0_DP]*orientation
              CALL Normalise(normalProjection,unitNormal,err,error,*999)
              CALL L2Norm(boundaryNormal-unitNormal,normalDifference,err,error,*999)
              IF(normalDifference>normalTolerance) EXIT
              CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_VALUES_SET_TYPE,boundaryNumber, &
                & dependentInterpolationParameters,err,error,*999)
            END IF

            !Get interpolated velocity and pressure
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & dependentInterpolatedPoint,ERR,ERROR,*999)
            velocity(1:numberOfDimensions)=dependentInterpolatedPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            IF(equationsSet%specification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & equationsSet%specification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              !Get interpolated mesh velocity
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                & independentInterpolatedPoint,ERR,ERROR,*999)
              meshVelocity(1:numberOfDimensions)=independentInterpolatedPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            ENDIF

            ! Stabilisation term to correct for possible retrograde flow divergence.
            ! See: Moghadam et al 2011 “A comparison of outlet boundary treatments for prevention of backflow divergence..." and
            !      Ismail et al 2014 "A stable approach for coupling multidimensional cardiovascular and pulmonary networks..."
            ! Note: beta is a relative scaling factor 0 <= beta <= 1; default 0.0
            stabilisationTerm = 0.0_DP
            normalFlow = DOT_PRODUCT(velocity-meshVelocity,unitNormal)
            IF(normalFlow < -ZERO_TOLERANCE) THEN
              stabilisationTerm = normalFlow - ABS(normalFlow)
            ELSE
              stabilisationTerm = 0.0_DP
            END IF

            ! Check for Neumann integrated boundary types rather than fixed pressure types
            IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
              & boundaryType==BOUNDARY_CONDITION_FIXED_CELLML .OR. &
              & boundaryType==BOUNDARY_CONDITION_PRESSURE) THEN
              !Get the pressure value interpolation parameters
              pressureInterpolationParameters=>equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
              pressureInterpolatedPoint=>equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              IF(numberOfDimensions==3) THEN
                CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,boundaryNumber, &
                  & pressureInterpolationParameters,err,error,*999)
              ELSE
                CALL FIELD_INTERPOLATION_PARAMETERS_LINE_GET(FIELD_PRESSURE_VALUES_SET_TYPE,boundaryNumber, &
                  & pressureInterpolationParameters,err,error,*999)
              ENDIF
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                & pressureInterpolatedPoint,ERR,ERROR,*999)
              pressure=pressureInterpolatedPoint%VALUES(numberOfDimensions+1,NO_PART_DERIV)
            END IF

            !Jacobian and Gauss weighting term
            jacobianGaussWeights=pointMetrics%JACOBIAN*quadratureScheme1%GAUSS_WEIGHTS(gaussIdx)
            !Loop over field components
            DO componentIdx=1,dependentVariable%NUMBER_OF_COMPONENTS-1
              !Work out the first index of the vector for this element - (i.e. the number of previous)
              elementBaseDofIdx=dependentBasis1%NUMBER_OF_ELEMENT_PARAMETERS*(componentIdx-1)
              DO nodeIdx=1,basis1%NUMBER_OF_NODES
                IF(numberOfDimensions == 3) THEN
                  elementNodeIdx=dependentBasis1%NODE_NUMBERS_IN_LOCAL_FACE(nodeIdx,boundaryIdx)
                ELSE
                  elementNodeIdx=dependentBasis1%NODE_NUMBERS_IN_LOCAL_LINE(nodeIdx,boundaryIdx)
                END IF
                DO nodeDerivativeIdx=1,basis1%NUMBER_OF_DERIVATIVES(nodeIdx)
                  IF(numberOfDimensions == 3) THEN
                    globalNodeDerivativeIdx=dependentBasis1%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(nodeDerivativeIdx,nodeIdx,boundaryIdx)
                  ELSE
                    globalNodeDerivativeIdx=dependentBasis1%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nodeIdx,boundaryIdx)
                  END IF
                  elementParameterIdx=dependentBasis1%ELEMENT_PARAMETER_INDEX(globalNodeDerivativeIdx,elementNodeIdx)
                  parameterIdx=basis1%ELEMENT_PARAMETER_INDEX(nodeDerivativeIdx,nodeIdx)
                  elementDof=elementBaseDofIdx+elementParameterIdx
                  phim = quadratureScheme1%GAUSS_BASIS_FNS(parameterIdx,NO_PART_DERIV,gaussIdx)

                  IF (.NOT. jacobianFlag) THEN
                    IF (boundaryType==BOUNDARY_CONDITION_PRESSURE .OR. &
                      & boundaryType==BOUNDARY_CONDITION_FIXED_CELLML .OR. &
                      & boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
                      ! Integrated boundary pressure term
                      nonlinearMatrices%elementResidual%VECTOR(elementDof)=nonlinearMatrices%elementResidual%vector(elementDof)-&
                        &  pressure*unitNormal(componentIdx)*phim*jacobianGaussWeights
                    END IF
                    ! Boundary stabilisation term (if necessary )
                    IF (ABS(beta) > ZERO_TOLERANCE) THEN
                      nonlinearMatrices%elementResidual%VECTOR(elementDof)=nonlinearMatrices%elementResidual%vector(elementDof)-&
                        & 0.5_DP*beta*density*phim*(velocity(componentIdx)-meshVelocity(componentIdx))*stabilisationTerm* &
                        & jacobianGaussWeights
                    END IF
                  ! Jacobian matrix term is the derivative of the nonlinear stabilisation term
                  !Loop over field components
                  ELSE
                    IF (ABS(beta) > ZERO_TOLERANCE) THEN
                      DO componentIdx2=1,dependentVariable%NUMBER_OF_COMPONENTS-1
                        !Work out the first index of the rhs vector for this element - (i.e. the number of previous)
                        elementBaseDofIdx2=dependentBasis2%NUMBER_OF_ELEMENT_PARAMETERS*(componentIdx2-1)
                        DO nodeIdx2=1,basis2%NUMBER_OF_NODES
                          IF(numberOfDimensions == 3) THEN
                            elementNodeIdx2=dependentBasis2%NODE_NUMBERS_IN_LOCAL_FACE(nodeIdx2,boundaryIdx)
                          ELSE
                            elementNodeIdx2=dependentBasis2%NODE_NUMBERS_IN_LOCAL_LINE(nodeIdx2,boundaryIdx)
                          END IF
                          DO nodeDerivativeIdx2=1,basis2%NUMBER_OF_DERIVATIVES(nodeIdx2)
                            IF(numberOfDimensions == 3) THEN
                              globalNodeDerivativeIdx2=dependentBasis2%DERIVATIVE_NUMBERS_IN_LOCAL_FACE(nodeDerivativeIdx2, &
                                & nodeIdx2,boundaryIdx)
                            ELSE
                              globalNodeDerivativeIdx2=dependentBasis2%DERIVATIVE_NUMBERS_IN_LOCAL_LINE(nodeIdx2,boundaryIdx)
                            END IF
                            elementParameterIdx2=dependentBasis2%ELEMENT_PARAMETER_INDEX(globalNodeDerivativeIdx2,elementNodeIdx2)
                            parameterIdx2=basis2%ELEMENT_PARAMETER_INDEX(nodeDerivativeIdx2,nodeIdx2)
                            elementDof2=elementBaseDofIdx2+elementParameterIdx2
                            phin = quadratureScheme2%GAUSS_BASIS_FNS(parameterIdx2,NO_PART_DERIV,gaussIdx)
                            DO ni=1,basis2%NUMBER_OF_XI
                              dPhinDXi(ni) = quadratureScheme2%GAUSS_BASIS_FNS(parameterIdx2, &
                                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),gaussIdx)
                            END DO
                            ! Jacobian term
                            IF (componentIdx == componentIdx2) THEN
                              ! note that (u_j.n_j - |u_j.n_j|) term derivative will be zero
                              jacobianMatrix%elementJacobian%matrix(elementDof,elementDof2)= &
                                & jacobianMatrix%elementJacobian%matrix(elementDof,elementDof2) - &
                                & 0.5_DP*beta*density*phim*phin*stabilisationTerm*jacobianGaussWeights
                            END IF
                          END DO
                        END DO
                      END DO
                    END IF
                  END IF

                END DO !nodeDerivativeIdx
              END DO !nodeIdx
            END DO !componentIdx
          END DO !gaussIdx
        END DO !boundaryIdx
      END IF

    CASE DEFAULT
      ! Do nothing for other equation set subtypes
    END SELECT

    EXITS("NavierStokes_FiniteElementBoundaryIntegrate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementBoundaryIntegrate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementBoundaryIntegrate

  !
  !================================================================================================================================
  !

  !> Calculate the fluid flux through 3D boundaries for use in problems with coupled solutions (e.g. multidomain)
  SUBROUTINE NavierStokes_CalculateBoundaryFlux(equationsSet,coupledEquationsSet,iteration3D1D, &
    & convergedFlag,absolute3D0DTolerance,relative3D0DTolerance,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EQUATIONS_SET_TYPE), POINTER :: coupledEquationsSet !<A pointer to the coupled equations set (for 3D-1D coupling)
    INTEGER(INTG), INTENT(IN) :: iteration3D1D !<iteration number for the 3D-1D loop if this is a coupled problem
    REAL(DP), INTENT(IN) :: absolute3D0DTolerance !<absolute convergence criteria for 3D-0D coupling
    REAL(DP), INTENT(IN) :: relative3D0DTolerance !<relative convergence criteria for 3D-0D coupling
    LOGICAL, INTENT(INOUT) :: convergedFlag !<convergence flag for 3D-0D coupling
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers !<A pointer to the solver mapping
    TYPE(EquationsType), POINTER :: equations3D
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping3D
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_TYPE), POINTER :: geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable3D
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable,geometricVariable
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition3D
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(BASIS_TYPE), POINTER :: dependentBasis2
    TYPE(BASIS_TYPE), POINTER :: geometricFaceBasis
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(BASIS_TYPE), POINTER :: faceBasis
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: faceQuadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: equationsEquationsSetField
    TYPE(FIELD_TYPE), POINTER :: equationsSetField3D
    TYPE(FIELD_TYPE), POINTER :: dependentField3D
    TYPE(FIELD_TYPE), POINTER :: dependentField1D
    TYPE(FIELD_TYPE), POINTER :: independentField1D
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    INTEGER(INTG) :: faceIdx, faceNumber,elementIdx,nodeNumber,versionNumber
    INTEGER(INTG) :: componentIdx,gaussIdx
    INTEGER(INTG) :: faceNodeIdx, elementNodeIdx
    INTEGER(INTG) :: faceNodeDerivativeIdx, meshComponentNumber
    INTEGER(INTG) :: normalComponentIdx
    INTEGER(INTG) :: boundaryID,numberOfBoundaries,boundaryType,coupledNodeNumber,numberOfGlobalBoundaries
    INTEGER(INTG) :: MPI_IERROR,numberOfComputationalNodes
    INTEGER(INTG) :: i,j,computationalNode
    REAL(DP) :: gaussWeight, normalProjection,elementNormal(3)
    REAL(DP) :: normalDifference,normalTolerance
    REAL(DP) :: courant,maxCourant,toleranceCourant
    REAL(DP) :: velocityGauss(3),faceNormal(3),unitNormal(3),boundaryValue,faceArea,faceVelocity,facePressure
    REAL(DP) :: pressureGauss,faceTraction,mu,muScale,normal,rho,viscousTerm(3)
    REAL(DP) :: dUDXi(3,3),dXiDX(3,3),gradU(3,3),cauchy(3,3),traction(3),normalWave(2)
    REAL(DP) :: localBoundaryFlux(10),localBoundaryArea(10),globalBoundaryFlux(10),globalBoundaryArea(10)
    REAL(DP) :: localBoundaryPressure(10),globalBoundaryPressure(10),globalBoundaryMeanPressure(10)
    REAL(DP) :: localBoundaryNormalStress(10),globalBoundaryNormalStress(10),globalBoundaryMeanNormalStress(10)
    REAL(DP) :: couplingFlow,couplingStress,p1D,q1D,a1D,p3D,stress3DPrevious,stress1DPrevious,tolerance,p0D,q0D
    REAL(DP) :: flowError,pressureError
    LOGICAL :: couple1DTo3D,couple3DTo1D,boundary3D0DFound(10),boundary3D0DConverged(10)
    LOGICAL, ALLOCATABLE :: globalConverged(:)

    REAL(DP), POINTER :: geometricParameters(:)

    ENTERS("NavierStokes_CalculateBoundaryFlux",err,error,*999)

    NULLIFY(decomposition3D)
    NULLIFY(geometricDecomposition)
    NULLIFY(geometricParameters)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis)
    NULLIFY(dependentBasis2)
    NULLIFY(geometricFaceBasis)
    NULLIFY(geometricVariable)
    NULLIFY(equations3D)
    NULLIFY(vectorMatrices)
    NULLIFY(face)
    NULLIFY(faceBasis)
    NULLIFY(faceQuadratureScheme)
    NULLIFY(fieldVariable)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(rhsVector)
    NULLIFY(dependentField3D)
    NULLIFY(dependentField1D)
    NULLIFY(independentField1D)
    NULLIFY(geometricField)
    NULLIFY(equationsEquationsSetField)
    NULLIFY(equationsSetField3D)

    couple1DTo3D = .FALSE.
    couple3DTo1D = .FALSE.
    boundary3D0DFound = .FALSE.

    SELECT CASE(equationsSet%Specification(3))
    ! 3 D   t y p e s :   I n t e g r a t e   b o u n d a r y   v a l u e s
    ! ------------------------------------------------------------------------
    CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)

      ! Get 3D field pointers
      IF(ASSOCIATED(equationsSet)) THEN
        equations3D=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations3D)) THEN
          dependentField3D=>equationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(.NOT.ASSOCIATED(dependentField3D)) THEN
            CALL FlagError("Dependent field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF
        equationsSetField3D=>equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
        IF(.NOT.ASSOCIATED(equationsSetField3D)) THEN
          CALL FlagError("Equations set field (EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF
      ! Check for a coupled equations set
      IF(ASSOCIATED(coupledEquationsSet)) THEN
        dependentField1D=>coupledEquationsSet%DEPENDENT%DEPENDENT_FIELD
        IF(.NOT.ASSOCIATED(dependentField1D)) THEN
          CALL FlagError("Coupled 1D Dependent field is not associated.",err,error,*999)
        END IF
        independentField1D=>coupledEquationsSet%INDEPENDENT%INDEPENDENT_FIELD
        IF(.NOT.ASSOCIATED(independentField1D)) THEN
          CALL FlagError("Coupled 1D Independent field is not associated.",err,error,*999)
        END IF
      END IF
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations3D,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)

      dependentVariable3D=>nonlinearMapping%residualVariables(1)%ptr
      !Get the mesh decomposition and mapping
      decomposition3D=>dependentVariable3D%FIELD%DECOMPOSITION
      elementsMapping3D=>decomposition3D%DOMAIN(decomposition3D%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
      ! Get constant max Courant (CFL) number (default 1.0)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField3D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & 2,toleranceCourant,err,error,*999)
      IF(equationsSet%Specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      END IF

      ! Loop over elements to locate boundary elements
      maxCourant = 0.0_DP
      numberOfBoundaries = 0
      localBoundaryFlux = 0.0_DP
      localBoundaryArea = 0.0_DP
      localBoundaryPressure = 0.0_DP
      localBoundaryNormalStress = 0.0_DP
      DO elementIdx=1,elementsMapping3D%NUMBER_OF_LOCAL
        meshComponentNumber=dependentVariable3D%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        dependentBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementIdx)%BASIS
        decompElement=>DECOMPOSITION3D%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)

        ! Note: if CFL tolerance = 0, we'll skip this step, which speeds things up a bit
        IF (toleranceCourant > ZERO_TOLERANCE) THEN
          ! C F L  c o n d i t i o n   c h e c k
          ! ------------------------------------
          ! Calculate element metrics (courant #, cell Reynolds number)
          CALL NavierStokes_CalculateElementMetrics(equationsSet,elementIdx,err,error,*999)
          ! Get element metrics
          CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & elementIdx,3,courant,err,error,*999)
          IF(courant < -ZERO_TOLERANCE) THEN
            CALL FLAG_WARNING("Negative Courant (CFL) number.",ERR,ERROR,*999)
          END IF
          IF(courant > maxCourant) maxCourant = courant
          ! Check if element CFL number below specified tolerance
          IF(courant > toleranceCourant) THEN
            localError="Element "//TRIM(NUMBER_TO_VSTRING(decompElement%user_number, &
              & "*",ERR,ERROR))//" has violated the CFL condition "//TRIM(NUMBER_TO_VSTRING(courant, &
              & "*",ERR,ERROR))//" <= "//TRIM(NUMBER_TO_VSTRING(toleranceCourant,"*",ERR,ERROR))// &
              & ". Decrease timestep or increase CFL tolerance for the 3D Navier-Stokes problem."
            CALL FlagError(localError,ERR,ERROR,*999)
          END IF
        END IF

        ! B o u n d a r y   n o r m a l   a n d   I D
        ! ----------------------------------------------
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,5,elementNormal(1),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,6,elementNormal(2),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,7,elementNormal(3),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,8,boundaryValue,err,error,*999)
        boundaryID=NINT(boundaryValue)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,9,boundaryValue,err,error,*999)
        boundaryType=NINT(boundaryValue)
        !Check if is a non-wall boundary element
        IF(boundaryID > numberOfBoundaries) numberOfBoundaries=boundaryID
        IF(boundaryID>1) THEN
          faceArea=0.0_DP
          faceVelocity=0.0_DP
          facePressure=0.0_DP
          faceTraction=0.0_DP
          !Get the dependent interpolation parameters
          dependentInterpolationParameters=>equations3D%INTERPOLATION%dependentInterpParameters( &
            & FIELD_U_VARIABLE_TYPE)%PTR
          dependentInterpolatedPoint=>equations3D%INTERPOLATION%dependentInterpPoint( &
            & dependentVariable3D%VARIABLE_TYPE)%PTR
          ! Loop over faces to determine the boundary face contribution
          DO faceIdx=1,dependentBasis%NUMBER_OF_LOCAL_FACES
            !Get the face normal and quadrature information
            IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
              faceNumber=decompElement%ELEMENT_FACES(faceIdx)
            ELSE
              CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
            END IF
            face=>decomposition3D%TOPOLOGY%FACES%FACES(faceNumber)
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(face%BOUNDARY_FACE)) CYCLE

            SELECT CASE(dependentBasis%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              normalComponentIdx=ABS(face%XI_DIRECTION)
            CASE(BASIS_SIMPLEX_TYPE)
              CALL FLAG_WARNING("Boundary flux calculation not yet set up for simplex element types.",err,error,*999)
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasis%TYPE,"*",err,error))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(localError,err,error,*999)
            END SELECT

            faceBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
            faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

            ! Loop over face gauss points
            DO gaussIdx=1,faceQuadratureScheme%NUMBER_OF_GAUSS
              !Use the geometric field to find the face normal and Jacobian for the face integral
              geometricInterpolationParameters=>equations3D%INTERPOLATION%geometricInterpParameters( &
                & FIELD_U_VARIABLE_TYPE)%PTR
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
                & geometricInterpolationParameters,err,error,*999)
              geometricInterpolatedPoint=>equations3D%INTERPOLATION%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                & geometricInterpolatedPoint,err,error,*999)
              pointMetrics=>equations3D%INTERPOLATION%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)

              ! TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
              SELECT CASE(dependentBasis%TYPE)
              CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                ! Make sure this is the boundary face that corresponds with boundaryID (could be a wall rather than inlet/outlet)
                DO componentIdx=1,dependentVariable3D%NUMBER_OF_COMPONENTS-1
                  normalProjection=DOT_PRODUCT(pointMetrics%GU(normalComponentIdx,:),pointMetrics%DX_DXI(componentIdx,:))
                  IF(face%XI_DIRECTION<0) THEN
                    normalProjection=-normalProjection
                  END IF
                  faceNormal(componentIdx)=normalProjection
                END DO !componentIdx
                CALL Normalise(faceNormal,unitNormal,err,error,*999)
                CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
                normalTolerance=0.1_DP
                IF(normalDifference>normalTolerance) EXIT
              CASE(BASIS_SIMPLEX_TYPE)
                faceNormal=unitNormal
              CASE DEFAULT
                localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasis%TYPE,"*",err,error))// &
                  & " is not yet implemented for Navier-Stokes."
                CALL FlagError(localError,err,error,*999)
              END SELECT

              ! C a l c u l a t e   C a u c h y   s t r e s s   o n   c o u p l e d    n o d e s
              ! -----------------------------------------------------------------------------------
              gaussWeight=faceQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)
              ! Get the first partial velocity derivatives and pressure; calculate the Cauchy stress tensor if necessary
              IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
               & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                couple3DTo1D = .TRUE.
                ! Get the constitutive law (non-Newtonian) viscosity based on shear rate if needed
                IF(equationsSet%Specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                  ! Get the gauss point based value returned from the CellML solver
                  CALL Field_ParameterSetGetLocalGaussPoint(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,gaussIdx,elementIdx,1,mu,err,error,*999)
                  mu=mu*muScale
                END IF
                !Get the pressure and velocity interpolation parameters
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
                  & dependentInterpolationParameters,ERR,ERROR,*999)
                CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                 & dependentInterpolatedPoint,ERR,ERROR,*999)

                velocityGauss=0.0_DP
                pressureGauss=0.0_DP
                !Interpolated values at gauss point
                velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)
                pressureGauss=dependentInterpolatedPoint%values(4,NO_PART_DERIV)
                dUDXi = 0.0_DP
                dUDXi(1:3,1)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S1)
                dUDXi(1:3,2)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S2)
                dUDXi(1:3,3)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S3)
                ! Assemble viscous term
                dXiDX=0.0_DP
                dXiDX=pointMetrics%DXI_DX(:,:)
                CALL MatrixProduct(dUDXi,dXiDX,gradU,err,error,*999)
                cauchy = 0.0_DP
                traction = 0.0_DP
                ! Calculate Cauchy stress tensor
                DO i = 1,3
                  DO j = 1,3
                    ! DEBUG: pressure only?
                    cauchy(i,j) = mu*(gradU(i,j)+gradU(j,i))
                  END DO
                END DO
                !DEBUG
                viscousTerm = MATMUL(cauchy,faceNormal)
              ELSE
                !Get interpolated velocity
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
                  & dependentInterpolationParameters,ERR,ERROR,*999)
                CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                 & dependentInterpolatedPoint,ERR,ERROR,*999)
                velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)
                pressureGauss=dependentInterpolatedPoint%values(4,NO_PART_DERIV)
              END IF

              ! I n t e g r a t e    f a c e   a r e a ,   v e l o c i t y   a n d   t r a c t i o n
              ! ----------------------------------------------------------------------------------------
              DO componentIdx=1,dependentVariable3D%NUMBER_OF_COMPONENTS-1
                faceArea=faceArea + ABS(faceNormal(componentIdx)*gaussWeight*pointMetrics%JACOBIAN)
                faceVelocity=faceVelocity+velocityGauss(componentIdx)*faceNormal(componentIdx)*gaussWeight*pointMetrics%JACOBIAN
                facePressure=facePressure+pressureGauss*faceNormal(componentIdx)*gaussWeight*pointMetrics%JACOBIAN
                IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
                 & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                  faceTraction = faceTraction + (viscousTerm(componentIdx) - pressureGauss)* &
                    & faceNormal(componentIdx)*gaussWeight*pointMetrics%JACOBIAN
                END IF
              END DO !componentIdx
            END DO !gaussIdx
          END DO !faceIdx
          localBoundaryFlux(boundaryID) = localBoundaryFlux(boundaryID) + faceVelocity
          localBoundaryArea(boundaryID) = localBoundaryArea(boundaryID) + faceArea
          localBoundaryPressure(boundaryID) = localBoundaryPressure(boundaryID) + facePressure
          localBoundaryNormalStress(boundaryID) = localBoundaryNormalStress(boundaryID)+ faceTraction
        END IF !boundaryIdentifier
      END DO !elementIdx
      ! Distribute any updated element fields
      CALL Field_ParameterSetUpdateStart(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

      ! G a t h e r   v a l u e s   o v e r   t h r e a d s
      ! ------------------------------------------------------
      ! Need to add boundary flux for any boundaries split accross computational nodes
      numberOfGlobalBoundaries = 0
      globalBoundaryFlux = 0.0_DP
      globalBoundaryArea = 0.0_DP
      globalBoundaryPressure = 0.0_DP
      globalBoundaryNormalStress = 0.0_DP
      numberOfComputationalNodes=computationalEnvironment%numberOfComputationalNodes
      IF(numberOfComputationalNodes>1) THEN !use mpi
        CALL MPI_ALLREDUCE(localBoundaryFlux,globalBoundaryFlux,10,MPI_DOUBLE_PRECISION,MPI_SUM,   &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryArea,globalBoundaryArea,10,MPI_DOUBLE_PRECISION,MPI_SUM,   &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(localBoundaryNormalStress,globalBoundaryNormalStress,10,MPI_DOUBLE_PRECISION,MPI_SUM,  &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(localBoundaryPressure,globalBoundaryPressure,10,MPI_DOUBLE_PRECISION,MPI_SUM,  &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(numberOfBoundaries,numberOfGlobalBoundaries,1,MPI_INTEGER,MPI_MAX,  &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      ELSE
        numberOfGlobalBoundaries = numberOfBoundaries
        globalBoundaryFlux = localBoundaryFlux
        globalBoundaryArea = localBoundaryArea
        globalBoundaryPressure = localBoundaryPressure
        globalBoundaryNormalStress = localBoundaryNormalStress
      END IF
      globalBoundaryArea=ABS(globalBoundaryArea)
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          globalBoundaryMeanNormalStress(boundaryID)=globalBoundaryNormalStress(boundaryID)/globalBoundaryArea(boundaryID)
          globalBoundaryMeanPressure(boundaryID)=globalBoundaryPressure(boundaryID)/globalBoundaryArea(boundaryID)
        END IF
      END DO
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          computationalNode = ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
          IF(computationalNode==0) THEN
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  flow:  ", &
              & globalBoundaryFlux(boundaryID),err,error,*999)
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  area:  ", &
              & globalBoundaryArea(boundaryID),err,error,*999)
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  mean normal stress:  ", &
              & globalBoundaryMeanNormalStress(boundaryID),err,error,*999)
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  mean pressure:  ", &
              & globalBoundaryMeanPressure(boundaryID),err,error,*999)
            IF (toleranceCourant > ZERO_TOLERANCE) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Max Courant (CFL) number: ",maxCourant,err,error,*999)
            END IF
          END IF
        ELSE
          localError="Zero or negative area boundary detected on boundary "// &
            & TRIM(NUMBER_TO_VSTRING(boundaryID,"*",ERR,ERROR))//"."
          CALL FlagError(localError,ERR,ERROR,*999)
        END IF
      END DO

    ! 1 D   t y p e s :   p r e p a r e   t o   c o p y    a n y   c o u p l e d    v a l u e s
    ! --------------------------------------------------------------------------------------------
    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
      couple1DTo3D = .TRUE.
      ! Get coupled 3D field pointers
      IF(ASSOCIATED(coupledEquationsSet)) THEN
        equations3D=>coupledEquationsSet%EQUATIONS
        IF(ASSOCIATED(equations3D)) THEN
          dependentField3D=>coupledEquationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(dependentField3D)) THEN
            ! Do nothing
          ELSE
            CALL FlagError("Coupled 3D dependent field is not associated.",err,error,*999)
          END IF
          decomposition3D=>nonlinearMapping%residualVariables(1)%PTR%FIELD%DECOMPOSITION
          IF(ASSOCIATED(decomposition3D)) THEN
            elementsMapping3D=>decomposition3D%DOMAIN(decomposition3D%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
            IF(.NOT.ASSOCIATED(elementsMapping3D)) THEN
              CALL FlagError("Coupled 3D elements mapping not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Coupled 3D decomposition not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF
        equationsSetField3D=>coupledEquationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
        IF(.NOT.ASSOCIATED(equationsSetField3D)) THEN
          CALL FlagError("Equations set field (EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
        END IF
        dependentVariable3D=>nonlinearMapping%residualVariables(1)%PTR
        IF(.NOT.ASSOCIATED(dependentVariable3D)) THEN
          CALL FlagError("Dependent Variable 3D is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Coupled 3D equations set is not associated.",err,error,*999)
      END IF
      ! Check for a 1D equations set
      IF(ASSOCIATED(equationsSet)) THEN
        dependentField1D=>equationsSet%DEPENDENT%DEPENDENT_FIELD
        IF(.NOT.ASSOCIATED(dependentField1D)) THEN
          CALL FlagError("Coupled 1D Dependent field is not associated.",err,error,*999)
        END IF
        independentField1D=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
        IF(.NOT.ASSOCIATED(independentField1D)) THEN
          CALL FlagError("Coupled 1D Independent field is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("1D equations set is not associated.",err,error,*999)
      END IF
    CASE DEFAULT
      localError="Boundary flux calcluation for equations type "//TRIM(NUMBER_TO_VSTRING(equationsSet%Specification(3),"*", &
        & ERR,ERROR))//" is not yet implemented for Navier-Stokes."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    ! ------------------------------------------------------------------------------------
    ! C o p y    i n t e g r a t e d   v a l u e s    t o    t a r g e t    f i e l d s
    ! ------------------------------------------------------------------------------------
    convergedFlag = .TRUE.
    ! Loop over elements again to allocate flux terms to boundary nodes
    DO elementIdx=1,elementsMapping3D%TOTAL_NUMBER_OF_LOCAL
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,5,elementNormal(1),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,6,elementNormal(2),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,7,elementNormal(3),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,8,boundaryValue,err,error,*999)
      boundaryID=NINT(boundaryValue)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,9,boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,11,boundaryValue,err,error,*999)
      coupledNodeNumber=NINT(boundaryValue)
      IF(boundaryID>1) THEN
        meshComponentNumber=2
        decompElement=>decomposition3D%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)
        dependentBasis2=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementIdx)%BASIS

        ! M a p   3 D - 1 D    b o u n d a r y    v a l u e s
        ! --------------------------------------------------------
        IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
          & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
          ! Note: the coupled node number here is the user number- not the local number
          normal = 0.0_DP
          DO componentIdx =1,2
             CALL FIELD_PARAMETER_SET_GET_NODE(independentField1D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & 1,1,coupledNodeNumber,componentIdx,normalWave(componentIdx),err,error,*999)
             IF(ABS(normalWave(componentIdx)) > ZERO_TOLERANCE) normal=normalWave(componentIdx)
          END DO
          IF(ABS(normal)<ZERO_TOLERANCE) THEN
            localError="Characteristic normal wave not specified for couping node " &
              & //TRIM(NUMBER_TO_VSTRING(coupledNodeNumber,"*",ERR,ERROR))//"."
            CALL FlagError(localError,ERR,ERROR,*999)
          END IF
          IF(couple3DTo1D) THEN
            ! Convert integrated 3D values to 1D flow and pressure.
            couplingStress = globalBoundaryMeanNormalStress(boundaryID)
            couplingFlow = globalBoundaryFlux(boundaryID)
            ! Map the coupling flow and stress values from the 3D equations set to the coupled 1D field
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,1,couplingFlow,err,error,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,2,couplingStress,err,error,*999)
          ELSE IF(couple1DTo3D) THEN
            ! Get the flow and pressure values from the coupled 1D equations set
            CALL FIELD_PARAMETER_SET_GET_NODE(dependentField1D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,1,q1D,err,error,*999)
            CALL FIELD_PARAMETER_SET_GET_NODE(dependentField1D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,2,a1D,err,error,*999)
            CALL FIELD_PARAMETER_SET_GET_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,1,p1D,err,error,*999)
            ! DEBUG
            CALL Field_ParameterSetEnsureCreated(dependentField1D,FIELD_U2_VARIABLE_TYPE, &
              & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,ERR,ERROR,*999)
            CALL FIELD_PARAMETER_SET_GET_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,2,stress1DPrevious,err,error,*999)
            CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
            ! Convert 1D flow and pressure to coupling stress and flow
            couplingStress = -p1D
            couplingFlow = q1D
            p3D = couplingStress
            ! Update the coupling flow and stress values from the 3D equations set to the coupled 1D field
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,2,couplingStress,err,error,*999)
            CALL FIELD_PARAMETER_SET_UPDATE_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1,coupledNodeNumber,3,couplingFlow,err,error,*999)
          END IF
        END IF

        ! B o u n d a r y   F a c e    N o r m a l s
        ! --------------------------------------------------
        DO faceIdx=1,dependentBasis2%NUMBER_OF_LOCAL_FACES
          !Get the face normal and quadrature information
          IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
            faceNumber=decompElement%ELEMENT_FACES(faceIdx)
          ELSE
            CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
          END IF
          face=>decomposition3D%TOPOLOGY%FACES%FACES(faceNumber)
          IF(.NOT.(face%BOUNDARY_FACE)) CYCLE
          faceBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
          faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

          ! TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
          SELECT CASE(dependentBasis2%TYPE)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            normalComponentIdx=ABS(face%XI_DIRECTION)
            geometricInterpolationParameters=>equations3D%INTERPOLATION%geometricInterpParameters( &
              & FIELD_U_VARIABLE_TYPE)%PTR
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & geometricInterpolationParameters,err,error,*999)
            geometricInterpolatedPoint=>equations3D%INTERPOLATION%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATE_LOCAL_FACE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,1, &
              & geometricInterpolatedPoint,err,error,*999)
            pointMetrics=>equations3D%INTERPOLATION%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)
            DO componentIdx=1,dependentVariable3D%NUMBER_OF_COMPONENTS-1
              normalProjection=DOT_PRODUCT(pointMetrics%GU(normalComponentIdx,:),pointMetrics%DX_DXI(componentIdx,:))
              IF(face%XI_DIRECTION<0) THEN
                normalProjection=-normalProjection
              END IF
              faceNormal(componentIdx)=normalProjection
            END DO !componentIdx
            CALL Normalise(faceNormal,unitNormal,err,error,*999)
          CASE(BASIS_SIMPLEX_TYPE)
            !still have faceNormal/unitNormal
          CASE DEFAULT
            localError="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis2%TYPE,"*",ERR,ERROR))// &
              & " is not yet implemented for Navier-Stokes."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
          CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
          normalTolerance=0.1_DP
          IF(normalDifference>normalTolerance) CYCLE

          ! U p d a t e    N o d a l   V a l u e s
          ! --------------------------------------------------
          ! Update local nodes with integrated boundary flow values
          DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
            elementNodeIdx=dependentBasis2%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
            DO faceNodeDerivativeIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
              nodeNumber=decomposition3D%DOMAIN(meshComponentNumber)%PTR% &
               & TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)%ELEMENT_NODES(elementNodeIdx)
              versionNumber=1
              IF(couple1DTo3D) THEN
                IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
                  ! Copy coupling stress from 1D to pressure values type on 3D
                  !CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField3D,FIELD_U_VARIABLE_TYPE, &
                  !  & FIELD_VALUES_SET_TYPE,1,1,nodeNumber,4,p3D,err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField3D,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,1,1,nodeNumber,4,p3D,err,error,*999)
                ELSE IF(boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                  localError="Coupling flow boundary type not yet implemented for 3D side of coupled 3D-1D problems."
                  CALL FlagError(localError,ERR,ERROR,*999)
                END IF
              ELSE
                IF(boundaryType==BOUNDARY_CONDITION_FIXED_CELLML) THEN
                  ! Check current values against those passed to the CellML solver
                  CALL Field_ParameterSetGetLocalNode(equationsSetField3D,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,1,1,nodeNumber,1,q0D,err,error,*999)
                  CALL Field_ParameterSetGetLocalNode(equationsSetField3D,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,1,1,nodeNumber,2,p0D,err,error,*999)
                  flowError = globalBoundaryFlux(boundaryID)-q0D
                  pressureError = globalBoundaryMeanPressure(boundaryID)+p0D
                  IF ((ABS(flowError) < absolute3D0DTolerance .AND. ABS(pressureError) < absolute3D0DTolerance)) THEN
                    ! CONVERGED ABSOLUTE TOLERANCE
                  ELSE IF (ABS(flowError/globalBoundaryFlux(boundaryID)) < relative3D0DTolerance .AND. &
                    &  ABS(pressureError/globalBoundaryMeanPressure(boundaryID)) < relative3D0DTolerance) THEN
                    ! CONVERGED RELATIVE TOLERANCE
                  ELSE
                    IF (.NOT. boundary3D0DFound(boundaryID)) THEN
                      CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  flow: ", &
                        & q0D,err,error,*999)
                      CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  pressure: ", &
                        & p0D,err,error,*999)
                      CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  flow error: ", &
                        & flowError,err,error,*999)
                      CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  pressure error: ", &
                        & pressureError,err,error,*999)
                      boundary3D0DFound(boundaryID) = .TRUE.
                    END IF
                    convergedFlag = .FALSE.
                  END IF
                END IF
                ! Update flow rates on 3D boundary equations set field
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,1,globalBoundaryFlux(boundaryID),err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,2,globalBoundaryMeanPressure(boundaryID),err,error,*999)
              END IF
            END DO !nodeDerivativeIdx
          END DO !faceNodeIdx
        END DO !faceIdx
      END IF !boundaryIdentifier
    END DO !elementIdx

    !allocate array for mpi communication
    IF(numberOfComputationalNodes>1) THEN !use mpi
      ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
      CALL MPI_ALLGATHER(convergedFlag,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
       & computationalEnvironment%mpiCommunicator,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
      IF(ALL(globalConverged)) THEN
        convergedFlag = .TRUE.
      ELSE
        convergedFlag = .FALSE.
      END IF
    END IF

   ! Distribute any updated ields
    IF (ASSOCIATED(equationsSetField3D)) THEN
      CALL Field_ParameterSetUpdateStart(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      IF(convergedFlag) THEN
        ! If converged, update 0D initial conditions for the next step.
        CALL Field_ParameterSetsCopy(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & FIELD_PREVIOUS_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
      END IF
    END IF
    IF (ASSOCIATED(dependentField3D)) THEN
      CALL Field_ParameterSetUpdateStart(dependentField3D,FIELD_U_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField3D,FIELD_U_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
    END IF
    IF (ASSOCIATED(dependentField1D)) THEN
      CALL Field_ParameterSetUpdateStart(dependentField1D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField1D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateStart(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateStart(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    END IF

    EXITS("NavierStokes_CalculateBoundaryFlux")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateBoundaryFlux",err,error)
    RETURN 1
  END SUBROUTINE NavierStokes_CalculateBoundaryFlux

  !
  !================================================================================================================================
  !

  !> Update the solution for the 1D solver with boundary conditions from a lumped parameter model defined by CellML.
  !> For more information please see chapter 11 of: L. Formaggia, A. Quarteroni, and A. Veneziani, Cardiovascular mathematics:
  !> modeling and simulation of the circulatory system. Milan; New York: Springer, 2009.
  SUBROUTINE NavierStokes_Couple1D0D(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: iterativeLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,materialsField,independentField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVER_TYPE), POINTER :: solver1D
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivativeIdx,versionIdx,componentIdx,numberOfLocalNodes1D
    INTEGER(INTG) :: solver1dNavierStokesNumber,solverNumber,MPI_IERROR,timestep,iteration
    INTEGER(INTG) :: boundaryNumber,numberOfBoundaries,numberOfComputationalNodes
    INTEGER(INTG) :: dependentDof,boundaryConditionType
    REAL(DP) :: A0_PARAM,E_PARAM,H_PARAM,beta,pCellML,normalWave(2)
    REAL(DP) :: qPrevious,pPrevious,aPrevious,q1d,a1d,qError,aError,couplingTolerance
    LOGICAL :: boundaryNode,boundaryConverged(30),localConverged,MPI_LOGICAL,coupled1D0DBoundary
    LOGICAL, ALLOCATABLE :: globalConverged(:)

    ENTERS("NavierStokes_Couple1D0D",err,error,*999)

    !Get solvers based on the problem type
    SELECT CASE(controlLoop%PROBLEM%specification(3))
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      solverNumber = solver%GLOBAL_NUMBER
      ! In the Navier-Stokes/Characteristic subloop, the Navier-Stokes solver should be the second solver
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      IF(solverNumber == solver1dNavierStokesNumber) THEN
        solver1D=>controlLoop%SUB_LOOPS(2)%ptr%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%ptr
        iterativeLoop=>controlLoop%WHILE_LOOP
        iteration = iterativeLoop%ITERATION_NUMBER
        timestep = controlLoop%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
      ELSE
        localError="The solver number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
         & " does not correspond with the Navier-Stokes solver number for 1D-0D fluid coupling."
        CALL FlagError(localError,err,error,*999)
      END IF
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      solverNumber = solver%GLOBAL_NUMBER
      ! In the Navier-Stokes/Characteristic subloop, the Navier-Stokes solver should be the second solver
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      IF(solverNumber == solver1dNavierStokesNumber) THEN
        solver1D=>controlLoop%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%PTR
        iterativeLoop=>controlLoop%WHILE_LOOP
        iteration = iterativeLoop%ITERATION_NUMBER
        timestep = controlLoop%PARENT_LOOP%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
      ELSE
        localError="The solver number of "//TRIM(NUMBER_TO_VSTRING(solverNumber,"*",err,error))// &
         & " does not correspond with the Navier-Stokes solver number for 1D-0D fluid coupling."
        CALL FlagError(localError,ERR,ERROR,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    couplingTolerance = iterativeLoop%ABSOLUTE_TOLERANCE

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver1D)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          solverEquations=>solver1D%SOLVER_EQUATIONS
          IF(ASSOCIATED(solverEquations)) THEN
            solverMapping=>solverEquations%SOLVER_MAPPING
            IF(ASSOCIATED(solverMapping)) THEN
              equationsSet=>solverMapping%EQUATIONS_SETS(1)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
                dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
                NULLIFY(equations)
                CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                NULLIFY(vectorMapping)
                CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
                NULLIFY(dynamicMapping)
                CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver mapping is not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Solver equations is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control Loop is not associated.",err,error,*999)
    END IF

    !Get the number of local nodes
    domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
      & TOPOLOGY%NODES
    IF(ASSOCIATED(domainNodes)) THEN
      numberOfLocalNodes1D=domainNodes%NUMBER_OF_NODES
    ELSE
      CALL FlagError("Domain nodes are not associated.",err,error,*999)
    END IF

    boundaryNumber = 0
    boundaryConverged = .TRUE.
    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes1D
      nodeNumber = domainNodes%NODES(nodeIdx)%local_number
      !Check for the boundary node
      boundaryNode=dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
        & TOPOLOGY%NODES%NODES(nodeNumber)%BOUNDARY_NODE

      !Get node characteristic wave direction (specifies inlet/outlet)
      coupled1D0DBoundary = .FALSE.
      DO componentIdx=1,2
        CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave(componentIdx),err,error,*999)

        ! Check that this is a 1D-0D boundary
        NULLIFY(fieldVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
        dependentDof = fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
          & NODE_PARAM2DOF_MAP%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
        boundaryConditions=>solver1D%SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable, &
          & boundaryConditionsVariable,err,error,*999)
        boundaryConditionType=boundaryConditionsVariable%CONDITION_TYPES(dependentDof)
        IF (boundaryConditionType == BOUNDARY_CONDITION_FIXED_CELLML) THEN
          coupled1D0DBoundary = .TRUE.
        END IF
      END DO

      !!!-- F i n d   B o u n d a r y   N o d e s --!!!
      IF(ABS(normalWave(1))>ZERO_TOLERANCE .AND. boundaryNode .AND. coupled1D0DBoundary) THEN

        boundaryNumber = boundaryNumber + 1
        boundaryConverged(boundaryNumber) = .FALSE.
        !Get material parameters
        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)
        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)
        CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,3,H_PARAM,err,error,*999)
        beta=(4.0_DP*SQRT(PI)*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)

        ! C u r r e n t   Q 1 D , A 1 D , p C e l l M L   V a l u e s
        ! ------------------------------------------------------------
        !Get Q1D
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,1,q1d,err,error,*999)
        !Get A1D
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,2,a1d,err,error,*999)
        !Get pCellML
        CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,1,pCellML,err,error,*999)

        ! C h e c k  1 D / 0 D   C o n v e r g e n c e   f o r   t h i s   n o d e
        ! -------------------------------------------------------------------------
        IF(iteration == 1 .AND. timestep == 0) THEN
          ! Create the previous iteration field values type on the dependent field if it does not exist
          NULLIFY(fieldVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE)%ptr)) THEN
            CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U_VARIABLE_TYPE, &
             & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
          END IF
          NULLIFY(fieldVariable)
          CALL Field_VariableGet(dependentField,FIELD_U1_VARIABLE_TYPE,fieldVariable,err,error,*999)
          IF(.NOT.ASSOCIATED(fieldVariable%PARAMETER_SETS%SET_TYPE(FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE)%ptr)) THEN
            CALL FIELD_PARAMETER_SET_CREATE(dependentField,FIELD_U1_VARIABLE_TYPE, &
             & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
          END IF
        ELSE
          !Get previous Q1D
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeNumber,1,qPrevious,err,error,*999)
          !Get previous A1D
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeNumber,2,aPrevious,err,error,*999)
          !Get previous pCellML value
          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeNumber,1,pPrevious,err,error,*999)
          ! Check if the boundary interface values have converged
          qError = ABS(qPrevious - q1d)
          aError = ABS(aPrevious - a1d)
          IF( qError < couplingTolerance .AND. aError < couplingTolerance) THEN
            boundaryConverged(boundaryNumber) = .TRUE.
          END IF
        END IF

        ! store current Q and p Boundary values as previous iteration value
        CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,1,q1d,err,error,*999)
        CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,2,a1d,err,error,*999)
        CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,1,pCellML,err,error,*999)

      END IF !Find boundary nodes
    END DO !Loop over nodes
    CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    numberOfBoundaries = boundaryNumber

    IF(solverNumber == solver1dNavierStokesNumber) THEN
      ! ------------------------------------------------------------------
      ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
      ! ------------------------------------------------------------------
      ! Check whether all boundaries on the local process have converged
      IF(numberOfBoundaries == 0 .OR. ALL(boundaryConverged(1:numberOfBoundaries))) THEN
        localConverged = .TRUE.
      ELSE
        localConverged = .FALSE.
      END IF
      ! Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
      numberOfComputationalNodes=computationalEnvironment%numberOfComputationalNodes
      IF(numberOfComputationalNodes>1) THEN !use mpi
        !allocate array for mpi communication
        ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
        CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
         & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,err,error,*999)
        IF(ALL(globalConverged)) THEN
          !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          iterativeLoop%CONTINUE_LOOP=.FALSE.
        END IF
        DEALLOCATE(globalConverged)
      ELSE
        IF(localConverged) THEN
          !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          iterativeLoop%CONTINUE_LOOP=.FALSE.
        END IF
      END IF

      ! If the solution hasn't converged, need to revert field values to pre-solve state
      ! before continued iteration. This will counteract the field updates that occur
      ! in SOLVER_DYNAMIC_MEAN_PREDICTED_CALCULATE. Ignore for initialisation
      IF(timestep == 0) THEN
        iterativeLoop%CONTINUE_LOOP=.FALSE.
      END IF
      IF(iterativeLoop%CONTINUE_LOOP .EQV. .TRUE. ) THEN
        CALL Field_ParameterSetsCopy(dependentField,dynamicMapping%dynamicVariableType,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
        CALL Field_ParameterSetsCopy(dependentField,dynamicMapping%dynamicVariableType,FIELD_PREVIOUS_RESIDUAL_SET_TYPE, &
          & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
      END IF
    END IF

    EXITS("NavierStokes_Couple1D0D")
    RETURN
999 ERRORSEXITS("NavierStokes_Couple1D0D",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_Couple1D0D

  !
  !================================================================================================================================
  !
  !> Check whether Coupled 3D and 1D solvers have converged at boundaries
  SUBROUTINE NavierStokes_Couple3D1D(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: iterativeLoop
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet1D,equationsSet3D
    TYPE(EquationsType), POINTER :: equations1D,equations3D
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping1D,vectorMapping3D
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping1D,dynamicMapping3D
    TYPE(EquationsVectorType), POINTER :: vectorEquations1D,vectorEquations3D
    TYPE(FIELD_TYPE), POINTER :: dependentField1D,dependentField3D,independentField
    TYPE(SOLVER_TYPE), POINTER :: solver1D
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivativeIdx,versionIdx,componentIdx,numberOfLocalNodes1D
    INTEGER(INTG) :: solver1dNavierStokesNumber,MPI_IERROR,timestep,iteration
    INTEGER(INTG) :: boundaryNumber,boundaryType1D,numberOfBoundaries,numberOfComputationalNodes
    INTEGER(INTG) :: solver3dNavierStokesNumber,userNodeNumber,localDof,globalDof,computationalNode
    REAL(DP) :: normalWave(2)
    REAL(DP) :: flow1D,stress1D,flow1DPrevious,stress1DPrevious,flow3D,stress3D,flowError,stressError
    REAL(DP) :: maxStressError,maxFlowError,flowTolerance,stressTolerance,absoluteCouplingTolerance
    REAL(DP) :: absoluteCouplingTolerance2,relativeCouplingTolerance
    LOGICAL :: boundaryNode,boundaryConverged(30),localConverged,globalConverged

    ENTERS("NavierStokes_Couple3D1D",ERR,ERROR,*999)

    !Get solvers based on the problem type
    SELECT CASE(controlLoop%PROBLEM%specification(3))
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      iterativeLoop=>controlLoop%WHILE_LOOP
      iteration = iterativeLoop%ITERATION_NUMBER
      timestep = controlLoop%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
      absoluteCouplingTolerance = iterativeLoop%ABSOLUTE_TOLERANCE
      absoluteCouplingTolerance2 = iterativeLoop%ABSOLUTE_TOLERANCE
      relativeCouplingTolerance = iterativeLoop%RELATIVE_TOLERANCE
      IF(.NOT. ASSOCIATED(iterativeLoop)) THEN
        localError="Iterative loop is not associated for 3D-1D fluid coupling."
        CALL FlagError(localError,ERR,ERROR,*999)
      END IF
      ! 1D solver & equations pointers
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      !TODO: make this more general!
      solver1D=>controlLoop%SUB_LOOPS(1)%PTR%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%PTR
      IF(ASSOCIATED(solver1D)) THEN
        boundaryConditions=>solver1D%SOLVER_EQUATIONS%BOUNDARY_CONDITIONS
        IF(ASSOCIATED(boundaryConditions)) THEN
          equationsSet1D => solver1D%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
          IF(ASSOCIATED(equationsSet1D)) THEN
            dependentField1D=>equationsSet1D%DEPENDENT%DEPENDENT_FIELD
            IF(ASSOCIATED(dependentField1D)) THEN
              independentField=>equationsSet1D%INDEPENDENT%INDEPENDENT_FIELD
              IF(.NOT. ASSOCIATED(independentField)) THEN
                localError="Independent field 1D is not associated."
                CALL FlagError(localError,ERR,ERROR,*999)
              END IF
            ELSE
              CALL FlagError("Dependent field 1D is not associated.",ERR,ERROR,*999)
            END IF
          ELSE
            CALL FlagError("Equations set 1D is not associated.",err,error,*999)
          END IF
        ELSE
          localError="Boundary conditions are not associated for 3D-1D fluid coupling."
          CALL FlagError(localError,ERR,ERROR,*999)
        END IF
      ELSE
        localError="The 1D Navier-Stokes solver is not associated for 3D-1D fluid coupling."
        CALL FlagError(localError,ERR,ERROR,*999)
      END IF
      ! 3D equations pointers
      solver3dNavierStokesNumber = 1
      ! TODO: make this more general!
      equationsSet3D=>controlLoop%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(solver3dNavierStokesNumber)%PTR% &
        & SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
      IF(ASSOCIATED(equationsSet3D)) THEN
        dependentField3D=>equationsSet3D%DEPENDENT%DEPENDENT_FIELD
        IF(.NOT. ASSOCIATED(dependentField3D)) THEN
          localError="Dependent field 3D is not associated."
          CALL FlagError(localError,ERR,ERROR,*999)
        END IF
      ELSE
        CALL FlagError("Equations set 3D is not associated.",err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(controlLoop%PROBLEM%specification(3),"*",err,error))// &
        & " is not valid for 3D-1D Navier-Stokes fluid coupling."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    NULLIFY(equations1D)
    CALL EquationsSet_EquationsGet(equationsSet1D,equations1D,err,error,*999)
    NULLIFY(vectorEquations1D)
    CALL Equations_VectorEquationsGet(equations1D,vectorEquations1D,err,error,*999)
    NULLIFY(vectorMapping1D)
    CALL EquationsVector_VectorMappingGet(vectorEquations1D,vectorMapping1D,err,error,*999)
    NULLIFY(dynamicMapping1D)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping1D,dynamicMapping1D,err,error,*999)
    NULLIFY(equations3D)
    CALL EquationsSet_EquationsGet(equationsSet3D,equations3D,err,error,*999)
    NULLIFY(vectorEquations3D)
    CALL Equations_VectorEquationsGet(equations3D,vectorEquations3D,err,error,*999)
    NULLIFY(vectorMapping3D)
    CALL EquationsVector_VectorMappingGet(vectorEquations3D,vectorMapping3D,err,error,*999)
    NULLIFY(dynamicMapping3D)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping3D,dynamicMapping3D,err,error,*999)

    !Get the number of local nodes

    domainNodes=>dependentField1D%DECOMPOSITION%DOMAIN(dependentField1D%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
      & TOPOLOGY%NODES
    IF(ASSOCIATED(domainNodes)) THEN
      numberOfLocalNodes1D=domainNodes%NUMBER_OF_NODES
    ELSE
      CALL FlagError("Domain nodes are not associated.",err,error,*999)
    END IF

    boundaryNumber = 0
    maxStressError = 0.0_DP
    maxFlowError = 0.0_DP
    boundaryConverged = .TRUE.
    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes1D
      nodeNumber = domainNodes%NODES(nodeIdx)%LOCAL_NUMBER
      userNodeNumber = domainNodes%NODES(nodeIdx)%USER_NUMBER
      !Check for the boundary node- go to next if not boundary
      boundaryNode=dependentField1D%DECOMPOSITION%DOMAIN(dependentField1D%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
        & TOPOLOGY%NODES%NODES(nodeNumber)%BOUNDARY_NODE
      IF(.NOT. boundaryNode) CYCLE

      !Get node characteristic wave direction (specifies inlet/outlet)
      DO componentIdx=1,2
        CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave(componentIdx),err,error,*999)
      END DO

      !!!-- F i n d   B o u n d a r y   N o d e s --!!!
      IF(ABS(normalWave(1))>ZERO_TOLERANCE .OR. ABS(normalWave(2))>ZERO_TOLERANCE) THEN
        ! Check that this is a coupled 3D-1D boundary
        boundaryType1D = 0
        NULLIFY(fieldVariable)
        CALL Field_VariableGet(dependentField1D,FIELD_U_VARIABLE_TYPE,fieldVariable,ERR,ERROR,*999)
        CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
        CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField1D,FIELD_U_VARIABLE_TYPE,versionIdx,derivativeIdx, &
          & userNodeNumber,1,localDof,globalDof,err,error,*999)
        IF(boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
          boundaryType1D = BOUNDARY_CONDITION_COUPLING_FLOW
        END IF
        CALL FIELD_COMPONENT_DOF_GET_USER_NODE(dependentField1D,FIELD_U_VARIABLE_TYPE,versionIdx,derivativeIdx, &
          & userNodeNumber,2,localDof,globalDof,err,error,*999)
        IF(boundaryConditionsVariable%CONDITION_TYPES(globalDof)==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
          IF(boundaryType1D==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
            localError="Boundary type for node number "//TRIM(NUMBER_TO_VSTRING(userNodeNumber,"*",err,error))// &
              & " is set as both FLOW and STRESS for 3D-1D Navier-Stokes fluid coupling."
            CALL FlagError(localError,ERR,ERROR,*999)
          END IF
          boundaryType1D = BOUNDARY_CONDITION_COUPLING_STRESS
        END IF
        ! If this is not a coupled 3D-1D boundary, go to the next node
        IF(boundaryType1D==0) CYCLE

        boundaryNumber = boundaryNumber + 1
        boundaryConverged(boundaryNumber) = .FALSE.
        CALL Field_ParameterSetEnsureCreated(dependentField1D,FIELD_U1_VARIABLE_TYPE, &
          & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetEnsureCreated(dependentField1D,FIELD_U2_VARIABLE_TYPE, &
          & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
        !Get 1D flow and stress from current iteration
        CALL Field_ParameterSetGetLocalNode(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,2,stress1D,err,error,*999)
        CALL Field_ParameterSetGetLocalNode(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,3,flow1D,err,error,*999)
        !Get 3D flow and stress
        CALL Field_ParameterSetGetLocalNode(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,1,flow3D,err,error,*999)
        CALL Field_ParameterSetGetLocalNode(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeNumber,2,stress3D,err,error,*999)
        IF (boundaryNumber ==1) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"------3D-1D----- iteration:  ",iteration,err,error,*999)
        END IF
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"  node:  ",userNodeNumber,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    1D flow:  ",flow1D,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    3D flow:  ",flow3D,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    1D stress: ",stress1D,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"    3D stress: ",stress3D,err,error,*999)

        ! C h e c k  3 D / 1 D   C o n v e r g e n c e   f o r   t h i s   n o d e
        ! -------------------------------------------------------------------------
        ! Check if the boundary interface values have converged
        flowError = flow1D - flow3D
        stressError = stress1D - stress3D
        flowTolerance = ABS(relativeCouplingTolerance*flow3D)
        stressTolerance = ABS(relativeCouplingTolerance*stress1D)

        !DEBUG check if within 1%- maybe set up a relative tolerance?
        IF(ABS(flowError) < flowTolerance .AND. ABS(stressError) < stressTolerance) THEN
          boundaryConverged(boundaryNumber) = .TRUE.
        ELSE IF( ABS(flowError) < absoluteCouplingTolerance .AND. ABS(stressError) < absoluteCouplingTolerance2) THEN
          boundaryConverged(boundaryNumber) = .TRUE.
        END IF
        maxFlowError = MAX(ABS(maxFlowError),ABS(flowError))
        maxStressError = MAX(ABS(maxStressError),ABS(stressError))

        ! Set current iteration values to previous
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField1D,FIELD_U1_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,1,flow3D,err,error,*999)
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField1D,FIELD_U1_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,2,stress3D,err,error,*999)
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,2,stress1D,err,error,*999)
        CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField1D,FIELD_U2_VARIABLE_TYPE, &
         & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,3,flow1D,err,error,*999)

      END IF !Find boundary nodes
    END DO !Loop over nodes
    CALL Field_ParameterSetUpdateStart(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField1D,FIELD_U1_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateStart(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField1D,FIELD_U2_VARIABLE_TYPE,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
      & err,error,*999)
    numberOfBoundaries = boundaryNumber

    ! ------------------------------------------------------------------
    ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
    ! ------------------------------------------------------------------
    ! Check whether all boundaries on the local process have converged
    localConverged=.FALSE.
    globalConverged=.FALSE.
    IF(numberOfBoundaries == 0 .OR. ALL(boundaryConverged(1:numberOfBoundaries))) THEN
      localConverged = .TRUE.
    END IF
    ! Need to check that boundaries have converged globally (on all domains) if this is a MPI problem
    numberOfComputationalNodes=computationalEnvironment%numberOfComputationalNodes
    IF(numberOfComputationalNodes>1) THEN !use mpi
      !allocate array for mpi communication

      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
      CALL MPI_ALLREDUCE(localConverged,globalConverged,1,MPI_LOGICAL,MPI_LAND,computationalEnvironment%mpiCommunicator,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      IF(globalConverged) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D coupling converged; # iterations: ", &
          & iteration,err,error,*999)
        iterativeLoop%CONTINUE_LOOP=.FALSE.
      ELSE
        computationalNode = ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"Rank ",computationalNode," 3D/1D max flow error:  ", &
          & maxFlowError,err,error,*999)
        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"Rank ",computationalNode," 3D/1D max stress error:  ", &
          & maxStressError,err,error,*999)
      END IF
    ELSE
      IF(localConverged) THEN
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D coupling converged; # iterations: ", &
          & iteration,err,error,*999)
        iterativeLoop%CONTINUE_LOOP=.FALSE.
      ELSE
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D max flow error:  ", &
          & maxFlowError,err,error,*999)
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D max stress error: ", &
          & maxStressError,err,error,*999)
      END IF
    END IF

    ! If the solution hasn't converged, need to revert field values to pre-solve state
    ! before continued iteration. This will counteract the field updates that occur
    ! in SOLVER_DYNAMIC_MEAN_PREDICTED_CALCULATE. Ignore for initialisation
    IF(timestep == 0) THEN
      iterativeLoop%CONTINUE_LOOP=.FALSE.
    END IF
    IF(iterativeLoop%CONTINUE_LOOP .EQV. .TRUE. ) THEN
      ! Reset 1D values
      CALL Field_ParameterSetsCopy(dependentField1D,dynamicMapping1D%dynamicVariableType, &
        & FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
      CALL Field_ParameterSetsCopy(dependentField1D,dynamicMapping1D%dynamicVariableType, &
        & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP,ERR,ERROR,*999)
      ! Reset 3D values
      CALL Field_ParameterSetsCopy(dependentField3D,dynamicMapping3D%dynamicVariableType, &
        & FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP,ERR,ERROR,*999)
      CALL Field_ParameterSetsCopy(dependentField3D,dynamicMapping3D%dynamicVariableType, &
        & FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP,ERR,ERROR,*999)
    END IF

    EXITS("NavierStokes_Couple3D1D")
    RETURN
999 ERRORSEXITS("NavierStokes_Couple3D1D",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_Couple3D1D

  !
  !================================================================================================================================
  !

  !> Check convergence of
  SUBROUTINE NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(SOLVER_TYPE), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(CONTROL_LOOP_WHILE_TYPE), POINTER :: iterativeLoop
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField,independentField,materialsField
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(SOLVER_TYPE), POINTER :: solver1DNavierStokes
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: nodeNumber,nodeIdx,derivativeIdx,versionIdx,componentIdx,i
    INTEGER(INTG) :: solver1dNavierStokesNumber,solverNumber
    INTEGER(INTG) :: branchNumber,numberOfBranches,numberOfComputationalNodes,numberOfVersions
    INTEGER(INTG) :: MPI_IERROR,timestep,iteration,outputIteration
    REAL(DP) :: couplingTolerance,l2ErrorW(30),wPrevious(2,7),wNavierStokes(2,7),wCharacteristic(2,7),wError(2,7)
    REAL(DP) :: l2ErrorQ(100),qCharacteristic(7),qNavierStokes(7),wNext(2,7)
    REAL(DP) :: totalErrorWPrevious,startTime,stopTime,currentTime,timeIncrement
    REAL(DP) :: l2ErrorA(100),aCharacteristic(7),aNavierStokes(7),totalErrorW,totalErrorQ,totalErrorA
    REAL(DP) :: totalErrorMass,totalErrorMomentum
    REAL(DP) :: rho,alpha,normalWave,A0_PARAM,E_PARAM,H_PARAM,beta,aNew,penaltyCoeff
    LOGICAL :: branchConverged(100),localConverged,MPI_LOGICAL,boundaryNode,fluxDiverged
    LOGICAL, ALLOCATABLE :: globalConverged(:)

    ENTERS("NavierStokes_CoupleCharacteristics",err,error,*999)

    SELECT CASE(controlLoop%PROBLEM%specification(3))
    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      solver1DNavierStokes=>controlLoop%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%ptr
      CALL CONTROL_LOOP_TIMES_GET(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
       & timestep,outputIteration,err,error,*999)
      iteration = controlLoop%WHILE_LOOP%ITERATION_NUMBER
      iterativeLoop=>controlLoop%WHILE_LOOP
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      solver1DNavierStokes=>controlLoop%PARENT_LOOP%SUB_LOOPS(2)%ptr%SOLVERS%SOLVERS(solver1dNavierStokesNumber)%ptr
      iterativeLoop=>controlLoop%WHILE_LOOP
      iteration = iterativeLoop%ITERATION_NUMBER
      timestep = controlLoop%PARENT_LOOP%PARENT_LOOP%TIME_LOOP%ITERATION_NUMBER
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    solverNumber = solver%GLOBAL_NUMBER
    couplingTolerance = iterativeLoop%ABSOLUTE_TOLERANCE

    IF(ASSOCIATED(controlLoop)) THEN
      IF(ASSOCIATED(solver1DNavierStokes)) THEN
        IF(ASSOCIATED(controlLoop%PROBLEM)) THEN
          solverEquations=>solver1DNavierStokes%SOLVER_EQUATIONS
          IF(ASSOCIATED(solverEquations)) THEN
            solverMapping=>solverEquations%SOLVER_MAPPING
            IF(ASSOCIATED(solverMapping)) THEN
              equationsSet=>solverMapping%EQUATIONS_SETS(1)%ptr
              IF(ASSOCIATED(equationsSet)) THEN
                dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
                independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
                materialsField=>equationsSet%MATERIALS%MATERIALS_FIELD
              ELSE
                CALL FlagError("Equations set is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver mapping is not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Solver equations is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Control Loop is not associated.",err,error,*999)
    END IF

    !Get the number of local nodes
    domainNodes=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
      & TOPOLOGY%NODES
    branchNumber = 0
    branchConverged = .TRUE.
    fluxDiverged = .FALSE.
    totalErrorQ = 0.0_DP
    totalErrorA = 0.0_DP
    totalErrorW = 0.0_DP
    totalErrorMass = 0.0_DP
    totalErrorMomentum = 0.0_DP
    totalErrorWPrevious = 0.0_DP
    l2ErrorQ = 0.0_DP
    l2ErrorA = 0.0_DP
    l2ErrorW = 0.0_DP

    ! Get material constants
    CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & 2,rho,err,error,*999)
    CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & 3,alpha,err,error,*999)
    CALL Field_ParameterSetGetConstant(equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
      & FIELD_VALUES_SET_TYPE,1,penaltyCoeff,err,error,*999)

    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
      nodeNumber = domainNodes%NODES(nodeIdx)%local_number
      derivativeIdx = 1
      numberOfVersions=domainNodes%NODES(nodeNumber)%DERIVATIVES(derivativeIdx)%numberOfVersions
      boundaryNode=domainNodes%NODES(nodeNumber)%BOUNDARY_NODE

      !DEBUG
      CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE, &
       & FIELD_VALUES_SET_TYPE,1,1,nodeNumber,1,normalWave,err,error,*999)
      !Find branch nodes
      IF(numberOfVersions>1) THEN
        branchNumber = branchNumber + 1
        branchConverged(branchNumber) = .FALSE.

        wError = 0.0_DP
        i = 0
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            i = i +1
            CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
            IF(ABS(normalWave)>ZERO_TOLERANCE) THEN

              ! Get the previously set characteristic (W) for this timestep-
              !  if this is the first iteration it will be based on extrapolated values
              !  otherwise it will come from the last iteration of this subroutine.
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,componentIdx,wPrevious(componentIdx,versionIdx),err,error,*999)

              !Get material parameters
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,1,A0_PARAM,err,error,*999)
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,2,E_PARAM,err,error,*999)
              CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,derivativeIdx,nodeNumber,3,H_PARAM,err,error,*999)
              beta=(4.0_DP*SQRT(PI)*E_PARAM*H_PARAM)/(3.0_DP*A0_PARAM)

              ! Get current Q,A values based on N-S solve
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,1,qNavierStokes(versionIdx),err,error,*999)
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,2,aNavierStokes(versionIdx),err,error,*999)

              ! Calculate the characteristic based on the values converged upon by the
              !  N-S solver at this iteration.
              wNavierStokes(componentIdx,versionIdx)= ((qNavierStokes(versionIdx)/aNavierStokes(versionIdx))+ &
               & normalWave*4.0_DP*SQRT(beta/(2.0_DP*rho))*(aNavierStokes(versionIdx)**(0.25_DP) - (A0_PARAM)**(0.25_DP)))

              IF(boundaryNode) THEN
                aNew = (1.0_DP/(beta/(2.0_DP*rho)))**2.0_DP*((wNavierStokes(componentIdx,versionIdx))/8.0_DP+ &
                 & SQRT(beta/(2.0_DP*rho))*((A0_PARAM)**0.25_DP))**4.0_DP
                CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_PREVIOUS_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber, &
                 & 2,aNew,err,error,*999)
              END IF

              ! Get characteristic (flux conserving) Q,A values
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,1,qCharacteristic(versionIdx),err,error,*999)
              CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
               & versionIdx,derivativeIdx,nodeNumber,2,aCharacteristic(versionIdx),err,error,*999)

              ! Calculate the characteristic based on the upwind values
              wCharacteristic(componentIdx,versionIdx)= ((qCharacteristic(versionIdx)/aCharacteristic(versionIdx))+ &
               & normalWave*4.0_DP*SQRT((beta/(2.0_DP*rho)))*(aCharacteristic(versionIdx)**(0.25_DP) - (A0_PARAM)**(0.25_DP)))
            END IF
          END DO
        END DO

        ! Evaluate error between current and previous Q,A values
        IF(numberOfVersions > 1 ) THEN
          CALL L2Norm(qNavierStokes-qCharacteristic,l2ErrorQ(branchNumber),err,error,*999)
          CALL L2Norm(aNavierStokes-aCharacteristic,l2ErrorA(branchNumber),err,error,*999)
        END IF
        ! Check if the branch values have converged
        IF((ABS(l2ErrorQ(branchNumber)) < couplingTolerance) .AND. (ABS(l2ErrorA(branchNumber)) < couplingTolerance)) THEN
          branchConverged(branchNumber) = .TRUE.
        END IF
        totalErrorQ = totalErrorQ + l2ErrorQ(branchNumber)
        totalErrorA = totalErrorA + l2ErrorA(branchNumber)

        wNext = ((wNavierStokes + wCharacteristic)/2.0_DP)
        ! If N-S/C w values did not converge re-solve with new w.
        IF(numberOfVersions > 1) THEN
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE, &
               & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
              IF(ABS(normalWave)>ZERO_TOLERANCE) THEN
                !Update W value
                CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                 & versionIdx,derivativeIdx,nodeNumber,componentIdx,wNext(componentIdx,versionIdx),err,error,*999)
              END IF
            END DO
          END DO
        END IF

      END IF !Find boundary nodes
    END DO !Loop over nodes
    CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateStart(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    numberOfBranches = branchNumber

    ! ------------------------------------------------------------------
    ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
    ! ------------------------------------------------------------------
    ! Check whether all branches on the local process have converged
    IF(numberOfBranches == 0 .OR. ALL(branchConverged(1:numberOfBranches))) THEN
      localConverged = .TRUE.
    ELSE
      localConverged = .FALSE.
    END IF
    ! Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
    numberOfComputationalNodes=computationalEnvironment%numberOfComputationalNodes
    IF(numberOfComputationalNodes>1) THEN !use mpi
      !allocate array for mpi communication
      ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
      CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
       & computationalEnvironment%mpiCommunicator,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,err,error,*999)
      IF(ALL(globalConverged)) THEN
        !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        controlLoop%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
      END IF
      DEALLOCATE(globalConverged)
    ELSE
      IF(localConverged) THEN
        !CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        controlLoop%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
      END IF
    END IF

    EXITS("NavierStokes_CoupleCharacteristics")
    RETURN
999 ERRORSEXITS("NavierStokes_CoupleCharacteristics",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_CoupleCharacteristics

  !
  !================================================================================================================================
  !

  !>Calculated the rate of deformation (shear rate) for a navier-stokes finite element equations set.
  SUBROUTINE NavierStokes_ShearRateCalculate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_TYPE), POINTER :: materialsField
    INTEGER(INTG) :: elementIdx,decompositionLocalElementNumber
    INTEGER(INTG) :: gaussIdx
    INTEGER(INTG) :: meshComponentNumber,numberOfDimensions,i,j,userElementNumber
    INTEGER(INTG) :: localElementNumber,startElement,stopElement
    REAL(DP) :: gaussWeight,shearRate,rateOfStrainMag,strainRate
    REAL(DP) :: dUdXi(3,3),dXidX(3,3),dUdX(3,3),dUdXTrans(3,3),D(3,3),velocityGauss(3)
    REAL(DP) :: shearRateDefault
    LOGICAL :: ghostElement,elementExists,defaultUpdate

    ENTERS("NavierStokes_ShearRateCalculate",err,error,*999)

    CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...Calculating shear rate...",err,error,*999)

    NULLIFY(decomposition)
    NULLIFY(dependentBasis)
    NULLIFY(equations)
    NULLIFY(quadratureScheme)
    NULLIFY(fieldVariable)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(dependentField)
    NULLIFY(materialsField)
    NULLIFY(vectorEquations)
    NULLIFY(vectorMapping)
    NULLIFY(nonlinearMapping)

    ! Some preliminary sanity checks
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
      dependentVariable=>nonlinearMapping%residualVariables(1)%ptr
      meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
      !Get the mesh decomposition and mapping
      decomposition=>dependentVariable%FIELD%DECOMPOSITION
      elementsMapping=>decomposition%DOMAIN(decomposition%MESH_COMPONENT_NUMBER)%ptr%MAPPINGS%ELEMENTS
      fieldVariable=>nonlinearMapping%residualVariables(1)%ptr
      numberOfDimensions=fieldVariable%NUMBER_OF_COMPONENTS - 1
      dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(dependentVariable%VARIABLE_TYPE)%ptr
      geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
      defaultUpdate=.FALSE.

      ! Loop over internal and boundary elements, skipping ghosts
      startElement = elementsMapping%INTERNAL_START
      stopElement = elementsMapping%BOUNDARY_FINISH
      ! Loop over internal and boundary elements
      DO elementIdx=startElement,stopElement
        localElementNumber=elementsMapping%DOMAIN_LIST(elementIdx)
        userElementNumber = elementsMapping%LOCAL_TO_GLOBAL_MAP(localElementNumber)
        !Check computational node for elementIdx
        elementExists=.FALSE.
        ghostElement=.TRUE.
        CALL DECOMPOSITION_TOPOLOGY_ELEMENT_CHECK_EXISTS(decomposition%TOPOLOGY, &
          & userElementNumber,elementExists,decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(ghostElement) THEN
          CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Ghost: ",userElementNumber,err,error,*999)
        END IF

        IF(elementExists) THEN
          dependentBasis=>decomposition%DOMAIN(meshComponentNumber)%ptr%TOPOLOGY%ELEMENTS%ELEMENTS(localElementNumber)%BASIS
          quadratureScheme=>dependentBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber,equations%interpolation% &
            & dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

          ! Loop over gauss points
          DO gaussIdx=1,quadratureScheme%NUMBER_OF_GAUSS
            !Get interpolated velocity
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & dependentInterpolatedPoint,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & geometricInterpolatedPoint,err,error,*999)
            pointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_VOLUME_TYPE,pointMetrics,err,error,*999)
            gaussWeight=quadratureScheme%GAUSS_WEIGHTS(gaussIdx)
            !Interpolated values at gauss point
            dXidX=0.0_DP
            dUdXi=0.0_DP
            velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)

            dUdXi(1:3,1)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S1)
            dUdXi(1:3,2)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S2)
            IF(numberOfDimensions == 3) THEN
              dUdXi(1:3,3)=dependentInterpolatedPoint%VALUES(1:3,PART_DERIV_S3)
            ELSE
              dUdXi(1:3,3)=0.0_DP
            END IF
            dXidX=pointMetrics%DXI_DX(:,:)
            dUdX=0.0_DP
            dUdXTrans=0.0_DP
            strainRate=0.0_DP

            CALL MatrixProduct(dUdXi,dXidX,dUdX,err,error,*999) !dU/dX = dU/dxi * dxi/dX (deformation gradient tensor)
            DO i=1,3
              DO j=1,3
                D(i,j) = (dUdX(i,j) + dUdX(j,i))/2.0_DP
              END DO
            END DO
            rateOfStrainMag= 0.5_DP*(D(1,1)*D(1,1)+D(1,2)*D(1,2)+D(1,3)*D(1,3)+D(2,1)*D(2,1)+ &
              & D(2,2)*D(2,2)+D(2,3)*D(2,3) + D(3,1)*D(3,1)+ D(3,2)*D(3,2)+D(3,3)*D(3,3))
            shearRate=SQRT(ABS(2.0_DP*rateOfStrainMag))

            ! TODO: hard coded fix for milliseconds to seconds conversion- need to fix this to be more general!
            shearRate = shearRate*1000.0_DP
            CALL Field_ParameterSetUpdateLocalGaussPoint(materialsField,FIELD_V_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,gaussIdx,localElementNumber,2,shearRate,err,error,*999)

          END DO !gaussIdx
        END IF ! check for ghost element
        IF(defaultUpdate .EQV. .TRUE.) THEN
          EXIT
        END IF
      END DO !elementIdx
      CALL Field_ParameterSetUpdateStart(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

      IF(defaultUpdate .EQV. .TRUE.) THEN
        shearRateDefault=1.0E-10_DP
        CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Setting default shear field values...", &
         & shearRateDefault,err,error,*999)
        CALL FIELD_COMPONENT_VALUES_INITIALISE(materialsField,FIELD_V_VARIABLE_TYPE, &
         & FIELD_VALUES_SET_TYPE,1,shearRateDefault,err,error,*999)
      END IF

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
        & " is not valid for shear rate calculation in a Navier-Stokes equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_ShearRateCalculate")
    RETURN
999 ERRORSEXITS("NavierStokes_ShearRateCalculate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_ShearRateCalculate

  !
  !================================================================================================================================
  !

  !>Pre-residual evaluation a navier-stokes finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementPreResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    LOGICAL :: convergedFlag

    ENTERS("NavierStokes_FiniteElementPreResidualEvaluate",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        ! Shear rate should either be calculated here to update at each minor iteration
        ! or during post solve so it is updated once per timestep
        !CALL NavierStokes_ShearRateCalculate(equationsSet,err,error,*999)
      CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
         convergedFlag = .FALSE.
         CALL NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*999)
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        ! & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
         & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        !Do nothing
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_FiniteElementPreResidualEvaluate")
    RETURN
999 ERRORS("NavierStokes_FiniteElementPreResidualEvaluate",err,error)
    EXITS("NavierStokes_FiniteElementPreResidualEvaluate")
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementPreResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE NavierStokes_ControlLoopPostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: subloop,subloop2,subloop3,iterativeWhileLoop2,iterativeWhileLoop3
    TYPE(SOLVER_TYPE), POINTER :: navierStokesSolver,navierStokesSolver3D,navierStokesSolver1D,solver
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping,solverMapping2
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet,equationsSet2,coupledEquationsSet
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: numberOfSolvers,solverIdx,solverIdx2,equationsSetIdx,equationsSetIdx2
    INTEGER(INTG) :: subloopIdx,subloopIdx2,subloopIdx3,iteration3D1D
    REAL(DP) :: absolute3D0DTolerance,relative3D0DTolerance
    LOGICAL :: convergedFlag
    character(70) :: label

    ENTERS("NavierStokes_ControlLoopPostLoop",err,error,*999)

    NULLIFY(equationsSet)
    NULLIFY(coupledEquationsSet)
    NULLIFY(dependentField)
    NULLIFY(fieldVariable)
    convergedFlag = .FALSE.
    absolute3D0DTolerance = 0.0_DP
    relative3D0DTolerance = 0.0_DP

    IF(ASSOCIATED(controlLoop)) THEN
      SELECT CASE(controlLoop%PROBLEM%specification(3))
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_PGM_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        ! Do nothing
      CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%PTR
          IF(ASSOCIATED(navierStokesSolver)) THEN
            NULLIFY(coupledEquationsSet)
            equationsSet=>navierStokesSolver%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
            IF(ASSOCIATED(equationsSet)) THEN
              dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
              IF(dependentField%DECOMPOSITION%CALCULATE_FACES) THEN
                CALL NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,ERR,ERROR,*999)
              END IF
            ELSE
              CALL FlagError("Equations set is not associated.",err,error,*999)
            ENDIF
            CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
          ELSE
            CALL FlagError("Navier-Stokes solver not associated.",ERR,ERROR,*999)
          END IF
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NUMBER_TO_VSTRING(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for the specified Navier-Stokes problem type."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          ! Do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          ! Global time loop - export data
          navierStokesSolver=>controlLoop%SUB_LOOPS(1)%ptr%SOLVERS%SOLVERS(2)%ptr
          CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%ptr
          CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,err,error,*999)
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NumberToVString(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          ! CellML simple loop - do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          ! Global time loop - export data
          navierStokesSolver=>controlLoop%SUB_LOOPS(1)%ptr%SUB_LOOPS(2)%ptr%SOLVERS%SOLVERS(2)%ptr
          CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(navierStokesSolver,err,error,*999)
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          ! Couple 1D/0D loop
          IF(controlLoop%CONTROL_LOOP_LEVEL==2) THEN
            navierStokesSolver=>controlLoop%SUB_LOOPS(2)%ptr%SOLVERS%SOLVERS(2)%ptr
            ! update 1D/0D coupling parameters and check convergence
            CALL NavierStokes_Couple1D0D(controlLoop,navierStokesSolver,err,error,*999)
          ! Couple Navier-Stokes/Characteristics loop
          ELSE IF(controlLoop%CONTROL_LOOP_LEVEL==3) THEN
            navierStokesSolver=>controlLoop%SOLVERS%SOLVERS(2)%ptr
            CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,err,error,*999)
          ELSE
            localError="The while loop level of "//TRIM(NumberToVString(controlLoop%CONTROL_LOOP_LEVEL,"*",err,error))// &
              & " is invalid for a Coupled 1D0D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END IF
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NumberToVString(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
         & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
        SELECT CASE(controlLoop%LOOP_TYPE)
        ! Simple loops- could be 3D or 0D
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          numberOfSolvers=controlLoop%SOLVERS%NUMBER_OF_SOLVERS
          DO solverIdx=1,numberOfSolvers
            solver=>controlLoop%SOLVERS%SOLVERS(solverIdx)%PTR
            SELECT CASE(solver%SOLVE_TYPE)
            CASE(SOLVER_DYNAMIC_TYPE)
              solverMapping=>solver%SOLVER_EQUATIONS%SOLVER_MAPPING
              IF(ASSOCIATED(solverMapping)) THEN
                DO equationsSetIdx = 1,solverMapping%NUMBER_OF_EQUATIONS_SETS
                  equationsSet=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%PTR
                  IF(ASSOCIATED(equationsSet)) THEN
                    SELECT CASE(equationsSet%Specification(3))
                    ! --- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
                    CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                       & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                      NULLIFY(coupledEquationsSet)
                      ! TODO: This is a bit of a hack- make more robust
                      NULLIFY(iterativeWhileLoop2)
                      NULLIFY(iterativeWhileLoop3)
                      ! 3D-1D iteration
                      iteration3D1D=controlLoop%PARENT_LOOP%WHILE_LOOP%ITERATION_NUMBER
                      ! 1D-0D
                      IF (controlLoop%PROBLEM%specification(3)==PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE) THEN
                        iterativeWhileLoop2=>controlLoop%PARENT_LOOP%SUB_LOOPS(1)%PTR
                        ! 1D
                        CALL CONTROL_LOOP_SUB_LOOP_GET(iterativeWhileLoop2,2,iterativeWhileLoop3,ERR,ERROR,*999)
                        IF (.NOT. ASSOCIATED(iterativeWhileLoop3)) THEN
                          CALL FlagError("Could not find 1D subloop!",err,error,*999)
                        END IF
                        IF(iterativeWhileLoop3%NUMBER_OF_SUB_LOOPS==0) THEN
                          DO solverIdx2=1,iterativeWhileLoop3%SOLVERS%NUMBER_OF_SOLVERS
                            solverMapping2=>iterativeWhileLoop3%SOLVERS%SOLVERS(solverIdx2)%PTR%SOLVER_EQUATIONS%SOLVER_MAPPING
                            IF(ASSOCIATED(solverMapping2)) THEN
                              DO equationsSetIdx2 = 1,solverMapping2%NUMBER_OF_EQUATIONS_SETS
                                equationsSet2=>solverMapping2%EQUATIONS_SETS(equationsSetIdx)%PTR
                                IF(ASSOCIATED(equationsSet2)) THEN
                                  SELECT CASE(equationsSet2%specification(3))
                                  CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                                    &  EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
                                    IF(ASSOCIATED(coupledEquationsSet)) THEN
                                      localError="Coupled 3D-1D equations set already found for multiscale Navier-Stokes problem."
                                      CALL FlagError(localError,err,error,*999)
                                    ELSE
                                      coupledEquationsSet=>equationsSet2
                                    END IF
                                  CASE DEFAULT
                                    ! Do nothing
                                  END SELECT
                                ELSE
                                  CALL FlagError("Equations set 2 is not associated.",err,error,*999)
                                END IF
                              END DO ! equations set loop 2
                            ELSE
                              CALL FlagError("Solver mapping 2 is not associated.",ERR,ERROR,*999)
                            END IF !
                          END DO ! solverIdx2
                        ELSE
                          CALL FlagError("Unrecognized subloop pattern for 3D-1D-0D!.",ERR,ERROR,*999)
                        END IF ! check for futher subloops
                      END IF
                      controlLoop%PARENT_LOOP%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
                    CASE DEFAULT
                      localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%Specification(3),"*",err,error))// &
                        & " is not valid for a dynamic solver in a simple loop for a multiscale Navier-Stokes problem."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  END IF
                END DO ! equations set loop
              ELSE
                CALL FlagError("Solver mapping is not associated.",ERR,ERROR,*999)
              END IF !
            CASE(SOLVER_DAE_TYPE)
              ! CellML solver simple loop - do nothing
            CASE DEFAULT
              localError="The solve type of "//TRIM(NUMBER_TO_VSTRING(solver%SOLVE_TYPE,"*",ERR,ERROR))// &
                & " is invalid for a simple loop in a Navier-Stokes multiscale problem."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
          END DO ! solverIdx
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          ! Global time loop - export data from all dynamic solvers and equations sets
          ! Check for subloops two layers down
          NULLIFY(subloop)
          NULLIFY(subloop2)
          NULLIFY(subloop3)
          IF(controlLoop%NUMBER_OF_SUB_LOOPS>0) THEN
            DO subloopIdx=1,controlLoop%NUMBER_OF_SUB_LOOPS
              subloop=>controlLoop%SUB_LOOPS(subloopIdx)%PTR
              IF(subloop%NUMBER_OF_SUB_LOOPS>0) THEN
                DO subloopIdx2=1,subloop%NUMBER_OF_SUB_LOOPS
                  subloop2=>subloop%SUB_LOOPS(subloopIdx2)%PTR

                  IF(subloop2%NUMBER_OF_SUB_LOOPS>0) THEN
                    ! 1D output
                    DO subloopIdx3=1,subloop2%NUMBER_OF_SUB_LOOPS
                      subloop3=>subloop2%SUB_LOOPS(subloopIdx3)%PTR
                      IF(subloop3%NUMBER_OF_SUB_LOOPS>0) THEN
                        CALL FlagError("Unrecognized number of subloops in 3D-1D-0D.",ERR,ERROR,*999)
                      ELSE
                        IF (ASSOCIATED(subloop3%SOLVERS)) THEN
                          DO solverIdx=1,subloop3%SOLVERS%NUMBER_OF_SOLVERS
                            solver=>subloop3%SOLVERS%SOLVERS(solverIdx)%PTR
                            IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                              solverMapping=>solver%SOLVER_EQUATIONS%SOLVER_MAPPING
                              IF(ASSOCIATED(solverMapping)) THEN
                                CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(solver,err,error,*999)
                              ELSE
                                CALL FlagError("1D Solver mapping is not associated where expected for multiscale.",ERR,ERROR,*999)
                              END IF
                            END IF !
                          END DO ! Solver loop
                        END IF
                      END IF
                    END DO

                  ELSE
                    ! 3D output
                    IF (ASSOCIATED(subloop2%SOLVERS)) THEN
                      DO solverIdx=1,subloop2%SOLVERS%NUMBER_OF_SOLVERS
                        solver=>subloop2%SOLVERS%SOLVERS(solverIdx)%PTR
                        IF(solver%SOLVE_TYPE==SOLVER_DYNAMIC_TYPE) THEN
                          solverMapping=>solver%SOLVER_EQUATIONS%SOLVER_MAPPING
                          IF(ASSOCIATED(solverMapping)) THEN
                            CALL NAVIER_STOKES_POST_SOLVE_OUTPUT_DATA(solver,err,error,*999)
                          ELSE
                            CALL FlagError("3D Solver mapping is not associated where expected for multiscale.",ERR,ERROR,*999)
                          END IF
                        END IF !
                      END DO ! solverIdx
                    END IF
                  END IF
                END DO ! subloop2
              ELSE
                CALL FlagError("Navier-Stokes Multiscale problem level 2 loop has no subloops.",ERR,ERROR,*999)
              END IF
            END DO ! subloop 1
          ELSE
            CALL FlagError("Navier-Stokes Multiscale problem time loop has no subloops.",ERR,ERROR,*999)
          END IF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          ! TODO: here we get loop type by label- may need to think of a more robust way to do this
          CALL CONTROL_LOOP_LABEL_GET(controlLoop,label,err,error,*999)
          SELECT CASE(label)
          CASE("3D-0D Iterative Loop")
            ! Will handle 3D-0D coupling in calc boundary flux
            ! update 3D/D coupling parameters and check convergence
            !CALL NavierStokes_Couple3D0D(controlLoop,err,error,*999)
          CASE("3D-1D Iterative Loop")
            ! update 1D/3D coupling parameters and check convergence
            CALL NavierStokes_Couple3D1D(controlLoop,err,error,*999)
          CASE("1D-0D Iterative Coupling Convergence Loop")
            NULLIFY(navierStokesSolver1D,navierStokesSolver3D)
            ! TODO: make this more general!e
            navierStokesSolver1D=>controlLoop%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(2)%PTR
            ! update 1D/0D coupling parameters and check convergence
            CALL NavierStokes_Couple1D0D(controlLoop,navierStokesSolver1D,err,error,*999)
            IF (controlLoop%WHILE_LOOP%CONTINUE_LOOP .EQV. .FALSE.) THEN
              ! TODO: make this more general!
              navierStokesSolver3D=>controlLoop%PARENT_LOOP%SUB_LOOPS(2)%PTR%SOLVERS%SOLVERS(1)%PTR
              iteration3D1D = controlLoop%PARENT_LOOP%WHILE_LOOP%ITERATION_NUMBER
              IF (ASSOCIATED(navierStokesSolver3D)) THEN
                equationsSet=>navierStokesSolver1D%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                coupledEquationsSet=>navierStokesSolver3D%SOLVER_EQUATIONS%SOLVER_MAPPING%EQUATIONS_SETS(1)%PTR
                CALL NavierStokes_CalculateBoundaryFlux(equationsSet,coupledEquationsSet,iteration3D1D, &
                  & convergedFlag,absolute3D0DTolerance,relative3D0DTolerance,ERR,ERROR,*999)
              ELSE
                CALL FlagError("Could not locate 3D solver from 1d-0d subloop",ERR,ERROR,*999)
              END IF
            END IF
          CASE("1D Iterative Loop")
            ! No longer using the NS/C coupling loop
            controlLoop%WHILE_LOOP%CONTINUE_LOOP=.FALSE.
          CASE DEFAULT
            localError="The iterative loop label of "//label// &
              & " does not correspond to a recognised loop type for a Navier-Stokes multiscale problem."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NUMBER_TO_VSTRING(controlLoop%LOOP_TYPE,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(controlLoop%PROBLEM%specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    END IF

    EXITS("NavierStokes_ControlLoopPostLoop")
    RETURN
999 ERRORSEXITS("NavierStokes_ControlLoopPostLoop",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_ControlLoopPostLoop

  !
  !================================================================================================================================
  !

  !>Updates boundary conditions for multiscale fluid problems
  SUBROUTINE NavierStokes_UpdateMultiscaleBoundary(equationsSet,boundaryConditions,timeIncrement,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    REAL(DP), INTENT(IN) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(DOMAIN_TYPE), POINTER :: dependentDomain
    TYPE(EquationsType), POINTER :: equations
    TYPE(FIELD_TYPE), POINTER :: dependentField,materialsField,independentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: rho,A0,H0,E,beta,pExternal,lengthScale,timeScale,massScale
    REAL(DP) :: pCellml,qCellml,ABoundary,QBoundary,W1,W2,ACellML,normalWave(2,7),norm
    REAL(DP) :: Q3D,A3D,p3D
    REAL(DP), POINTER :: Impedance(:),Flow(:)
    INTEGER(INTG) :: nodeIdx,versionIdx,derivativeIdx,componentIdx,numberOfVersions,numberOfLocalNodes
    INTEGER(INTG) :: dependentDof,boundaryConditionType,k

    ENTERS("NavierStokes_UpdateMultiscaleBoundary",err,error,*999)

    NULLIFY(dependentDomain)
    NULLIFY(equations)
    NULLIFY(geometricField)
    NULLIFY(dependentField)
    NULLIFY(independentField)
    NULLIFY(materialsField)
    NULLIFY(fieldVariable)

    ! Preliminary checks; get field and domain pointers
    IF(ASSOCIATED(equationsSet)) THEN
      SELECT CASE(equationsSet%Specification(3))
      CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_STREE1D0D_SUBTYPE, &
        &  EQUATIONS_SET_STREE1D0D_ADV_SUBTYPE, &
        &  EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        equations=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations)) THEN
          geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
          independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
          dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(dependentField)) THEN
            dependentDomain=>dependentField%DECOMPOSITION%DOMAIN(dependentField% &
              & DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR
            IF(.NOT.ASSOCIATED(dependentDomain)) THEN
              CALL FlagError("Dependent domain is not associated.",err,error,*999)
            END IF
          ELSE
            CALL FlagError("Geometric field is not associated.",err,error,*999)
          END IF
          materialsField=>equationsSet%materials%MATERIALS_FIELD
          IF(.NOT.ASSOCIATED(materialsField)) THEN
            CALL FlagError("Materials field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        END IF
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(equationsSet%Specification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes multiscale boundary update."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    END IF

    SELECT CASE(equationsSet%specification(3))
    !!!-- 1 D    E q u a t i o n s   S e t --!!!
    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

      independentField=>equationsSet%INDEPENDENT%INDEPENDENT_FIELD
      numberOfLocalNodes=dependentDomain%TOPOLOGY%NODES%NUMBER_OF_NODES
      derivativeIdx=1
      !Get constant material parameters
      CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
      CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,4,pExternal,err,error,*999)
      !Get materials scale factors
      CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,5,lengthScale,err,error,*999)
      CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,6,timeScale,err,error,*999)
      CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE, &
        & FIELD_VALUES_SET_TYPE,7,massScale,err,error,*999)

      !!!--  L o o p   o v e r   l o c a l    n o d e s  --!!!
      DO nodeIdx=1,numberOfLocalNodes
        numberOfVersions=dependentDomain%TOPOLOGY%NODES%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%numberOfVersions

        !Get normal wave direction
        normalWave=0.0_DP
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            CALL Field_ParameterSetGetLocalNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)
          END DO
        END DO
        !!!-- F i n d   b o u n d a r y    n o d e s --!!!
        IF(ABS(normalWave(1,1)) > ZERO_TOLERANCE .OR. ABS(normalWave(2,1))> ZERO_TOLERANCE) THEN
          CALL L2Norm(normalWave(:,1),norm,err,error,*999)
          IF(numberOfVersions == 1 .AND. norm > ZERO_TOLERANCE) THEN
            versionIdx = 1
            !Get material parameters
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,1,A0,err,error,*999)
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,2,E,err,error,*999)
            CALL Field_ParameterSetGetLocalNode(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,3,H0,err,error,*999)
            beta=(4.0_DP*(SQRT(PI))*E*H0)/(3.0_DP*A0)
            NULLIFY(fieldVariable)
            CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
            ! Get the boundary condition type for the dependent field primitive variables (Q,A)
            DO componentIdx=1,2
              dependentDof = fieldVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP% &
               & NODE_PARAM2DOF_MAP%NODES(nodeIdx)%DERIVATIVES(derivativeIdx)%VERSIONS(versionIdx)
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(boundaryConditions,fieldVariable,boundaryConditionsVariable,err,error,*999)
              boundaryConditionType=boundaryConditionsVariable%CONDITION_TYPES(dependentDof)
              SELECT CASE(boundaryConditionType)

              ! N o n - r e f l e c t i n g   B o u n d a r y
              ! ----------------------------------------------------
              CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING)
                ! Outlet - set W2 to 0, get W1 from the extrapolated value
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,1,W1,err,error,*999)
                  W2 = 0.0_DP
                ! Inlet - set W1 to 0, get W2 from the extrapolated value
                ELSE
                  W1 = 0.0_DP
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,2,W2,err,error,*999)
                END IF
                ! Calculate new area value based on W1, W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                 & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area 1D non-reflecting boundary detected at node "//TRIM(NUMBER_TO_VSTRING(nodeIdx, &
                   & "*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,ABoundary,err,error,*999)

              ! C o u p l e d   C e l l M L  ( 0 D )   B o u n d a r y
              ! ------------------------------------------------------------
              CASE(BOUNDARY_CONDITION_FIXED_CELLML)
                !Get qCellML used in pCellML calculation
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,1,QCellML,err,error,*999)
                !Get pCellML if this is a coupled problem
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,2,pCellml,err,error,*999)
                ! Convert pCellML from SI base units specified in CellML file to scaled units (e.g., kg/(m.s^2) --> g/(mm.ms^2))
                ! pCellml = pCellml*massScale/(lengthScale*(timeScale**2.0_DP))
                ! Convert pCellML --> A0D
                ACellML=((pCellml-pExternal)/beta+SQRT(A0))**2.0_DP
                !  O u t l e t
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,1,W1,err,error,*999)
                  ! Calculate W2 from 0D domain
                  W2 = QCellml/ACellml - 4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP - A0**0.25_DP)
                !  I n l e t
                ELSE
                  ! Calculate W1 from 0D domain
                  W1 = QCellml/ACellml + 4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP - A0**0.25_DP)
                  ! Calculate W2 from 1D domain
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,2,W2,err,error,*999)
                END IF
                ! Calculate new area value based on W1,W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                  & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area coupled 1D0D boundary detected at node "//TRIM(NUMBER_TO_VSTRING(nodeIdx, &
                   & "*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,ABoundary,err,error,*999)

              ! C o u p l e d    3 D    B o u n d a r y
              ! ------------------------------------------------------------
              CASE(BOUNDARY_CONDITION_COUPLING_FLOW)
                !Get q3D
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,1,Q3D,err,error,*999)
                QBoundary = Q3D
                ! Set new Q value based on 3D value
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,1,QBoundary,err,error,*999)
              CASE(BOUNDARY_CONDITION_COUPLING_STRESS)
                !Get q3D
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,1,Q3D,err,error,*999)
                !Get p3D if this is a coupled problem
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,2,p3D,err,error,*999)
                ! Convert p3D --> A3D
                A3D=((p3D-pExternal)/beta+SQRT(A0))**2.0_DP
                !  O u t l e t
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,1,W1,err,error,*999)
                  ! Calculate W2 from 0D domain
                  W2 = Q3D/A3D - 4.0_DP*SQRT(beta/(2.0_DP*rho))*(A3D**0.25_DP - A0**0.25_DP)
                !  I n l e t
                ELSE
                  ! Calculate W1 from 0D domain
                  W1 = Q3D/A3D + 4.0_DP*SQRT(beta/(2.0_DP*rho))*(A3D**0.25_DP - A0**0.25_DP)
                  ! Calculate W2 from 1D domain
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                   & versionIdx,derivativeIdx,nodeIdx,2,W2,err,error,*999)
                END IF
                ! Calculate new area value based on W1,W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                  & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                !DEBUG
                !ABoundary=A3D
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area coupled 3D1D boundary detected at node "//TRIM(NUMBER_TO_VSTRING(nodeIdx, &
                   & "*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                END IF
                ! Set new A value based on 3D boundary and update dof
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,ABoundary,err,error,*999)

              ! S t r u c t u r e d   T r e e   B o u n d a r y
              ! ------------------------------------------------------------
              CASE(BOUNDARY_CONDITION_FIXED_STREE)
                NULLIFY(Impedance)
                NULLIFY(Flow)
                !Get qCellML used in pCellML calculation
                CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,1,QCellML,err,error,*999)
                !Get flow function
                CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & Impedance,err,error,*999)
                !Get impedance function
                CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & Flow,err,error,*999)
                pCellml = 0.0_DP
                DO k=1,size(Flow)
                  pCellml=pCellml+Flow(k)*Impedance(k)*timeIncrement
                END DO
                ! Convert pCellML --> A0D
                ACellML=((pCellml-pExternal)/beta+SQRT(A0))**2.0_DP
                !  O u t l e t
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,derivativeIdx,nodeIdx,1,W1,err,error,*999)
                  ! Calculate W2 from 0D domain
                  W2 = QCellml/ACellml-4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP-A0**0.25_DP)
                !  I n l e t
                ELSE
                  ! Calculate W1 from 0D domain
                  W1 = QCellml/ACellml+4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP-A0**0.25_DP)
                  ! Calculate W2 from 1D domain
                  CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & versionIdx,derivativeIdx,nodeIdx,2,W2,err,error,*999)
                END IF
                ! Calculate new area value based on W1, W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                 & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField,FIELD_U_VARIABLE_TYPE, &
                 & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx,2,ABoundary,err,error,*999)

              CASE(BOUNDARY_CONDITION_FREE, &
                 & BOUNDARY_CONDITION_FIXED, &
                 & BOUNDARY_CONDITION_FIXED_INLET, &
                 & BOUNDARY_CONDITION_FIXED_OUTLET, &
                 & BOUNDARY_CONDITION_FIXED_FITTED)
                ! Do nothing

              CASE DEFAULT
                localError="The boundary conditions type "//TRIM(NUMBER_TO_VSTRING(boundaryConditionType,"*",err,error))// &
                 & " is not valid for a coupled characteristic problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            END DO ! componentIdx

          END IF ! boundary node
        END IF ! branch or boundary node
      END DO !Loop over nodes
      ! Update distributed fields
      CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

    !!!-- 3 D    E q u a t i o n s   S e t --!!!
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_PGM_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
      ! Do nothing

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
       & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_UpdateMultiscaleBoundary")
    RETURN
999 ERRORSEXITS("NavierStokes_UpdateMultiscaleBoundary",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_UpdateMultiscaleBoundary

  !
  !================================================================================================================================
  !

  !> Calculate the fluid flux through 3D boundaries for use in problems with coupled solutions (e.g. multidomain)
  SUBROUTINE NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations3D
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping3D
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable3D
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition3D
    TYPE(DECOMPOSITION_ELEMENT_TYPE), POINTER :: decompElement
    TYPE(BASIS_TYPE), POINTER :: dependentBasis
    TYPE(BASIS_TYPE), POINTER :: dependentBasis2
    TYPE(DECOMPOSITION_FACE_TYPE), POINTER :: face
    TYPE(BASIS_TYPE), POINTER :: faceBasis
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: faceQuadratureScheme
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: geometricInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: pointMetrics
    TYPE(FIELD_TYPE), POINTER :: equationsSetField3D
    TYPE(FIELD_TYPE), POINTER :: dependentField3D
    INTEGER(INTG) :: faceIdx, faceNumber,elementIdx,nodeNumber,versionNumber
    INTEGER(INTG) :: componentIdx,gaussIdx
    INTEGER(INTG) :: faceNodeIdx, elementNodeIdx
    INTEGER(INTG) :: faceNodeDerivativeIdx, meshComponentNumber
    INTEGER(INTG) :: boundaryID,numberOfBoundaries,boundaryType,coupledNodeNumber,numberOfGlobalBoundaries
    INTEGER(INTG) :: MPI_IERROR,numberOfComputationalNodes
    INTEGER(INTG) :: computationalNode,xiDirection(3),orientation
    REAL(DP) :: gaussWeight, elementNormal(3)
    REAL(DP) :: normalDifference,normalTolerance
    REAL(DP) :: courant,maxCourant,toleranceCourant,boundaryValueTemp
    REAL(DP) :: velocityGauss(3),faceNormal(3),unitNormal(3),boundaryValue,faceArea,faceVelocity,facePressure
    REAL(DP) :: pressureGauss,faceTraction,mu,muScale
    REAL(DP) :: localBoundaryFlux(10),localBoundaryArea(10),globalBoundaryFlux(10),globalBoundaryArea(10)
    REAL(DP) :: localBoundaryPressure(10),globalBoundaryPressure(10),globalBoundaryMeanPressure(10)
    REAL(DP) :: localBoundaryNormalStress(10),globalBoundaryMeanNormalStress(10)
    REAL(DP) :: p0D,q0D
    LOGICAL :: boundary3D0DFound(10)
    LOGICAL :: convergedFlag !<convergence flag for 3D-0D coupling
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_CalculateBoundaryFlux3D0D",err,error,*999)

    NULLIFY(decomposition3D)
    NULLIFY(decompElement)
    NULLIFY(dependentBasis)
    NULLIFY(dependentBasis2)
    NULLIFY(equations3D)
    NULLIFY(face)
    NULLIFY(faceBasis)
    NULLIFY(faceQuadratureScheme)
    NULLIFY(dependentInterpolatedPoint)
    NULLIFY(dependentInterpolationParameters)
    NULLIFY(geometricInterpolatedPoint)
    NULLIFY(geometricInterpolationParameters)
    NULLIFY(dependentField3D)
    NULLIFY(equationsSetField3D)

    boundary3D0DFound = .FALSE.

    SELECT CASE(equationsSet%Specification(3))
    ! 3 D   t y p e s :   I n t e g r a t e   b o u n d a r y   v a l u e s
    ! ------------------------------------------------------------------------
    CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      ! Get 3D field pointers
      IF(ASSOCIATED(equationsSet)) THEN
        equations3D=>equationsSet%EQUATIONS
        IF(ASSOCIATED(equations3D)) THEN
          dependentField3D=>equationsSet%DEPENDENT%DEPENDENT_FIELD
          IF(.NOT.ASSOCIATED(dependentField3D)) THEN
            CALL FlagError("Dependent field is not associated.",err,error,*999)
          END IF
        ELSE
          CALL FlagError("Equations set equations is not associated.",err,error,*999)
        END IF
        equationsSetField3D=>equationsSet%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
        IF(.NOT.ASSOCIATED(equationsSetField3D)) THEN
          CALL FlagError("Equations set field (EQUATIONS_SET_FIELD_FIELD) is not associated.",err,error,*999)
        END IF
      ELSE
        CALL FlagError("Solvers is not associated.",err,error,*999)
      END IF

      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations3D,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      dependentVariable3D=>nonlinearMapping%residualVariables(1)%ptr
      !Get the mesh decomposition and mapping
      decomposition3D=>dependentVariable3D%FIELD%DECOMPOSITION
      elementsMapping3D=>decomposition3D%DOMAIN(decomposition3D%MESH_COMPONENT_NUMBER)%PTR%MAPPINGS%ELEMENTS
      ! Get constant max Courant (CFL) number (default 1.0)
      CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSetField3D,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & 2,toleranceCourant,err,error,*999)
      IF(equationsSet%Specification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FIELD_PARAMETER_SET_GET_CONSTANT(equationsSet%MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
          & FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      END IF

      ! Loop over elements to locate boundary elements
      maxCourant = 0.0_DP
      numberOfBoundaries = 0
      localBoundaryFlux = 0.0_DP
      localBoundaryArea = 0.0_DP
      localBoundaryPressure = 0.0_DP
      localBoundaryNormalStress = 0.0_DP
      DO elementIdx=1,elementsMapping3D%NUMBER_OF_LOCAL
        meshComponentNumber=dependentVariable3D%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        dependentBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementIdx)%BASIS
        decompElement=>DECOMPOSITION3D%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)

        ! Note: if CFL tolerance = 0, we'll skip this step, which speeds things up a bit
        IF (toleranceCourant > ZERO_TOLERANCE) THEN
          ! C F L  c o n d i t i o n   c h e c k
          ! ------------------------------------
          ! Calculate element metrics (courant #, cell Reynolds number)
          CALL NavierStokes_CalculateElementMetrics(equationsSet,elementIdx,err,error,*999)
          ! Get element metrics
          CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
           & elementIdx,3,courant,err,error,*999)
          IF(courant < -ZERO_TOLERANCE) THEN
            CALL FLAG_WARNING("Negative Courant (CFL) number.",ERR,ERROR,*999)
          END IF
          IF(courant > maxCourant) maxCourant = courant
          ! Check if element CFL number below specified tolerance
          IF(courant > toleranceCourant) THEN
            localError="Element "//TRIM(NUMBER_TO_VSTRING(decompElement%user_number, &
              & "*",ERR,ERROR))//" has violated the CFL condition "//TRIM(NUMBER_TO_VSTRING(courant, &
              & "*",ERR,ERROR))//" <= "//TRIM(NUMBER_TO_VSTRING(toleranceCourant,"*",ERR,ERROR))// &
              & ". Decrease timestep or increase CFL tolerance for the 3D Navier-Stokes problem."
            CALL FlagError(localError,ERR,ERROR,*999)
          END IF
        END IF

        ! B o u n d a r y   n o r m a l   a n d   I D
        ! ----------------------------------------------
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,5,elementNormal(1),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,6,elementNormal(2),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,7,elementNormal(3),err,error,*999)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,8,boundaryValueTemp,err,error,*999)
        boundaryID=NINT(boundaryValueTemp)
        CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & elementIdx,9,boundaryValue,err,error,*999)
        boundaryType=NINT(boundaryValue)
        !Check if is a non-wall boundary element
        IF(boundaryID > numberOfBoundaries) numberOfBoundaries=boundaryID
        IF(boundaryID>1) THEN
          faceArea=0.0_DP
          faceVelocity=0.0_DP
          facePressure=0.0_DP
          faceTraction=0.0_DP
          ! Loop over faces to determine the boundary face contribution
          DO faceIdx=1,dependentBasis%NUMBER_OF_LOCAL_FACES
            !Get the face normal and quadrature information
            IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
              faceNumber=decompElement%ELEMENT_FACES(faceIdx)
            ELSE
              CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
            END IF
            face=>decomposition3D%TOPOLOGY%FACES%FACES(faceNumber)
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(face%BOUNDARY_FACE)) CYCLE

            xiDirection = 0.0_DP
            SELECT CASE(dependentBasis%TYPE)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              xiDirection(3)=ABS(face%XI_DIRECTION)
            CASE(BASIS_SIMPLEX_TYPE)
              CALL FLAG_WARNING("Boundary flux calculation not yet set up for simplex element types.",ERR,ERROR,*999)
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis%TYPE,"*",ERR,ERROR))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(localError,ERR,ERROR,*999)
            END SELECT
            faceBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
            faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR

            !Use the geometric field to find the face normal and Jacobian for the face integral
            geometricInterpolationParameters=>equations3D%INTERPOLATION%geometricInterpParameters( &
              & FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber, &
              & geometricInterpolationParameters,err,error,*999)
            geometricInterpolatedPoint=>equations3D%INTERPOLATION%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR

            dependentInterpolationParameters=>equations3D%INTERPOLATION%dependentInterpParameters( &
              & FIELD_U_VARIABLE_TYPE)%PTR
            CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber, &
              & dependentInterpolationParameters,err,error,*999)
            dependentInterpolatedPoint=>equations3D%INTERPOLATION%dependentInterpPoint( &
              & dependentVariable3D%VARIABLE_TYPE)%PTR

            xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
            xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
            orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*face%XI_DIRECTION)
            pointMetrics=>equations3D%INTERPOLATION%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR
            ! Loop over face gauss points
            DO gaussIdx=1,faceQuadratureScheme%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
                & geometricInterpolatedPoint,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE,pointMetrics,err,error,*999)

              ! Make sure this is the boundary face that corresponds with boundaryID (could be a wall rather than inlet/outlet)
              CALL CrossProduct(pointMetrics%DX_DXI(:,1),pointMetrics%DX_DXI(:,2),faceNormal,ERR,ERROR,*999)
              faceNormal = faceNormal*orientation
              CALL Normalise(faceNormal,unitNormal,err,error,*999)
              CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
              normalTolerance=0.1_DP
              IF(normalDifference>normalTolerance) EXIT

              !Get interpolated velocity and pressure
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
               & dependentInterpolatedPoint,ERR,ERROR,*999)
              velocityGauss=dependentInterpolatedPoint%values(1:3,NO_PART_DERIV)
              pressureGauss=dependentInterpolatedPoint%values(4,NO_PART_DERIV)

              gaussWeight=faceQuadratureScheme%GAUSS_WEIGHTS(gaussIdx)
              ! I n t e g r a t e    f a c e   a r e a ,   v e l o c i t y   a n d   p r e s s u r e
              ! ----------------------------------------------------------------------------------------
              faceArea=faceArea + gaussWeight*pointMetrics%JACOBIAN
              facePressure=facePressure + pressureGauss*gaussWeight*pointMetrics%JACOBIAN
              DO componentIdx=1,dependentVariable3D%NUMBER_OF_COMPONENTS-1
                faceVelocity=faceVelocity+velocityGauss(componentIdx)*unitNormal(componentIdx)*gaussWeight*pointMetrics%JACOBIAN
              END DO !componentIdx
            END DO !gaussIdx
          END DO !faceIdx
          localBoundaryFlux(boundaryID) = localBoundaryFlux(boundaryID) + faceVelocity
          localBoundaryArea(boundaryID) = localBoundaryArea(boundaryID) + faceArea
          localBoundaryPressure(boundaryID) = localBoundaryPressure(boundaryID) + facePressure
        END IF !boundaryIdentifier
      END DO !elementIdx
      ! Distribute any updated element fields
      CALL Field_ParameterSetUpdateStart(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)

      ! G a t h e r   v a l u e s   o v e r   t h r e a d s
      ! ------------------------------------------------------
      ! Need to add boundary flux for any boundaries split accross computational nodes
      numberOfGlobalBoundaries = 0
      globalBoundaryFlux = 0.0_DP
      globalBoundaryArea = 0.0_DP
      globalBoundaryPressure = 0.0_DP
      numberOfComputationalNodes=computationalEnvironment%numberOfComputationalNodes
      IF(numberOfComputationalNodes>1) THEN !use mpi
        CALL MPI_ALLREDUCE(localBoundaryFlux,globalBoundaryFlux,10,MPI_DOUBLE_PRECISION,MPI_SUM,   &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(localBoundaryArea,globalBoundaryArea,10,MPI_DOUBLE_PRECISION,MPI_SUM,   &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(localBoundaryPressure,globalBoundaryPressure,10,MPI_DOUBLE_PRECISION,MPI_SUM,  &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
        CALL MPI_ALLREDUCE(numberOfBoundaries,numberOfGlobalBoundaries,1,MPI_INTEGER,MPI_MAX,  &
	 & computationalEnvironment%mpiCommunicator,MPI_IERROR)
        CALL MPI_ERROR_CHECK("MPI_ALLREDUCE",MPI_IERROR,ERR,ERROR,*999)
      ELSE
        numberOfGlobalBoundaries = numberOfBoundaries
        globalBoundaryFlux = localBoundaryFlux
        globalBoundaryArea = localBoundaryArea
        globalBoundaryPressure = localBoundaryPressure
      END IF
      globalBoundaryArea=ABS(globalBoundaryArea)
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          globalBoundaryMeanPressure(boundaryID)=globalBoundaryPressure(boundaryID)/globalBoundaryArea(boundaryID)
        END IF
      END DO
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          computationalNode = ComputationalEnvironment_NodeNumberGet(ERR,ERROR)
          IF(computationalNode==0) THEN
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  flow:  ", &
              & globalBoundaryFlux(boundaryID),err,error,*999)
            CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  mean pressure:  ", &
              & globalBoundaryMeanPressure(boundaryID),err,error,*999)
            IF (toleranceCourant > ZERO_TOLERANCE) THEN
              CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Max Courant (CFL) number: ",maxCourant,err,error,*999)
            END IF
          END IF
        ELSE
          localError="Zero or negative area boundary detected on boundary "// &
            & TRIM(NUMBER_TO_VSTRING(boundaryID,"*",ERR,ERROR))//"."
          CALL FlagError(localError,ERR,ERROR,*999)
        END IF
      END DO

    CASE DEFAULT
      localError="Boundary flux calcluation for equations type "//TRIM(NUMBER_TO_VSTRING(equationsSet%Specification(3),"*", &
        & ERR,ERROR))//" is not yet implemented for Navier-Stokes."
      CALL FlagError(localError,ERR,ERROR,*999)
    END SELECT

    ! C o p y    i n t e g r a t e d   v a l u e s    t o    t a r g e t    f i e l d s
    ! ------------------------------------------------------------------------------------
    convergedFlag = .TRUE.
    ! Loop over elements again to allocate flux terms to boundary nodes
    DO elementIdx=1,elementsMapping3D%TOTAL_NUMBER_OF_LOCAL
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,5,elementNormal(1),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,6,elementNormal(2),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,7,elementNormal(3),err,error,*999)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,8,boundaryValue,err,error,*999)
      boundaryID=NINT(boundaryValue)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,9,boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      CALL Field_ParameterSetGetLocalElement(equationsSetField3D,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
       & elementIdx,11,boundaryValue,err,error,*999)
      coupledNodeNumber=NINT(boundaryValue)
      IF(boundaryID>1) THEN
        meshComponentNumber=2
        decompElement=>decomposition3D%TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)
        dependentBasis2=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%ELEMENTS% &
          & ELEMENTS(elementIdx)%BASIS

        ! B o u n d a r y   F a c e    N o r m a l s
        ! --------------------------------------------------
        DO faceIdx=1,dependentBasis2%NUMBER_OF_LOCAL_FACES
          !Get the face normal and quadrature information
          IF(ALLOCATED(decompElement%ELEMENT_FACES)) THEN
            faceNumber=decompElement%ELEMENT_FACES(faceIdx)
          ELSE
            CALL FlagError("Decomposition element faces is not allocated.",err,error,*999)
          END IF
          face=>decomposition3D%TOPOLOGY%FACES%FACES(faceNumber)
          IF(.NOT.(face%BOUNDARY_FACE)) CYCLE

          ! TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
          xiDirection = 0.0_DP
          SELECT CASE(dependentBasis%TYPE)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            xiDirection(3)=ABS(face%XI_DIRECTION)
          CASE(BASIS_SIMPLEX_TYPE)
            CALL FLAG_WARNING("Boundary flux calculation not yet set up for simplex element types.",ERR,ERROR,*999)
          CASE DEFAULT
            localError="Face integration for basis type "//TRIM(NUMBER_TO_VSTRING(dependentBasis%TYPE,"*",ERR,ERROR))// &
              & " is not yet implemented for Navier-Stokes."
            CALL FlagError(localError,ERR,ERROR,*999)
          END SELECT

          faceBasis=>decomposition3D%DOMAIN(meshComponentNumber)%PTR%TOPOLOGY%FACES%FACES(faceNumber)%BASIS
          faceQuadratureScheme=>faceBasis%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          !Use the geometric field to find the face normal and Jacobian for the face integral
          geometricInterpolationParameters=>equations3D%INTERPOLATION%geometricInterpParameters( &
            & FIELD_U_VARIABLE_TYPE)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_FACE_GET(FIELD_VALUES_SET_TYPE,faceNumber, &
            & geometricInterpolationParameters,err,error,*999)
          geometricInterpolatedPoint=>equations3D%INTERPOLATION%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR

          xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
          xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
          orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*face%XI_DIRECTION)
          pointMetrics=>equations3D%INTERPOLATION%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR
          CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,1, &
            & geometricInterpolatedPoint,err,error,*999)
          CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(COORDINATE_JACOBIAN_AREA_TYPE,pointMetrics,err,error,*999)
          CALL CrossProduct(pointMetrics%DX_DXI(:,1),pointMetrics%DX_DXI(:,2),faceNormal,ERR,ERROR,*999)
          faceNormal = faceNormal*orientation
          CALL Normalise(faceNormal,unitNormal,err,error,*999)
          CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
          normalTolerance=0.1_DP
          IF(normalDifference>normalTolerance) CYCLE

          ! U p d a t e    N o d a l   V a l u e s
          ! --------------------------------------------------
          ! Update local nodes with integrated boundary flow values
          DO faceNodeIdx=1,faceBasis%NUMBER_OF_NODES
            elementNodeIdx=dependentBasis2%NODE_NUMBERS_IN_LOCAL_FACE(faceNodeIdx,faceIdx)
            DO faceNodeDerivativeIdx=1,faceBasis%NUMBER_OF_DERIVATIVES(faceNodeIdx)
              nodeNumber=decomposition3D%DOMAIN(meshComponentNumber)%PTR% &
               & TOPOLOGY%ELEMENTS%ELEMENTS(elementIdx)%ELEMENT_NODES(elementNodeIdx)
              versionNumber=1
              IF(boundaryType==BOUNDARY_CONDITION_FIXED_CELLML) THEN
                ! Check current values against those passed to the CellML solver
                CALL Field_ParameterSetGetLocalNode(equationsSetField3D,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1,nodeNumber,1,q0D,err,error,*999)
                CALL Field_ParameterSetGetLocalNode(dependentField3D,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PRESSURE_VALUES_SET_TYPE,1,1,nodeNumber,4,p0D,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(dependentField3D,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionNumber,faceNodeDerivativeIdx,nodeNumber,4, &
                  & p0D,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(equationsSetField3D,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,versionNumber,faceNodeDerivativeIdx,nodeNumber,1, &
                  & q0D,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_NODE(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,1,globalBoundaryFlux(boundaryID),err,error,*999)
              END IF
            END DO !nodeDerivativeIdx
          END DO !faceNodeIdx
        END DO !faceIdx
      END IF !boundaryIdentifier
    END DO !elementIdx

    !allocate array for mpi communication
    IF(numberOfComputationalNodes>1) THEN !use mpi
      ALLOCATE(globalConverged(numberOfComputationalNodes),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",ERR,ERROR,*999)
      CALL MPI_ALLGATHER(convergedFlag,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL, &
       & computationalEnvironment%mpiCommunicator,MPI_IERROR)
      CALL MPI_ERROR_CHECK("MPI_ALLGATHER",MPI_IERROR,ERR,ERROR,*999)
      IF(ALL(globalConverged)) THEN
        convergedFlag = .TRUE.
      ELSE
        convergedFlag = .FALSE.
      END IF
    END IF

   ! Distribute any updated fields
    IF (ASSOCIATED(equationsSetField3D)) THEN
      CALL Field_ParameterSetUpdateStart(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(equationsSetField3D,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
    END IF
    IF (ASSOCIATED(dependentField3D)) THEN
      CALL Field_ParameterSetUpdateStart(dependentField3D,FIELD_U_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField3D,FIELD_U_VARIABLE_TYPE,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
    END IF

    EXITS("NavierStokes_CalculateBoundaryFlux3D0D")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateBoundaryFlux3D0D",err,error)
    RETURN 1
  END SUBROUTINE

  !
  !================================================================================================================================
  !

  !>Calculates the wall shear stress for fluid problems at all the boundary nodes
  SUBROUTINE NavierStokes_WallShearStressCalculate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<The equations set to calculate the wall shear stress for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryIdx,coordinateIdx,elementNumber,faceIdx,lineIdx,localNodeIdx,meshComponentNumber,nodeNumber, &
      & nodeIdx,numberOfBoundaries,numberOfDimensions,numberOfXi,startNode,stopNode,xiIdx
    REAL(DP) :: boundaryXi(2),deludelxi(3),fullXi(3),gradu,mu,normal(3),position(3),tangents(3,2),wss
    LOGICAL :: found
    TYPE(BASIS_TYPE), POINTER :: boundaryBasis,velocityBasis
    TYPE(DECOMPOSITION_TYPE), POINTER :: decomposition
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: domainElements
    TYPE(DOMAIN_FACES_TYPE), POINTER :: domainFaces
    TYPE(DOMAIN_FACE_TYPE), POINTER :: face
    TYPE(DOMAIN_LINES_TYPE), POINTER :: domainLines
    TYPE(DOMAIN_LINE_TYPE), POINTER :: line
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: nodesMappings
    TYPE(DOMAIN_MAPPINGS_TYPE), POINTER :: domainMappings
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: dependentInterpolationParameters,geometricInterpolationParameters, &
      & materialsInterpolationParameters
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: dependentInterpolatedPoint,geometricInterpolatedPoint, &
      & materialsInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: geometricInterpPointMetrics
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_WallShearStressCalculate",err,error,*999)


    ! Preliminary checks; get field and domain pointers
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    SELECT CASE(equationsSet%specification(3))
    CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      IF(ASSOCIATED(dependentField%VARIABLE_TYPE_MAP(FIELD_W_VARIABLE_TYPE)%ptr)) THEN
        NULLIFY(geometricField)
        NULLIFY(materialsField)
        NULLIFY(equations)
        NULLIFY(vectorEquations)
        NULLIFY(vectorMapping)
        NULLIFY(nonlinearMapping)
        NULLIFY(dependentVariable)
        NULLIFY(decomposition)
        NULLIFY(decompositionTopology)
        NULLIFY(domain)
        NULLIFY(domainMappings)
        NULLIFY(domainTopology)
        NULLIFY(domainElements)
        NULLIFY(nodesMappings)
        NULLIFY(domainNodes)
        NULLIFY(domainLines)
        NULLIFY(domainFaces)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
        CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,1,1,dependentVariable,err,error,*999)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        meshComponentNumber=dependentVariable%COMPONENTS(1)%MESH_COMPONENT_NUMBER
        CALL Decomposition_TopologyGet(decomposition,decompositionTopology,err,error,*999)
        CALL Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*999)
        CALL Domain_MappingsGet(domain,domainMappings,err,error,*999)
        CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
        CALL DomainTopology_ElementsGet(domainTopology,domainElements,err,error,*999)
        CALL DomainMappings_NodesGet(domainMappings,nodesMappings,err,error,*999)
        CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
        numberOfDimensions=dependentVariable%NUMBER_OF_COMPONENTS-1

        IF(numberOfDimensions==2) THEN
          IF(.NOT.decomposition%CALCULATE_LINES) CALL FlagError("Decomposition calculate lines is not set.",err,error,*999)
          CALL DomainTopology_LinesGet(domainTopology,domainLines,err,error,*999)
        ELSE IF(numberOfDimensions==3) THEN
          IF(.NOT.decomposition%CALCULATE_FACES) CALL FlagError("Decomposition calculate faces is not set.",err,error,*999)
          CALL DomainTopology_FacesGet(domainTopology,domainFaces,err,error,*999)
        ELSE
          localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        ENDIF

        dependentInterpolationParameters=>equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        geometricInterpolationParameters=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        materialsInterpolationParameters=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
        dependentInterpolatedPoint=>equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        materialsInterpolatedPoint=>equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
        geometricInterpPointMetrics=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr

        !Loop over internal and boundary nodes
        startNode = nodesMappings%INTERNAL_START
        stopNode = nodesMappings%BOUNDARY_FINISH
        !Loop over internal and boundary nodes
        nodes: DO nodeIdx=startNode,stopNode
          nodeNumber=nodesMappings%DOMAIN_LIST(nodeIdx)
          IF(.NOT.domainNodes%nodes(nodeNumber)%BOUNDARY_NODE) CYCLE nodes
          !Node is on the boundary. Loop over surrounding faces/lines
          IF(numberOfDimensions == 2) THEN
            numberOfBoundaries = domainNodes%nodes(nodeNumber)%NUMBER_OF_NODE_LINES
          ELSE
            numberOfBoundaries = domainNodes%nodes(nodeNumber)%NUMBER_OF_NODE_FACES
          ENDIF
          IF(numberOfBoundaries<=0) THEN
            localError="The number of boundaries for node "//TRIM(NumberToVString(nodeNumber,"*",err,error))// &
              & " of "//TRIM(NumberToVString(numberOfBoundaries,"*",err,error))// &
              & " is invalid. The number of boundaries must be > 0."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          wss=0.0_DP
          DO boundaryIdx=1,numberOfBoundaries
            boundaryXi=0.0_DP
            IF(numberOfDimensions==2) THEN
              lineIdx=domainNodes%nodes(nodeNumber)%NODE_LINES(boundaryIdx)
              NULLIFY(line)
              CALL DomainLines_LineGet(domainLines,lineIdx,line,err,error,*999)
              boundaryBasis=>line%basis
              IF(.NOT.ASSOCIATED(boundaryBasis)) THEN
                localError="Boundary basis is not associated for line number "// &
                  & TRIM(NumberToVString(lineIdx,"*",err,error))// "."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              elementNumber=line%ELEMENT_NUMBER
              NULLIFY(velocityBasis)
              CALL DomainElements_BasisGet(domainElements,elementNumber,velocityBasis,err,error,*999)
              !Find node position
              found=.FALSE.
              DO localNodeIdx=1,boundaryBasis%NUMBER_OF_NODES
                IF(line%NODES_IN_LINE(localNodeIdx)==nodeNumber) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !localNodeIdX
              IF(.NOT.found) THEN
                localError="Could not find node number "//TRIM(NumberToVString(nodeNumber,"*",err,error))// &
                  & " in line number "//TRIM(NumberToVString(lineIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
                & err,error,*999)
            ELSE
              faceIdx=domainNodes%nodes(nodeNumber)%NODE_FACES(boundaryIdx)
              NULLIFY(face)
              CALL DomainFaces_FaceGet(domainFaces,faceIdx,face,err,error,*999)
              boundaryBasis=>face%basis
              IF(.NOT.ASSOCIATED(boundaryBasis)) THEN
                localError="Boundary basis is not associated for face number "// &
                  & TRIM(NumberToVString(faceIdx,"*",err,error))// "."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              elementNumber=face%ELEMENT_NUMBER
              NULLIFY(velocityBasis)
              CALL DomainElements_BasisGet(domainElements,elementNumber,velocityBasis,err,error,*999)
              !Find node position
              found=.FALSE.
              DO localNodeIdx=1,boundaryBasis%NUMBER_OF_NODES
                IF(face%NODES_IN_FACE(localNodeIdx)==nodeNumber) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !localNodeIdX
              IF(.NOT.found) THEN
                localError="Could not find node number "//TRIM(NumberToVString(nodeNumber,"*",err,error))// &
                  & " in face number "//TRIM(NumberToVString(faceIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
                & err,error,*999)
            ENDIF
            CALL Field_InterpolateXi(FIRST_PART_DERIV,boundaryXi,geometricInterpolatedPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(boundaryBasis%NUMBER_OF_XI,geometricInterpPointMetrics,err,error,*999)
            CALL Field_PositionNormalTangentsCalculateIntPtMetric(geometricInterpPointMetrics,.FALSE.,position,normal,tangents, &
              & err,error,*999)
            numberOfXi=velocityBasis%NUMBER_OF_XI
            CALL Basis_LocalNodeXiCalculate(boundaryBasis,localNodeIdx,boundaryXi,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpolationParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpolationParameters, &
              & err,error,*999)
            CALL Basis_BoundaryXiToXi(velocityBasis,boundaryIdx,boundaryXi,fullXi,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,fullXi,geometricInterpolatedPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(velocityBasis%NUMBER_OF_XI,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,fullXi,dependentInterpolatedPoint,err,error,*999)
            CALL Field_InterpolateXi(NO_PART_DERIV,fullXi,materialsInterpolatedPoint,err,error,*999)
            mu=materialsInterpolatedPoint%VALUES(1,NO_PART_DERIV)
            DO coordinateIdx=1,numberOfDimensions
              DO xiIdx=1,numberOfXi
                deludelxi(xiIdx)=dependentInterpolatedPoint%values(coordinateIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
              ENDDO !xiIdx
              gradu=DOT_PRODUCT(deludelxi(1:numberOfXi),geometricInterpPointMetrics%DXI_DX(1:numberOfXi,coordinateIdx))
              wss=wss+mu*gradu*normal(coordinateIdx)
            ENDDO !coordinateIdx
          ENDDO !boundaryIdx
          wss=wss/REAL(numberOfBoundaries)
          CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_W_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,nodeNumber,1, &
            & wss,err,error,*999)
        ENDDO nodes !nodeIdx
      ENDIF
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
       & " is not valid for a Navier-Stokes equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_WallShearStressCalculate")
    RETURN
999 ERRORSEXITS("NavierStokes_WallShearStressCalculate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_WallShearStressCalculate

  !
  !================================================================================================================================
  !

END MODULE NAVIER_STOKES_EQUATIONS_ROUTINES
