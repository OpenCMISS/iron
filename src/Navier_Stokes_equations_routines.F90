!> \file  
!> \author Sebastian Krittian
!> \brief This module handles all Navier-Stokes fluid routines.
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
!> Contributor(s): Sebastian Krittian, David Ladd, Soroush Safaei, Chris Bradley
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

!>This module handles all Navier-Stokes fluid routines.
MODULE NavierStokesEquationsRoutines

  USE AdvectionEquationsRoutines
  USE AnalyticAnalysisRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE CharacteristicEquationsRoutines
  USE CmissMPI 
  USE CmissPetsc
  USE CmissPetscTypes
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemAccessRoutines
  USE DecompositionRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FLUID_MECHANICS_IO_ROUTINES
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lapack
  USE Maths
  USE MatrixVector
#ifndef NOMPIMOD
  USE MPI
#endif
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE SolverMatricesAccessRoutines
  USE StreeEquationsRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PUBLIC NavierStokes_AnalyticFunctionsEvaluate

  PUBLIC NavierStokes_EquationsSetSpecificationSet

  PUBLIC NavierStokes_EquationsSetSolutionMethodSet

  PUBLIC NavierStokes_EquationsSetSetup

  PUBLIC NavierStokes_PreSolveALEUpdateParameters,NavierStokes_PreSolveUpdateBoundaryConditions, &
    & NavierStokes_PreSolveALEUpdateMesh

  PUBLIC NavierStokes_PreSolve,NavierStokes_PostSolve

  PUBLIC NavierStokes_ProblemSpecificationSet

  PUBLIC NavierStokes_ProblemSetup

  PUBLIC NavierStokes_FiniteElementResidualEvaluate,NavierStokes_FiniteElementJacobianEvaluate

  PUBLIC NavierStokes_BoundaryConditionsAnalyticCalculate

  PUBLIC NavierStokes_ResidualBasedStabilisation

  PUBLIC NavierStokes_Couple1D0D

  PUBLIC NavierStokes_CoupleCharacteristics

  PUBLIC NavierStokes_FiniteElementPreResidualEvaluate

  PUBLIC NavierStokes_PostLoop

  PUBLIC NavierStokes_UpdateMultiscaleBoundary

  PUBLIC NavierStokes_WallShearStressCalculate

CONTAINS

!
!================================================================================================================================
!

  !>Sets/changes the solution method for a Navier-Stokes flow equation type of an fluid mechanics equations set class.
  SUBROUTINE NavierStokes_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
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
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(solutionMethod)
      CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
      CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
        equationsSet%solutionMethod=EQUATIONS_SET_NODAL_SOLUTION_METHOD
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
        localError="The specified solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes flow equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("NavierStokes_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already associated.",err,error,*999)
    IF(SIZE(specification,1)<3) &
      & CALL FlagError("Equations set specification must have three entries for a Navier-Stokes type equations set.", &
      & err,error,*999)
    
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
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE)
      CALL FlagError("Not implemented yet.",err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Set full specification
    ALLOCATE(equationsSet%specification(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE,subtype]

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
  SUBROUTINE NavierStokes_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,constantBasedComponents,elementBasedComponents,esSpecification(3),geometricComponentNumber, &
      & geometricMeshComponent,geometricScalingType,nodeBasedComponents,numberOfAnalyticComponents,numberOfDependentComponents, &
      & numberOfDimensions,numberOfIndependentComponents,numberOfMaterialsComponents1,numberOfMaterialsComponents2, &
      & solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: geometricField,materialsField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
       
    SELECT CASE(esSpecification(3))
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
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(equationsSetSetup%setupType)
      CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        !-----------------------------------------------------------------
        ! I n i t i a l   s e t u p
        !-----------------------------------------------------------------
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL NavierStokes_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
            CALL EquationsSet_LabelSet(equationsSet,"Navier-Stokes equations set",err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is not implemented for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL NavierStokes_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
            CALL EquationsSet_LabelSet(equationsSet,"Navier-Stokes equations set",err,error,*999)
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              !Create the auto created equations set field field for SUPG element metrics
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField,err,error,*999)
              CALL Field_LabelSet(equationsField%equationsSetField,"Equations Set Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsField%equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesSet(equationsField%equationsSetField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsField%equationsSetField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,"Penalty Coefficient", &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
              !Default the penalty coefficient value to 1E4
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 1,1.0E4_DP,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is not implemented for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)          
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          nodeBasedComponents = 2  ! boundary flow, pressure
          elementBasedComponents = 11 ! 4 element metrics, 3 boundary normal components, boundaryID, boundaryType, 
          ! C1, coupledNodeNumber
          constantBasedComponents = 4 ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL NavierStokes_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
            CALL EquationsSet_LabelSet(equationsSet,"Navier-Stokes equations set",err,error,*999)
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              !Create the auto created equations set field field for SUPG element metrics
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField,err,error,*999)
              CALL Field_LabelSet(equationsField%equationsSetField,"Equations Set Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsField%equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesSet(equationsField%equationsSetField,3,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsField%equationsSetField,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                & FIELD_U1_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,"BoundaryFlow",err,error,*999)
              CALL Field_VariableLabelSet(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE,"ElementMetrics",err,error,*999)
              CALL Field_VariableLabelSet(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE,"EquationsConstants", &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
                & nodeBasedComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE, &
                & elementBasedComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                & constantBasedComponents,err,error,*999)
            ELSE
              localError="User-specified fields are not yet implemented for an equations set field field setup type of "// &
                & TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))//" for a Navier-Stokes fluid."
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
              !Default the Element Metrics parameter values 0.0
              !Init boundary flux to 0
              DO componentIdx=1,nodeBasedComponents
                CALL Field_ComponentValuesInitialise(equationsField%equationsSetField, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
              ENDDO !componentIdx
              !Init Element Metrics to 0 (except C1)
              DO componentIdx=1,elementBasedComponents-1
                CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
              ENDDO !componentIdx
              !Default C1 to -1 for now, will be calculated in ResidualBasedStabilisation if not specified by user
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,elementBasedComponents,-1.0_DP,err,error,*999)
              !Boundary stabilisation scale factor (beta): default to 0
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
              !Max Courant (CFL) number: default to 0.0 (do not use)
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
              !Init Time increment to 0
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
              !Stabilisation type: default to 1 for RBS (0=none, 2=RBVM)
              CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)
            ELSE
              localError="User-specified fields are not yet implemented for an equations set field field setup type of "// &
                & TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))//" for a Navier-Stokes fluid."
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is not implemented for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        !-----------------------------------------------------------------
        ! G e o m e t r i c   f i e l d
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsField%equationsSetField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsField%equationsSetField,geometricField,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsField%equationsSetField,geometricScalingType,err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            ! do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            nodeBasedComponents = 2 !boundary flow, pressure
            elementBasedComponents = 11 !4 element metrics, 3 boundary normal components, boundaryID, boundaryType,
            !C1, coupledNodeNumber
            constantBasedComponents = 4 ! maxCFL, boundaryStabilisationBeta, timeIncrement, stabilisationType
            IF(equationsField%equationsSetFieldAutoCreated) THEN
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsField%equationsSetField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsField%equationsSetField,geometricField,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              !Node-based fields
              DO componentIdx=1,nodeBasedComponents
                CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Element-based fields
              DO componentIdx=1,elementBasedComponents
                CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Constant fields: boundary stabilisation scale factor and max courant #
              DO componentIdx=1,constantBasedComponents
                CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                  & componentIdx,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_U1_VARIABLE_TYPE, &
                  & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsField%equationsSetField,geometricScalingType,err,error,*999)
            ELSE
              !Do nothing
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
              CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,3,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_W_VARIABLE_TYPE],err,error,*999)
              ELSE
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              ENDIF
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              !Calculate number of components with one component for each dimension and one for pressure
              numberOfDependentComponents=numberOfDimensions+1
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              DO componentIdx=1,numberOfDependentComponents
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE,"Wss",err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE,FIELD_DP_TYPE, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE,1, &
                  & err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE,1, &
                  & geometricMeshComponent,err,error,*999)
              ENDIF
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfDependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_W_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDIF
                !Default geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,3,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_W_VARIABLE_TYPE],err,error,*999)
              ELSE
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                  & err,error,*999)
              ENDIF
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension and one for pressure
              numberOfDependentComponents=numberOfDimensions+1
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
                & err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_W_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_W_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_W_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfDependentComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                  & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_W_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_ParameterSetEnsureCreated(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
                DO componentIdx=1,numberOfDependentComponents
                  CALL Field_ComponentValuesInitialise(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !calculate number of components (Q,A) for U and dUdN
            numberOfDependentComponents=2
            IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              !start field creation with name 'DEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
              !define new created field to be dependent
              CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
              !set number of variables to 6 (U,DELUDELN,V,U1,U2)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,5,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              !Set data type
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              !Calculate number of components
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              !2 component (W1,W2) for V
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                & geometricMeshComponent,err,error,*999)
              !Default to the geometric interpolation setup
              DO componentIdx=1,numberOfDependentComponents
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfDependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U2_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field- Characteristic equations
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,5,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE, FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              !Calculate number of components (Q,A) for U and dUdN
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
                & err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              !2 component (W1,W2) for V
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%dependent%dependentFieldAutoCreated) &
              & CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)              
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              !start field creation with name 'DEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
              !define new created field to be dependent
              CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
              !set number of variables to 5 (U,DELUDELN,V,U1,U2)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,5,err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_U3_VARIABLE_TYPE],err,error,*999)
              ELSE
                CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                  & err,error,*999)
              ENDIF
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              !calculate number of components (Q,A)
              numberOfDependentComponents=2
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                & geometricMeshComponent,err,error,*999)
              !Default to the geometric interpolation setup
              DO componentIdx=1,numberOfDependentComponents
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U1_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U2_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfDependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U1_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                    & FIELD_U2_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%dependent%dependentField,geometricScalingType,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field- Characteristic equations
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,5,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components (Q,A) for U and dUdN
              numberOfDependentComponents=2
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE, &
                & numberOfDependentComponents,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U2_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%dependent%dependentFieldAutoCreated) &
              & CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsIndependent)
          CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              !Create the auto created independent field
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
              CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
              CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                & err,error,*999)
              !Calculate number of components with one component for each dimension
              numberOfIndependentComponents=numberOfDimensions
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfIndependentComponents,err,error,*999)
              !Default to the geometric interpolation setup
              DO componentIdx=1,numberOfIndependentComponents
                CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                  & err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
              END DO !componentIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                !Default geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
                !Other solutions not defined yet
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              !calculate number of components with one component for each dimension
              numberOfIndependentComponents=numberOfDimensions
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfIndependentComponents, &
                & err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a standard Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsIndependent)
          CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !normalDirection for wave relative to node for W1,W2
            numberOfIndependentComponents=2
            !Create the auto created independent field
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              !start field creation with name 'independentField'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
              !define new created field to be independent
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              !set number of variables to 1
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, FIELD_DP_TYPE, &
                & err,error,*999)
              !calculate number of components with one component for each dimension
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfIndependentComponents,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              !Default to the geometric interpolation setup
              DO componentIdx=1,numberOfIndependentComponents
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO !componentIdx
                CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
              CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO !componentIdx
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field- Characteristic equation
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfIndependentComponents, &
                & err,error,*999)
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) &
              & CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a standard Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsIndependent)
          CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !normalDirection for wave relative to node for W1,W2
            numberOfIndependentComponents=2
            !Create the auto created independent field
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              ! Do nothing? independent field should be set up by characteristic equation routines
            ELSE
              !Check the user specified field- Characteristic equation
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfIndependentComponents, &
                & err,error,*999)
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) &
              & CALL Field_CreateFinish(equationsSet%INDEPENDENT%independentField,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a standard Navier-Stokes fluid"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
          NULLIFY(equationsIndependent)
          CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              !Create the auto created independent field
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
              CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
              CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              !calculate number of components with one component for each dimension
              numberOfIndependentComponents=numberOfDimensions
              CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & numberOfIndependentComponents,err,error,*999)
              !Default to the geometric interpolation setup
              DO componentIdx=1,numberOfIndependentComponents
                CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                  & err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO !componentIdx
                !Default geometric field scaling
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
                !Other solutions not defined yet
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              numberOfIndependentComponents=numberOfDimensions
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                & numberOfIndependentComponents,err,error,*999)
              CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfIndependentComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO !componentIdx
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsIndependent%independentFieldAutoCreated) THEN
              CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a standard Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
        !-----------------------------------------------------------------
        ! A n a l y t i c   t y p e
        !-----------------------------------------------------------------
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(esSpecification(3))
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
          NULLIFY(equationsAnalytic)
          CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
          SELECT CASE(equationsSetSetup%actionType)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
            CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
            CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE
              !Check that domain is 2D
              IF(numberOfDimensions/=2) THEN
                localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                  & " is invalid. The analytic function type of "// &
                  & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                  & " requires that there be 2 geometric dimensions."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              !Check the materials values are constant
              NULLIFY(materialsField)
              CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
              CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
                & err,error,*999)
              CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION, &
                & err,error,*999)
              !Set analytic function type
              equationsAnalytic%analyticFunctionType=equationsSetSetup%analyticFunctionType
              numberOfAnalyticComponents=4
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
              & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
              & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
              & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
              !Check that this is a 1D equations set
              IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
                & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
                !Set analytic function type
                equationsAnalytic%analyticFunctionType=equationsSetSetup%analyticFunctionType
                !Set numbrer of components- Q,A (same as N-S depenedent field)
                numberOfAnalyticComponents=2
              ELSE
                localError="The third equations set specification must by a TRANSIENT1D or COUPLED1D0D "// &
                  & "to use an analytic function of type "// &
                  & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
              !Check that domain is 2D/3D
              IF(numberOfDimensions<2 .OR. numberOfDimensions>3) THEN
                localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                  & " is invalid. The analytic function type of "// &
                  & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                  & " requires that there be 2 or 3 geometric dimensions."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              !Set analytic function type
              equationsAnalytic%analyticFunctionType=equationsSetSetup%analyticFunctionType
              !Set numbrer of components
              numberOfAnalyticComponents=10
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN
              !Check that domain is 2D
              IF(numberOfDimensions/=2) THEN
                localError="The number of geometric dimensions of "// &
                  & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                  & " is invalid. The analytic function type of "// &
                  & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                  & " requires that there be 2 geometric dimensions."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              !Check the materials values are constant
              NULLIFY(materialsField)
              CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
              CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
                & err,error,*999)
              CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION, &
                & err,error,*999)
              !Set analytic function type
              equationsAnalytic%analyticFunctionType=equationsSetSetup%analyticFunctionType
              numberOfAnalyticComponents=2
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
              !Set analtyic function type
              equationsAnalytic%analyticFunctionType=EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1
            CASE DEFAULT
              localError="The specified analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for an analytic Navier-Stokes problem."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Create analytic field if required
            IF(numberOfAnalyticComponents>=1) THEN
              IF(equationsAnalytic%analyticFieldAutoCreated) THEN
                !Create the auto created analytic field
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsAnalytic%analyticField,err,error,*999)
                CALL Field_LabelSet(equationsAnalytic%analyticField,"Analytic Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsAnalytic%analyticField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(equationsAnalytic%analyticField,FIELD_INDEPENDENT_TYPE,err,error,*999)
                NULLIFY(geometricDecomposition)
                CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsAnalytic%analyticField,geometricDecomposition,err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsAnalytic%analyticField,geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsAnalytic%analyticField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsAnalytic%analyticField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,"Analytic",err,error,*999)
                CALL Field_DimensionSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                !Set the number of analytic components
                CALL Field_NumberOfComponentsSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfAnalyticComponents,err,error,*999)
                !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
                DO componentIdx=1,numberOfAnalyticComponents
                  CALL Field_ComponentMeshComponentSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & geometricMeshComponent,err,error,*999)
                  IF(equationsSetSetup%analyticFunctionType== EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                    CALL Field_ComponentInterpolationSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ELSE
                    CALL Field_ComponentInterpolationSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  ENDIF
                ENDDO !componentIdx
                !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsAnalytic%analyticField,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                IF(numberOfAnalyticComponents==1) THEN
                  CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                ENDIF
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                  & numberOfAnalyticComponents,err,error,*999)
              ENDIF
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
            IF(equationsAnalytic%analyticFieldAutoCreated) THEN
              !Finish creating the analytic field
              CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
              !Set the default values for the analytic field
              SELECT CASE(esSpecification(3))
              CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
                SELECT CASE(equationsAnalytic%analyticFunctionType)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
                  !Default the analytic parameter values (L, H, U_mean, Pout) to 0.0
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,4,0.0_DP,err,error,*999)
                CASE DEFAULT
                  localError="The analytic function type of "// &
                    & TRIM(NumberToVString(equationsAnalytic%analyticFunctionType,"*",err,error))// &
                    & " is invalid for an analytical static Navier-Stokes equation."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                SELECT CASE(equationsAnalytic%analyticFunctionType)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
                  !Default the analytic parameter values (U_characteristic, L) to 0.0
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
                  !Default the analytic parameter values to 0
                  numberOfAnalyticComponents = 10
                  DO componentIdx = 1,numberOfAnalyticComponents
                    CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                  ENDDO !componentIdx
                CASE DEFAULT
                  localError="The analytic function type of "// &
                    & TRIM(NumberToVString(equationsAnalytic%analyticFunctionType,"*",err,error))// &
                    & " is invalid for an analytical transient Navier-Stokes equation."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
                SELECT CASE(equationsAnalytic%analyticFunctionType)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                  & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN, &
                  & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART, &
                  & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                  !Default the analytic parameter period values to 0
                  CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                CASE DEFAULT
                  localError="The analytic function type of "// &
                    & TRIM(NumberToVString(equationsAnalytic%analyticFunctionType,"*",err,error))// &
                    & " is invalid for a 1D Navier-Stokes equation."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The third equations set specification of "// &
                  & TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                  & " is invalid for an analytical Navier-Stokes equation set."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for an analytic Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The third equations set specification of "// &
            & TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CALL EquationsSet_AssertMaterialsIsCreated(equationsSet,err,error,*999)
        NULLIFY(equationsMaterials)
        CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfCOmponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
          numberOfMaterialsComponents1=2 !viscosity, density
          SELECT CASE(equationsSetSetup%actionType)
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsMaterials%materialsFieldAutoCreated) THEN
              !Create the auto created materials field
              !start field creation with name 'MATERIAL_FIELD'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
              CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,1,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
              CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & numberOfMaterialsComponents1,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & 1,geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
            ENDIF
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsMaterials%materialsFieldAutoCreated) THEN
              !Finish creating the materials field
              CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
              !Set the default values for the materials field
              ! viscosity,density=1
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          numberOfMaterialsComponents1=2 !U_var (constant)  : viscosity scale, density
          numberOfMaterialsComponents2=2 !V_var (gaussBased): viscosity, shear rate
          SELECT CASE(equationsSetSetup%actionType)
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsMaterials%materialsFieldAutoCreated) THEN
              !Create the auto created materials field
              !start field creation with name 'MATERIAL_FIELD'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
              CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsMaterials%materialsField,"MaterialsField",err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE],err,error,*999)
              !Set up U_VARIABLE (constants)
              CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"MaterialsConstants", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & numberOfMaterialsComponents1,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              DO componentIdx=1,numberOfMaterialsComponents1
                CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Set up V_VARIABLE (gauss-point based, CellML in/out parameters)
              CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,"ConstitutiveValues", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & numberOfMaterialsComponents2,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              DO componentIdx=1,numberOfMaterialsComponents2
                CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,geometricComponentNumber,err,error,*999)
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                  & componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              !Check the U_VARIABLE
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents1, &
                & err,error,*999)
              !Check the U_VARIABLE
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
                & numberOfMaterialsComponents2,err,error,*999)
            ENDIF
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsMaterials%materialsFieldAutoCreated) THEN
              !Finish creating the materials field
              CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
              !Set the default values for the materials constants (viscosity scale, density)
              DO componentIdx=1,numberOfMaterialsComponents2
                CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
              ENDDO !componentIdx
              !Set the default values for the materials consitutive parameters (viscosity scale, density)
              DO componentIdx=1,numberOfMaterialsComponents2
                CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
              ENDDO !componentIdx
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          ! 1 variables for the 1D Navier-Stokes materials
          numberOfMaterialsComponents1=8 !(MU,RHO,alpha,pressureExternal,LengthScale,TimeScale,MassScale)
          numberOfMaterialsComponents2=3 !(A0,E,H0)
          SELECT CASE(equationsSetSetup%actionType)
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsMaterials%materialsFieldAutoCreated) THEN
              !Create the auto created materials field
              !start field creation with name 'MATERIAL_FIELD'
              CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
              CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
              NULLIFY(geometricDecomposition)
              CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
              CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,2,err,error,*999)
              !2 U,V materials field
              CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField, &
                & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              !Set up Navier-Stokes materials parameters
              CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & numberOfMaterialsComponents1,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & numberOfMaterialsComponents2,err,error,*999)
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,1, &
                & geometricComponentNumber,err,error,*999)
              DO componentIdx=1,numberOfMaterialsComponents1 
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              DO componentIdx=1,numberOfMaterialsComponents2
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              ! Set up coupling materials parameters
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
              !Default the field scaling to that of the geometric field
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              !Check N-S field variable
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfMaterialsComponents2, &
                & err,error,*999)
            ENDIF
            !Specify start action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Finish creating the materials field
            IF(equationsMaterials%materialsFieldAutoCreated) &
              & CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          !\todo: Think about gravity
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
            !? Maybe set finished flag????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            &  " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            &  " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            NULLIFY(equations)
            CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
            CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
            CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              NULLIFY(vectorMapping)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
              CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              NULLIFY(vectorMatrices)
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              !Use the analytic Jacobian calculation
              CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
                & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
            CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
              !Finish the creation of the equations
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              NULLIFY(vectorMapping)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
              CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
               CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              NULLIFY(vectorMatrices)
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(equations%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              !Use the analytic Jacobian calculation
              CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
                & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            NULLIFY(equations)
            CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
            CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
            CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              NULLIFY(vectorMapping)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
              CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              NULLIFY(vectorMatrices)
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE, &
                  & MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE, &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1, &
                  & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1, &
                  & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              !Use the analytic Jacobian calculation
              CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
                & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            NULLIFY(equations)
            CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
            CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
            CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              NULLIFY(vectorMapping)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
              CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              NULLIFY(vectorMatrices)
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              !Use the analytic Jacobian calculation
              CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
                & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " does not equal a Navier-Stokes fluid subtype."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("NavierStokes_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem pre solve.
  SUBROUTINE NavierStokes_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,esSpecification(3),iterationNumber,maxIterations,numberOfEquationsSets, &
      & numberOfSolverMatrices,pSpecification(3),solverGlobalNumber,solverMatrixIdx,solveType
    REAL(DP) :: absoluteTolerance,currentTime,relativeTolerance,timeIncrement
    LOGICAL :: continueLoop
    TYPE(ControlLoopType), POINTER :: controlLoop,parentLoop
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldParameterSetType), POINTER :: inputParameterSet,upwindParameterSet
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverMatricesType), POINTER :: solverMatrices
    TYPE(SolverMatrixType), POINTER :: solverMatrix
    TYPE(SolverType), POINTER :: solver2
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("NavierStokes_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    !Since we can have a fluid mechanics navier stokes equations set in a coupled problem setup we do not necessarily
    !have pSpecification(1)==FLUID_MECHANICS_CLASS
    SELECT CASE(pSpecification(1))
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)        
        ! TODO: Set up for multiple equations sets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          !Update boundary conditions and any analytic values
          CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
        ENDIF
      CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
        !Update transient boundary conditions and any analytic values
        CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
      CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        SELECT CASE(solveType)
        CASE(SOLVER_DYNAMIC_TYPE)
          ! --- D y n a m i c    S o l v e r s ---
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
              ! --- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
              IF(pSpecification(3) == PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
                NULLIFY(parentLoop)
                CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
                CALL ControlLoop_CurrentWhileInformationGet(parentLoop,iterationNumber,maxIterations,absoluteTolerance, &
                  & relativeTolerance,continueLoop,err,error,*999)
                IF(iterationNumber == 1) THEN
                  ! Only update fixed BCs once per timestep
                  CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
                ENDIF
              ELSE
                !Update boundary conditions and any analytic values
                CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              ! --- 1 D    N a v i e r - S t o k e s   E q u a t i o n s ---
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
                & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
                & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
              !Update boundary conditions and any analytic values
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
              ! --- A d v e c t i o n   S o l v e r ---
            CASE(EQUATIONS_SET_ADVECTION_SUBTYPE)
              CALL Advection_PreSolve(solver,err,error,*999)
            CASE DEFAULT
              localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*", &
                & err,error))//" is not valid for a nonlinear Navier-Stokes solver."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsSetIdx
        CASE(SOLVER_NONLINEAR_TYPE)
          ! --- N o n l i n e a r    S o l v e r s ---
          CALL ControlLoop_CurrentWhileInformationGet(controlLoop,iterationNumber,maxIterations,absoluteTolerance, &
            & relativeTolerance,continueLoop,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            SELECT CASE(esSpecification(3))
              ! --- C h a r a c t e r i s t i c   E q u a t i o n s ---
            CASE(EQUATIONS_SET_CHARACTERISTIC_SUBTYPE)
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,err,error,*999)
              CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
                & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,err,error,*999)
              IF(iterationNumber == 1) THEN
                CALL Field_ParameterSetEnsureCreated(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE, &
                  & err,error,*999)
                ! Extrapolate new W from Q,A if this is the first timestep
                ! (otherwise will be calculated based on Navier-Stokes values)
                CALL Characteristic_Extrapolate(solver,currentTime,timeIncrement,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                & " is not valid for a nonlinear Navier-Stokes solver."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsSetIdx
        CASE(SOLVER_DAE_TYPE)
          ! --- C e l l M L   S o l v e r ---
          CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime + timeIncrement,err,error,*999)
          CALL Solver_DAETimeStepSet(solver,timeIncrement/1000.0_DP,err,error,*999)
        CASE(SOLVER_LINEAR_TYPE)
          ! --- L i n e a r   S o l v e r s ---
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            SELECT CASE(esSpecification(3))
              ! --- S t r u c t u r e d   T r e e  E q u a t i o n s ---
            CASE(EQUATIONS_SET_STREE1D0D_SUBTYPE, &
              & EQUATIONS_SET_STREE1D0D_ADV_SUBTYPE)
              CALL Stree_PreSolve(SOLVER,err,error,*999)
            CASE DEFAULT
              localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                & " is not valid for a linear Navier-Stokes solver."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsSetIdx
        CASE DEFAULT
          localError="The solve type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
            & " is invalid for a multiscale Navier-Stokes problem type."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        &  PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        SELECT CASE(solveType)
          ! This switch takes advantage of the uniqueness of the solver types to do pre-solve operations
          ! for each of solvers in the various possible 1D subloops

          ! --- C h a r a c t e r i s t i c   S o l v e r ---
        CASE(SOLVER_NONLINEAR_TYPE)
          CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
          CALL ControlLoop_CurrentWhileInformationGet(controlLoop,iterationNumber,maxIterations,absoluteTolerance, &
            & relativeTolerance,continueLoop,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          ! Characteristic solver effectively solves for the mass/momentum conserving fluxes at the
          ! *NEXT* timestep by extrapolating current field values and then solving a system of nonlinear
          ! equations: cons mass, continuity of pressure, and the characteristics.
          NULLIFY(fieldVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,fieldVariable,err,error,*999)
          NULLIFY(inputParameterSet)
          CALL FieldVariable_ParameterSetCheck(fieldVariable,FIELD_INPUT_DATA1_SET_TYPE,inputParameterSet,err,error,*999)
          IF(.NOT.ASSOCIATED(inputParameterSet)) THEN
            CALL Field_ParameterSetCreate(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
          ENDIF
          CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & FIELD_INPUT_DATA1_SET_TYPE,1.0_DP,err,error,*999)
          CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_RESIDUAL_SET_TYPE, &
            & FIELD_INPUT_DATA2_SET_TYPE,1.0_DP,err,error,*999)

          IF(iterationNumber == 1) THEN
            NULLIFY(upwindParameterSet)
            CALL FieldVariable_ParameterSetCheck(fieldVariable,FIELD_UPWIND_VALUES_SET_TYPE,upwindParameterSet,err,error,*999)
            IF(.NOT.ASSOCIATED(upwindParameterSet)) &
              & CALL Field_ParameterSetCreate(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
            ! Extrapolate new W from Q,A if this is the first timestep (otherwise will be calculated based on Navier-Stokes
            ! values)
            CALL Characteristic_Extrapolate(solver,currentTime,timeIncrement,err,error,*999)
          ENDIF

          ! --- 1 D   N a v i e r - S t o k e s   S o l v e r ---
        CASE(SOLVER_DYNAMIC_TYPE)
          CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
          IF(solverGlobalNumber==2) THEN
            ! update solver matrix
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            NULLIFY(solverMapping)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            NULLIFY(solverMatrices)
            CALL SolverEquations_SolverMatricesGet(solverEquations,solverMatrices,err,error,*999)
            CALL SolverMatrices_NumberOfSolverMatricesGet(solverMatrices,numberOfSolverMatrices,err,error,*999)
            DO solverMatrixIdx=1,numberOfSolverMatrices
              NULLIFY(solverMatrix)
              CALL SolverMatrices_SolverMatrixGet(solverMatrices,solverMatrixIdx,solverMatrix,err,error,*999)
              solverMatrix%updateMatrix=.TRUE.
            ENDDO !solverMatrixIdx
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
              & FIELD_VALUES_SET_TYPE,1.0_DP,err,error,*999)
            CALL Field_ParameterSetsCopy(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE, &
              & FIELD_RESIDUAL_SET_TYPE,1.0_DP,err,error,*999)
          ELSE
            ! --- A d v e c t i o n   S o l v e r ---
            CALL Advection_PreSolve(solver,err,error,*999)
          ENDIF
          ! Update boundary conditions
          CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)

          ! --- C e l l M L    S o l v e r ---
        CASE(SOLVER_DAE_TYPE)
          ! DAE solver-set time
          CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
          CALL Solver_DAETimesSet(solver,currentTime,currentTime + timeIncrement,err,error,*999)
          CALL Solver_DAETimeStepSet(solver,timeIncrement/1000.0_DP,err,error,*999)

          ! --- S T R E E    S o l v e r ---
        CASE(SOLVER_LINEAR_TYPE)
          CALL Stree_PreSolve(solver,err,error,*999)

        CASE DEFAULT
          localError="The solve type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
            & " is invalid for a 1D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
        CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
      CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        !Pre solve for the linear solver
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        IF(solveType==SOLVER_LINEAR_TYPE) THEN
          !Update boundary conditions for mesh-movement
          CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
          NULLIFY(solvers)
          CALL Solver_SolversGet(solver,solvers,err,error,*999)
          NULLIFY(solver2)
          CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
          NULLIFY(dynamicSolver)
          CALL Solver_DynamicSolverGet(solver2,dynamicSolver,err,error,*999)
          dynamicSolver%ale=.FALSE.
          !Update material properties for Laplace mesh movement
          CALL NavierStokes_PreSolveALEUpdateParameters(solver,err,error,*999)
          !Pre solve for the linear solver
        ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
          NULLIFY(dynamicSolver)
          CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
          IF(dynamicSolver%ale) THEN
            !First update mesh and calculates boundary velocity values
            CALL NavierStokes_PreSolveALEUpdateMesh(solver,err,error,*999)
            !Then apply both normal and moving mesh boundary conditions
            CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
          ELSE
            CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
          ENDIF
        ELSE
          localError="The solve type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
            & " is invalid for an ALE Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
          CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
          !Pre solve for the linear solver
          IF(solveType==SOLVER_LINEAR_TYPE) THEN
            !TODO if first time step smooth imported mesh with respect to absolute nodal position?
            !Update boundary conditions for mesh-movement
            CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
            NULLIFY(solvers)
            CALL Solver_SolversGet(solver,solvers,err,error,*999)
            NULLIFY(solver2)
            IF(pSpecification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
              CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
            ELSE
              CALL Solvers_SolverGet(solvers,4,solver2,err,error,*999)
            ENDIF
            NULLIFY(dynamicSolver)
            CALL Solver_DynamicSolverGet(solver2,dynamicSolver,err,error,*999)
            dynamicSolver%ale=.FALSE.
            !Update material properties for Laplace mesh movement
            CALL NavierStokes_PreSolveALEUpdateParameters(solver,err,error,*999)
            !Pre solve for the dynamic solver which deals with the coupled FiniteElasticity-NavierStokes problem
          ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
            NULLIFY(dynamicSolver)
            CALL Solver_DynamicSolverGet(solver,dynamicSolver,err,error,*999)
            IF(dynamicSolver%ale) THEN
              !Apply both normal and moving mesh boundary conditions
              CALL NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
            ELSE
              CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
            ENDIF
          ELSE
            localError="The solve type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
              & " is invalid for an Elasticity ALE Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem class."
          CALL FlagError(localError,Err,Error,*999)
        END SELECT
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for NavierStokes_PreSolve of a multi physics problem class."
        CALL FlagError(localError,Err,Error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem class "//TRIM(NumberToVString(pSpecification(1),"*",err,error))// &
        & " is not valid for Navier-Stokes fluid types."
      CALL FlagError(localError,Err,Error,*999)
    END SELECT
      
    EXITS("NavierStokes_PreSolve")
    RETURN
999 ERRORSEXITS("NavierStokes_PreSolve",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PreSolve

!
!================================================================================================================================
!

  !>Sets/changes the problem subtype for a Navier-Stokes fluid type.
  SUBROUTINE NavierStokes_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("NavierStokes_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
     IF(SIZE(problemSpecification,1)<3) &
      & CALL FlagError("Navier-Stokes problem specification must have three entries.",err,error,*999)
    
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
      & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
      !All ok
    CASE(PROBLEM_OPTIMISED_NAVIER_STOKES_SUBTYPE)
      CALL FlagError("Not implemented yet.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(problem%specification(3),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_NAVIER_STOKES_EQUATION_TYPE,problemSubtype]

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
  SUBROUTINE NavierStokes_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a Navier-Stokes fluid on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(CellMLEquationsType), POINTER :: cellMLEquations
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(ControlLoopType), POINTER :: iterativeWhileLoop,iterativeWhileLoop2,iterativeWhileLoop3,simpleLoop
    TYPE(SolverEquationsType), POINTER :: solverEquations,meshSolverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(SolverType), POINTER :: solver,meshSolver,cellmlSolver
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_ProblemSetup",err,error,*999)
    
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
      !All steady state cases of Navier-Stokes
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
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a simple control loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
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
            & " is invalid for a Navier-Stokes fluid."
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
          !Set the solver to be a nonlinear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
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
          !Get the solver equations
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Transient cases and moving mesh
    CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
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
            & " is invalid for a transient Navier-Stokes fluid."
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
            & " is invalid for a transient Navier-Stokes fluid."
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
          !Set the first solver to be an CellML Evaluator for time varying boundary conditions
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Navier-Stokes boundary condition CellML evaluation solver",err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          !Set the second solver to be a first order dynamic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Navier-Stokes dynamic nonlinear solver",err,error,*999)
          CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          !setup CellML evaluator for constitutive law
          IF(pSpecification(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
            !Create the CellML evaluator solver
            NULLIFY(cellmlSolver)
            CALL Solver_NewtonCellMLEvaluatorCreate(solver,cellmlSolver,err,error,*999)
            !Link the CellML evaluator solver to the solver
            CALL Solver_LinkedSolverAdd(solver,cellmlSolver,SOLVER_CELLML_EVALUATOR_TYPE,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver equations
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
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
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          !Get the boundary condition solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the CellML equations
          NULLIFY(cellMLEquations)
          CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          !Set the time dependence
          CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
          IF(pSpecification(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
            !Get the multiscale solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            !Get the CellML evaluator solver
            NULLIFY(cellMLSolver)
            CALL Solver_NewtonLinkedCellMLSolverGet(solver,cellMLSolver,err,error,*999)
            !Create the CellML equations
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(cellMLSolver,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Get the CellML boundary condition solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(cellMLEquations)
          CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
          !Finish the CellML equations creation
          CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          IF(pSpecification(3)==PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE) THEN
            !Get the multiscale CellML solver
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            !Get the CellML evaluator solver
            NULLIFY(cellMLSolver)
            CALL Solver_NewtonLinkedCellMLSolverGet(solver,cellMLSolver,err,error,*999)
            !Get the CellML equations for the CellML evaluator solver
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(cellmlSolver,cellMLEquations,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a CellML setup for a  transient Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a transient Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
      ! Multiscale: 3D/1D/0D coupled
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
            & " is invalid for a transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          NULLIFY(iterativeWhileLoop)
          ! The 3D-1D boundary value iterative coupling loop
          CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
          NULLIFY(iterativeWhileLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
          CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,100,err,error,*999)
          CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E-6_DP,err,error,*999)
          CALL ControlLoop_RelativeToleranceSet(iterativeWhileLoop,0.001_DP,err,error,*999)
          IF (pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"3D-0D Iterative Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The 3D Navier-Stokes simple loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"3D Navier-Stokes",err,error,*999)
          ELSE
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"3D-1D Iterative Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The 1D-0D boundary value iterative coupling loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop2,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop2,100,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,0.1_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop2,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop2,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes coupling loop
            NULLIFY(iterativeWhileLoop3)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,2,iterativeWhileLoop3,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop3,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop3,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop3,1.0E10_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop3,"1D Iterative Loop",err,error,*999)
            ! The 3D Navier-Stokes simple loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"3D Navier-Stokes",err,error,*999)
          ENDIF
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
            & " is invalid for a transient Navier-Stokes fluid."
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
          ! The 3D-1D iterative coupling loop
          NULLIFY(iterativeWhileLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
          IF (pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"DAE Solver",err,error,*999)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes 3D Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          ELSE
            ! Iterative loop 2 couples 1D and 0D, checking convergence
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"DAE Solver",err,error,*999)
            ! Iterative loop 3 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop3)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,2,iterativeWhileLoop3,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop3,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- 1 D   N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes 1D Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes 3D Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          NULLIFY(iterativeWhileLoop)
          IF(pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            !Finish the 0D solvers
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            !Finish the 3D solver
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
            !Finish the 0D solvers
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            !Finish the 1D solvers
            NULLIFY(iterativeWhileLoop3)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,2,iterativeWhileLoop3,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop3,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            !Finish the 3D solver
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        !Get the control loop
        IF(pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          ! 3D0D subloop
          NULLIFY(iterativeWhileLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
          ! 3D subloop
          NULLIFY(simpleLoop)
          CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !Create the solver equations
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
              & err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          ! 3D1D subloop
          NULLIFY(iterativeWhileLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
          ! 1D0D subloop
          NULLIFY(iterativeWhileLoop2)
          CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
          ! 1D subloop
          NULLIFY(iterativeWhileLoop3)
          CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,2,iterativeWhileLoop3,err,error,*999)
          ! 3D subloop
          NULLIFY(simpleLoop)
          CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,simpleLoop,err,error,*999)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            ! 1D NS/C subloop
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop3,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- 1 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !Create the solver equations
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,&
              & err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            ! 1D NS/C subloop
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop3,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- 1 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- 3 D  N A V I E R   S T O K E S --
            !
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a Navier-Stokes fluid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        !Create the CELLML solver equations
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          IF (pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            ! 3D-1D loop
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! 0D loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          ELSE
            ! 3D-1D loop
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! 1D-0D loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
            ! 0D loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          IF (pSpecification(3)==PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(controlLoopRoot)
            CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
            NULLIFY(controlLoop)
            CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
            ! 3D-0D
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! 0D
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ELSE
            NULLIFY(controlLoopRoot)
            CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
            NULLIFY(controlLoop)
            CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
            ! 3D-1D
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! 1D-0D
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,iterativeWhileLoop2,err,error,*999)
            ! 0D
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a 3D-1d-0D Navier-Stokes fluid type."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &    !1D Navier-Stokes
      & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &    !  with coupled 0D boundaries
      & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, & !  with coupled advection
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, & !  with coupled 0D boundaries and advection
      & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &      !  with stree
      & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)     !  with stree and advection

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
            & " is invalid for Coupled1dDaeNavierStokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Time Loop
          NULLIFY(controlLoop)
          CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          IF(pSpecification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
            ! The 1D-0D boundary value iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop2,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop2,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
              & err,error,*999)
          ELSE IF(pSpecification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            ! The 1D-0D boundary value iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,2,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop2,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop2,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
              & err,error,*999)
            ! The simple Advection solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"Advection",err,error,*999)
          ELSE IF(pSpecification(3) == PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E3_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
          ELSE IF(pSpecification(3) == PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,2,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,1.0E6_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
            ! The simple Advection solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"Advection",err,error,*999)
          ELSE IF(pSpecification(3) == PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
            ! The 1D-0D boundary value iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop2,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop2,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop",err,error,*999)
          ELSE IF(pSpecification(3) == PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            ! The 1D-0D boundary value iterative coupling loop
            CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,2,err,error,*999)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop,0.1_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop,"1D-0D Iterative Coupling Convergence Loop",err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsSet(iterativeWhileLoop,2,err,error,*999)
            ! The simple CellML solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"0D CellML solver Loop",err,error,*999)
            ! The Characteristics branch solver/ Navier-Stokes iterative coupling loop
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            CALL ControlLoop_TypeSet(iterativeWhileLoop2,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
            CALL ControlLoop_MaximumIterationsSet(iterativeWhileLoop2,1000,err,error,*999)
            CALL ControlLoop_AbsoluteToleranceSet(iterativeWhileLoop2,1.0E6_DP,err,error,*999)
            CALL ControlLoop_LabelSet(iterativeWhileLoop2,"1D Characteristic/NSE branch value convergence Loop", &
              & err,error,*999)
            ! The simple Advection solver loop
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            CALL ControlLoop_TypeSet(simpleLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
            CALL ControlLoop_LabelSet(simpleLoop,"Advection",err,error,*999)
          ENDIF
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a 1d transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Create the solvers
      CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        NULLIFY(controlLoop)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            ! Simple loop 1 contains the Advection solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Advection Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_LINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"DAE Solver",err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Advection Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_LINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Linear Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Advection Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_LINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DAE_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"DAE Solver",err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Simple loop 1 contains the 0D/CellML DAE solver
            ! (this subloop holds 1 solver)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(simpleLoop,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Linear Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL Solvers_CreateStart(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Characteristic Solver",err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL Solver_LabelSet(solver,"Navier-Stokes Solver",err,error,*999)
            CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          IF(pSpecification(3)==PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            CALL Solvers_CreateFinish(solvers,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a 1d transient Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Create the solver equations
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solvers) !WHY IS THE SOLVERS GO HERE WHEN IT WAS JUST GOT ABOVE????
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE(PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            ! Iterative loop couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,2,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- A D V E C T I O N --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE(PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
            ! Iterative loop 1 couples 1D and 0D problems, checking convergence of boundary values
            ! (this subloop holds 2 subloops)
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
            !
            !-- D A E --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            ! Iterative loop 2 couples N-S and characteristics, checking convergence of branch values
            ! (this subloop holds 2 solvers)
            NULLIFY(iterativeWhileLoop2)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,2,iterativeWhileLoop2,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(iterativeWhileLoop2,solvers,err,error,*999)
            !
            !-- C H A R A C T E R I S T I C --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
            !
            !-- N A V I E R   S T O K E S --
            !
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for a Navier-Stokes equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Create the CELLML solver equations
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          IF(pSpecification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
            & pSpecification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
          ELSE
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          ENDIF
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL CellMLEquations_CreateStart(solver,cellMLEquations,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(cellMLEquations,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(cellMLEquations,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for cellML equations setup Navier-Stokes equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          IF(pSpecification(3) == PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
            & pSpecification(3) == PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
            NULLIFY(iterativeWhileLoop)
            CALL ControlLoop_SubLoopGet(controlLoop,1,iterativeWhileLoop,err,error,*999)
            NULLIFY(simpleLoop)
            CALL ControlLoop_SubLoopGet(iterativeWhileLoop,1,simpleLoop,err,error,*999)
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(simpleLoop,solvers,err,error,*999)
          ELSE
            NULLIFY(solvers)
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          ENDIF
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
            & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
            NULLIFY(solver)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            NULLIFY(cellMLEquations)
            CALL Solver_CellMLEquationsGet(solver,cellMLEquations,err,error,*999)
            CALL CellMLEquations_CreateFinish(cellMLEquations,err,error,*999)
          CASE DEFAULT
            localError="The third problem specification of "// &
              & TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is not valid for cellML equations setup Navier-Stokes fluid mechanics problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a CellML setup for a 1D Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a 1d transient Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
      !Quasi-static Navier-Stokes
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
            & " is invalid for a quasistatic Navier-Stokes fluid."
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
            & " is invalid for a quasistatic Navier-Stokes fluid."
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
          !Set the solver to be a nonlinear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a quasistatic Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver equations
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a quasistatic Navier-Stokes equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a quasistatic Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Navier-Stokes ALE cases
    CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(problemSetup%setupType)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing????
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a ALE Navier-Stokes fluid."
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
            & " is invalid for a ALE Navier-Stokes fluid."
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
          !Set the first solver to be a linear solver for the Laplace mesh movement problem
          NULLIFY(meshSolver)
          CALL Solvers_SolverGet(solvers,1,meshSolver,err,error,*999)
          CALL Solver_TypeSet(meshSolver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(meshSolver,SOLVER_PETSC_LIBRARY,err,error,*999)
          !Set the solver to be a first order dynamic solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          !Finish the solvers creation
          CALL Solvers_CreateFinish(solvers,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a ALE Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(problemSetup%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(meshSolver)
          CALL Solvers_SolverGet(solvers,1,meshSolver,err,error,*999)
          !Create the solver equations
          NULLIFY(meshSolverEquations)
          CALL SolverEquations_CreateStart(meshSolver,meshSolverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(meshSolverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(meshSolverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(meshSolverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          NULLIFY(controlLoopRoot)
          CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          !Get the solver equations
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(meshSolver)
          CALL Solvers_SolverGet(solvers,1,meshSolver,err,error,*999)
          NULLIFY(meshSolverEquations)
          CALL Solver_SolverEquationsGet(meshSolver,meshSolverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(meshSolverEquations,err,error,*999)
          !Get the solver equations
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a Navier-Stokes fluid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a ALE Navier-Stokes fluid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_ProblemSetup")
    RETURN
999 ERRORSEXITS("NavierStokes_ProblemSetup",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnMeshComponent, &
      & columnXiIdx,componentIdx,coordinateIdx,derivativeIdx,elementVersionNumber,esSpecification(3),firstNode,gaussPointIdx, &
      & lastNode,maxColumnElementDOFIdx,minColumnElementDOFIdx,maxRowElementDOFIdx,minRowElementDOFIdx,nodeNumber,nodeIdx, &
      & numberOfColumnElementParameters,numberOfDimensions,numberOfElements,numberOfElementNodes,numberOfGauss,numberOfParameters, &
      & numberOfResidualComponents,numberOfRowsComponents,numberOfRowElementParameters,numberOfVersions,numberOfXi, &
      & residualVariableType,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowMeshComponent,rowXiIdx,rowsVariableType, &
      & versionIdx,xiIdx
    REAL(DP) :: aDeriv,a0Deriv,a0Param,alpha,area,AUpwind,aValue,beta,columnPhi,dColumnPhidXi(3),dRowPhidXi(3),dXidX(3,3), &
      & eDeriv,eParam,g0Param,gaussWeight,hDeriv,hParam,jacobian,jacobianGaussWeight,kappa,Lref,mass,momentum,Mref,muParam, &
      & muScale,normal,normalWave,pExternal,pressure,qDeriv,QUpwind,qValue,rhoParam,rowPhi,sum,Tref,uDeriv(3,3),uValue(3), &
      & wValue(3),x(3)
    LOGICAL :: firstAssemblyRHS,firstAssemblyStiffness,update,updateDamping,updateLinear,updateMatrices,updateResidual, &
      & updateRHS,updateStiffness
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,independentBasis,rowBasis
    TYPE(CoordinateSystemType), POINTER :: coordinateSystem
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition,independentDecomposition
    TYPE(DecompositionElementsType), POINTER :: dependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: dependentDecompositionTopology
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,independentDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,&
      & independentDomainElements,rowDomainElements
    TYPE(DomainNodesType), POINTER :: dependentDomainNodes
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology, &
      & independentDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix,jacobianMatrix
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField,independentField
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,dependentInterpPoint,independentInterpPoint, &
      & uMaterialsInterpPoint,vMaterialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldInterpolationParametersType), POINTER :: geometricInterpParameters,dependentInterpParameters, &
      & independentInterpParameters,uMaterialsInterpParameters,vMaterialsInterpParameters
    TYPE(FieldVariableType), POINTER :: geometricVariable,residualVariable,rowsVariable,uDependentVariable, &
      & uIndependentVariable,uMaterialsVariable,vMaterialsVariable
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,quadratureScheme,rowQuadratureScheme
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_FiniteElementResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
 
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateResidual,err,error,*999)
    NULLIFY(rhsMapping)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    firstAssemblyRHS=.FALSE.
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
      CALL EquationsMatricesRHS_FirstAssemblyGet(rhsVector,firstAssemblyRHS,err,error,*999)
    ENDIF
    
    NULLIFY(linearMatrices)
    NULLIFY(dynamicMatrices)
    NULLIFY(stiffnessMatrix)
    updateStiffness=.FALSE.
    NULLIFY(dampingMatrix)
    updateDamping=.FALSE.
   
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      NULLIFY(stiffnessMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      NULLIFY(dampingMatrix)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateResidual.OR.updateRHS.OR.updateMatrices)

    IF(update) THEN

      uValue=0.0_DP
      uDeriv=0.0_DP
      dXidX=0.0_DP
    
      !Set general and specific pointers
      NULLIFY(region)
      CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
      NULLIFY(coordinateSystem)
      CALL Region_CoordinateSystemGet(region,coordinateSystem,err,error,*999)
      CALL CoordinateSystem_DimensionGet(coordinateSystem,numberOfDimensions,err,error,*999)
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(residualVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,residualVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(residualVariable,residualVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(residualVariable,numberOfResidualComponents,err,error,*999)
      
      NULLIFY(dependentField)
      CALL EquationsInterpolation_DependentFieldGet(equationsInterpolation,dependentField,err,error,*999)
      NULLIFY(uDependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,uDependentVariable,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsInterpolation_GeometricFieldGet(equationsInterpolation,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      NULLIFY(materialsField)      
      CALL EquationsInterpolation_MaterialsFieldGet(equationsInterpolation,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      NULLIFY(geometricDomain)
      CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
      NULLIFY(geometricDomainTopology)
      CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
      NULLIFY(geometricDomainElements)
      CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentDomainNodes)
      CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      
      NULLIFY(jacobianMatrix)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
        NULLIFY(linearMatrices)
        CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
      CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
        NULLIFY(linearMatrices)
        CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
      CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
        NULLIFY(dampingMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
        NULLIFY(dampingMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsInterpolation_IndependentFieldGet(equationsInterpolation,independentField,err,error,*999)
        NULLIFY(uIndependentVariable)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,uIndependentVariable,err,error,*999)
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,3,alpha,err,error,*999)
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,8,g0Param,err,error,*999)
        CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
        NULLIFY(vMaterialsInterpParameters)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & vMaterialsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vMaterialsInterpParameters,err,error,*999)
      CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
        NULLIFY(dampingMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
        NULLIFY(independentField)
        CALL EquationsInterpolation_IndependentFieldGet(equationsInterpolation,independentField,err,error,*999)
        NULLIFY(independentDecomposition)
        CALL Field_DecompositionGet(independentField,independentDecomposition,err,error,*999)
        NULLIFY(independentDomain)
        CALL Decomposition_DomainGet(independentDecomposition,0,independentDomain,err,error,*999)
        NULLIFY(independentDomainTopology)
        CALL Domain_DomainTopologyGet(independentDomain,independentDomainTopology,err,error,*999)
        NULLIFY(independentDomainElements)
        CALL DomainTopology_DomainElementsGet(independentDomainTopology,independentDomainElements,err,error,*999)
        NULLIFY(independentBasis)
        CALL DomainElements_ElementBasisGet(independentDomainElements,elementNumber,independentBasis,err,error,*999)
        NULLIFY(dynamicMatrices)
        CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
        NULLIFY(stiffnessMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
        CALL EquationsMatrix_FirstAssemblyGet(stiffnessMatrix,firstAssemblyStiffness,err,error,*999)
        NULLIFY(dampingMatrix)
        CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
        CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
        NULLIFY(independentInterpParameters)
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        NULLIFY(independentInterpPoint)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,independentInterpPoint, &
          & err,error,*999)        
        CALL Field_InterpolationParametersElementGet(FIELD_MESH_VELOCITY_SET_TYPE,elementNumber,independentInterpParameters, &
          & err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,residualVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
        & err,error,*999)
      
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      
      NULLIFY(uMaterialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uMaterialsInterpParameters, &
        & err,error,*999)
      NULLIFY(uMaterialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uMaterialsInterpPoint, &
        & err,error,*999)
      
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,uMaterialsInterpParameters, &
        & err,error,*999)

      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
      IF(ASSOCIATED(equationsAnalytic)) THEN
        CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
      ENDIF
      
      !Loop over Gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,uMaterialsInterpPoint, &
          & err,error,*999)

        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        DO xiIdx=1,numberOfXi
          DO coordinateIdx=1,numberOfDimensions
            dXidX(xiIdx,coordinateIdx)=geometricInterpPointMetrics%dXidX(xiIdx,coordinateIdx)
          ENDDO !coordinateIdx
        ENDDO !xiIdx
        
        wValue=0.0_DP
        IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
            & err,error,*999)
          wValue(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        ! Get the constitutive law (non-Newtonian) viscosity based on shear rate
        IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
          ! Note the constant from the U_VARIABLE is a scale factor
          muScale = uMaterialsInterpPoint%values(1,NO_PART_DERIV)
          ! Get the gauss point based value returned from the CellML solver
          CALL Field_ParameterSetGetLocalGaussPoint(materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & gaussPointIdx,elementNumber,1,muParam,err,error,*999)
          muParam=muParam*muScale
        ELSE
          muParam = uMaterialsInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        rhoParam = uMaterialsInterpPoint%values(2,NO_PART_DERIV)

        !Start with matrix calculations
        IF(esSpecification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          !Loop over field components
          rowElementDOFIdx=0
          DO rowComponentIdx=1,numberOfDimensions
            CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
            NULLIFY(rowDomain)
            CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
            NULLIFY(rowDomainTopology)
            CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
            NULLIFY(rowDomainElements)
            CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
            NULLIFY(rowBasis)
            CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
            CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
            NULLIFY(rowQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_GaussWeightGet(rowQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            jacobianGaussWeight=jacobian*gaussWeight
            DO rowElementParameterIdx=1,numberOfRowElementParameters
              rowElementDOFIdx=rowElementDOFIdx+1
              !Calculate some general values
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussPointIdx,rowPhi,err,error,*999)
              DO xiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                  & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,dRowPhidXi(xiIdx),err,error,*999)
              ENDDO !xiIdx
              columnElementDOFIdx=0
              IF(updateMatrices) THEN
                !Loop over element columns
                DO columnComponentIdx=1,numberOfResidualComponents
                  CALL FieldVariable_ComponentMeshComponentGet(residualVariable,columnComponentIdx,columnMeshComponent, &
                    & err,error,*999)
                  NULLIFY(columnDomain)
                  CALL Decomposition_DomainGet(dependentDecomposition,columnMeshComponent,columnDomain,err,error,*999)
                  NULLIFY(columnDomainTopology)
                  CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                  NULLIFY(columnDomainElements)
                  CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                  NULLIFY(columnBasis)
                  CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                  CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                  NULLIFY(columnQuadratureScheme)
                  CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                    & err,error,*999)
                  DO columnElementParameterIdx=1,numberOfColumnElementParameters
                    columnElementDOFIdx=columnElementDOFIdx+1
                    !Calculate some general values
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                    DO xiIdx=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,dColumnPhidXi(xiIdx),err,error,*999)
                    ENDDO !xiIdx
                    !Laplace only matrix
                    IF(updateStiffness) THEN
                      !LAPLACE TYPE
                      IF(columnComponentIdx==rowComponentIdx) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO coordinateIdx=1,numberOfDimensions
                          DO rowXiIdx=1,numberOfXi
                            DO columnXiIdx=1,numberOfXi
                              sum=sum+muParam*dColumnPhidXi(columnXiIdx)*dXidX(columnXiIdx,coordinateIdx)* &
                                & dRowPhidXi(rowXiIdx)*dXidX(rowXiIdx,coordinateIdx)
                            ENDDO !columnXiIdx
                          ENDDO !rowXiIdx
                        ENDDO !coordinateIdx
                        !Calculate MATRIX
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                    ENDIF
                    !General matrix
                    IF(updateStiffness) THEN
                      !GRADIENT TRANSPOSE TYPE
                      IF(esSpecification(3)/=EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE) THEN
                        IF(columnComponentIdx<numberOfResidualComponents) THEN
                          sum=0.0_DP
                          !Calculate sum
                          DO rowXiIdx=1,numberOfXi
                            DO columnXiIdx=1,numberOfXi
                              !note rowComponentIdx/columnComponentIdx derivative in dXidX
                              sum=sum+muParam*dColumnPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)* &
                                & dRowPhidXi(columnXiIdx)*dXidX(columnXiIdx,columnComponentIdx)
                            ENDDO !columnXiIdx
                          ENDDO !rowXiIdx
                          !Calculate MATRIX
                          stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                            & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) &
                            & +sum*jacobianGaussWeight
                        ENDIF
                      ENDIF
                    ENDIF
                    !Contribution through ALE
                    IF(updateStiffness) THEN
                      !GRADIENT TRANSPOSE TYPE
                      IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                        & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
                        IF(columnComponentIdx==rowComponentIdx) THEN
                          sum=0.0_DP
                          !Calculate sum
                          DO coordinateIdx=1,numberOfDimensions
                            DO xiIdx=1,numberOfXi
                              sum=sum-rhoParam*wValue(coordinateIdx)*dColumnPhidXi(xiIdx)*dXidX(xiIdx,coordinateIdx)*rowPhi
                            ENDDO !xiIdx
                          ENDDO !coordinateIdx
                          !Calculate MATRIX
                          stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                            & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                            & sum*jacobianGaussWeight
                        ENDIF
                      ENDIF
                    ENDIF
                    !Pressure contribution (B transpose)
                    IF(updateStiffness) THEN
                      !LAPLACE TYPE
                      IF(columnComponentIdx==numberOfResidualComponents) THEN
                        sum=0.0_DP
                        !Calculate sum
                        DO rowXiIdx=1,numberOfXi
                          sum=sum-columnPhi*dRowPhidXi(rowXiIdx)*dXidX(rowXiIdx,rowComponentIdx)
                        ENDDO !rowXiIdx
                        !Calculate MATRIX
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                    ENDIF
                    !Damping matrix
                    IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
                      & esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
                      & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                      & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
                      & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE .OR. &
                      & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                      IF(updateDamping) THEN
                        IF(columnComponentIdx==rowComponentIdx) THEN
                          sum=0.0_DP
                          !Calculate sum
                          sum=rowPhi*columnPhi*rhoParam
                          !Calculate MATRIX
                          dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                            & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO !columnElementParameterIdx
                ENDDO !columnComponentIdx
              ENDIF
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
          !Analytic RHS vector
          IF(firstAssemblyRHS) THEN
            IF(updateRHS) THEN
              IF(ASSOCIATED(equationsAnalytic)) THEN
                IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                  x(1:numberOfDimensions)=geometricInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
                  rowElementDOFIdx=0
                  DO rowComponentIdx=1,numberOfDimensions
                    CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
                    NULLIFY(rowDomain)
                    CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
                    NULLIFY(rowDomainTopology)
                    CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
                    NULLIFY(rowDomainElements)
                    CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
                    NULLIFY(rowBasis)
                    CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
                    CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
                    NULLIFY(rowQuadratureScheme)
                    CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
                    CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
                    jacobianGaussWeight=jacobian*gaussWeight
                    DO rowElementParameterIdx=1,NumberOfRowElementParameters
                      rowElementDOFIdx=rowElementDOFIdx+1
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                        & gaussPointIdx,rowPhi,err,error,*999)
                      !note rowComponentIdx value derivative
                      sum=0.0_DP
                      IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(-2.0_DP/3.0_DP*(x(1)**3*rhoParam+3.0_DP*muParam*10.0_DP**2- &
                            & 3.0_DP*rhoParam*x(2)**2*x(1))/(10.0_DP**4))
                        ENDIF
                      ELSE IF(equationsAnalytic%analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(-4.0_DP*muParam/10.0_DP/10.0_DP*EXP((X(1)-X(2))/10.0_DP))
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(16.0_DP*muParam*PI**2/10.0_DP**2*COS(2.0_DP*PI*x(2)/10.0_DP)* &
                            & COS(2.0_DP*PI*x(1)/10.0_DP)- &
                            & 2.0_DP*COS(2.0_DP*PI*x(2)/10.0_DP)*SIN(2.0_DP*PI*x(2)/10.0_DP)*rhoParam*PI/10.0_DP)
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=rowPhi*(2.0_DP*SIN(x(1))*COS(x(2)))*muParam
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(-2.0_DP*COS(x(1))*SIN(x(2)))*muParam
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5) THEN
                        !do nothing
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(-2.0_DP/3.0_DP*(rhoParam*x(1)**3+6.0_DP*rhoParam*x(1)*x(3)*x(2)+ &
                            & 6.0_DP*muParam*10.0_DP**2- &
                            & 3.0_DP*rhoParam*X(2)**2*x(1)-3.0_DP*rhoParam*x(3)*x(1)**2-3.0_DP*rhoParam*x(3)*x(2)**2)/ &
                            & (10.0_DP**4))
                        ELSE IF(rowComponentIdx==3) THEN
                          !Calculate sum
                          sum=rowPhi*(-2.0_DP/3.0_DP*(6.0_DP*rhoParam*x(1)*x(3)*x(2)+rhoParam*X(1)**3+ &
                            & 6.0_DP*muParam*10.0_DP**2- &
                            & 3.0_DP*rhoParam*x(1)*x(3)**2-3.0_DP*rhoParam*x(2)*x(1)**2-3.0_DP*rhoParam*x(2)*x(3)**2)/ &
                            & (10.0_DP**4))
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*((-4.0_DP*muParam*EXP((x(1)-x(2))/10.0_DP)-2.0_DP*muParam*EXP((x(2)-x(3))/10.0_DP)+ &
                            & rhoParam*EXP((x(3)-x(2))/10.0_DP)*10.0_DP)/10.0_DP**2)
                        ELSE IF(rowComponentIdx==3) THEN
                          !Calculate sum
                          sum=rowPhi*(-(4.0_DP*muParam*EXP((x(3)-x(1))/10.0_DP)+2.0_DP*muParam*EXP((x(2)-x(3))/10.0_DP)+ &
                            & rhoParam*EXP((x(3)-x(2))/10.0_DP)*10.0_DP)/10.0_DP** 2)
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          sum=0.0_DP
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          sum=rowPhi*(2.0_DP*COS(2.0_DP*PI*x(2)/10.0_DP)*(18.0_DP*COS(2.0_DP*PI*x(1)/10.0_DP)* &
                            & muParam*PI*SIN(2.0_DP*PI*x(3)/10.0_DP)-3.0_DP*rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)**2* &
                            & SIN(2.0_DP*PI*x(2)/10.0_DP)*10.0_DP-2.0_DP*rhoParam*SIN(2.0_DP*PI*x(2)/10.0_DP)*10.0_DP+ &
                            & 2.0_DP*rhoParam*SIN(2.0_DP*PI*x(2)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*x(3)/10.0_DP)**2)*PI/ &
                            & 10.0_DP**2)
                        ELSE IF(rowComponentIdx==3) THEN
                          !Calculate sum
                          sum=rowPhi*(-2.0_DP*PI*COS(2.0_DP*PI*x(3)/10.0_DP)*rhoParam*SIN(2.0_DP*PI*x(3)/10.0_DP)* &
                            & (-1.0_DP+COS(2.0_DP*PI*x(2)/10.0_DP)**2)/10.0_DP)
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4) THEN
                        IF(rowComponentIdx==1) THEN
                          !Calculate sum
                          !sum=rowPhi*(2.0_DP*SIN(x(1))*COS(x(2)))*muParam
                        ELSE IF(rowComponentIdx==2) THEN
                          !Calculate sum
                          !sum=rowPhi*(-2.0_DP*COS(x(1))*SIN(x(2)))*muParam
                        ENDIF
                      ELSE IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5) THEN
                        !do nothing
                      ENDIF
                      !Calculate RHS vector
                      rhsVector%elementVector%vector(rowElementDOFIdx)= &
                        & rhsVector%elementVector%vector(rowElementDOFIdx)+sum*jacobianGaussWeight
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ELSE
                  rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
                ENDIF
              ENDIF
            ENDIF
          ENDIF

          !Calculate nonlinear vector
          IF(updateResidual) THEN
            ! Get interpolated velocity and velocity gradient values for nonlinear term
            uValue(1:numberOfDimensions)=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            DO xiIdx=1,numberOfXi
              uDeriv(1:numberOfDimensions,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))= &
                & dependentInterpPoint%values(1:numberOfDimensions,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
            ENDDO !xiIdx
            !Here wValues must be ZERO if ALE part of linear matrix
            wValue=0.0_DP
            rowElementDOFIdx=0
            DO rowComponentIdx=1,numberOfDimensions
              CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
              NULLIFY(rowDomain)
              CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
              NULLIFY(rowDomainTopology)
              CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
              NULLIFY(rowDomainElements)
              CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
              NULLIFY(rowBasis)
              CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
              CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
              NULLIFY(rowQuadratureScheme)
              CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
              CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
              jacobianGaussWeight=jacobian*gaussWeight
              DO rowElementParameterIdx=1,numberOfRowElementParameters
                rowElementDOFIdx=rowElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                  & gaussPointIdx,rowPhi,err,error,*999)
                !note rowComponentIdx value derivative
                sum=0.0_DP
                ! Convective form
                DO xiIdx=1,numberOfXi
                  sum=sum+rhoParam*rowPhi*( &
                    & (uValue(1))*(uDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,1))+ &
                    & (uValue(2))*(uDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,2))+ &
                    & (uValue(3))*(uDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,3)))
                ENDDO !xiIdx

                residualVector%elementResidual%vector(rowElementDOFIdx)= &
                  & residualVector%elementResidual%vector(rowElementDOFIdx)+sum*jacobianGaussWeight

              ENDDO !rowElementParameterIdx
            ENDDO !rowComponentIdx
          ENDIF
        ENDIF

        !------------------------------------------------------------------
        ! R e s i d u a l - b a s e d    S t a b i l i s a t i o n
        !------------------------------------------------------------------
        IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
          CALL NavierStokes_ResidualBasedStabilisation(equationsSet,elementNumber,gaussPointIdx, &
            & muParam,rhoParam,.FALSE.,err,error,*999)
        ENDIF
        
        !------------------------------------------------------------------
        ! 1 D  T r a n s i e n t
        !------------------------------------------------------------------
        !Start with matrix calculations
        IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
          qValue=dependentInterpPoint%values(1,NO_PART_DERIV)
          qDeriv=dependentInterpPoint%values(1,FIRST_PART_DERIV)
          aValue=dependentInterpPoint%values(2,NO_PART_DERIV)
          aDeriv=dependentInterpPoint%values(2,FIRST_PART_DERIV)
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,vMaterialsInterpPoint, &
            & err,error,*999)
          a0Param=vMaterialsInterpPoint%values(1,NO_PART_DERIV)
          a0Deriv=vMaterialsInterpPoint%values(1,FIRST_PART_DERIV)
          eParam=vMaterialsInterpPoint%values(2,NO_PART_DERIV)
          eDeriv=vMaterialsInterpPoint%values(2,FIRST_PART_DERIV)
          hParam=vMaterialsInterpPoint%values(3,NO_PART_DERIV)
          hDeriv=vMaterialsInterpPoint%values(3,FIRST_PART_DERIV)
          beta = (4.0_DP*(SQRT(PI))*eParam*hParam)/(3.0_DP*a0Param)  !(kg/m2/s2)
          kappa = 8.0_DP*PI*muParam/rhoParam ! viscous resistance operator

          ! If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
          IF(aValue < a0Param*0.001_DP) aValue = a0Param*0.001_DP
 
          rowElementDOFIdx=0
          !Loop Over Element Rows
          DO rowComponentIdx=1,numberOfResidualComponents
            CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
            NULLIFY(rowDomain)
            CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
            NULLIFY(rowDomainTopology)
            CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
            NULLIFY(rowDomainElements)
            CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
            NULLIFY(rowBasis)
            CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
            CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
            NULLIFY(rowQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
            dXidX=0.0_DP
            !Calculate dxi_dx in 3D
            DO xiIdx=1,numberOfXi
              DO coordinateIdx=1,numberOfDimensions
                dXidX(1,1)=dXidX(1,1)+(geometricInterpPointMetrics%dXidX(xiIdx,coordinateIdx))**2.0_DP
              ENDDO !coordinateIdx
            ENDDO !xiIdx
            dXidX(1,1)=SQRT(dXidX(1,1))
            !Loop Over Element rows
            DO rowElementParameterIdx=1,numberOfRowElementParameters
              rowElementDOFIdx=rowElementDOFIdx+1
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussPointIdx,rowPhi,err,error,*999)
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,FIRST_PART_DERIV, &
                & gaussPointIdx,dRowPhidXi(1),err,error,*999)
              columnElementDOFIdx=0
              IF(updateMatrices) THEN
                !Loop Over Element Columns
                DO columnComponentIdx=1,numberOfResidualComponents
                  CALL FieldVariable_ComponentMeshComponentGet(residualVariable,columnComponentIdx,columnMeshComponent, &
                    & err,error,*999)
                  NULLIFY(columnDomain)
                  CALL Decomposition_DomainGet(dependentDecomposition,columnMeshComponent,columnDomain,err,error,*999)
                  NULLIFY(columnDomainTopology)
                  CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                  NULLIFY(columnDomainElements)
                  CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                  NULLIFY(columnBasis)
                  CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                  NULLIFY(columnQuadratureScheme)
                  CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                    & err,error,*999)
                  DO columnElementParameterIdx=1,columnBasis%numberOfElementParameters
                    columnElementDOFIdx=columnElementDOFIdx+1
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & FIRST_PART_DERIV,gaussPointIdx,dColumnPhidXi(1),err,error,*999)
                    !
                    !-- D A M P I N G  M A T R I X --
                    !
                    IF(updateDamping) THEN
                      !Momentum Equation, dQ/dt
                      IF(rowComponentIdx==1.AND.columnComponentIdx==1) THEN
                        sum=rowPhi*columnPhi
                        dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                      !Mass Equation, dA/dt
                      IF(rowComponentIdx==2.AND.columnComponentIdx==2) THEN
                        sum=rowPhi*columnPhi
                        dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                    ENDIF
                    !
                    !-- S T I F F N E S S  M A T R I X --
                    !
                    IF(updateStiffness) THEN
                      IF(rowComponentIdx==1.AND.columnComponentIdx==2) THEN
                        !Momentum Equation, linearisable A0 terms
                        sum=-rowPhi*columnPhi*(beta*SQRT(a0Param)/rhoParam)*(hDeriv/hParam + eDeriv/eParam)*dXidX(1,1)
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                        !Momentum Equation, gravitational force
                        sum=rowPhi*columnPhi*g0Param
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                      !Mass Equation, dQ/dX, flow derivative
                      IF(rowComponentIdx==2 .AND. columnComponentIdx==1) THEN
                        sum=rowPhi*dColumnPhidXi(1)*dXidX(1,1)
                        stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                      ENDIF
                    ENDIF                    
                  ENDDO !columnElementParameterIdx
                ENDDO !columnComponentIdx
              ENDIF
              !
              !-- N O N L I N E A R  V E C T O R --
              !
              IF(updateResidual) THEN
                !Momentum Equation
                IF(rowComponentIdx==1) THEN
                  sum=((2.0_DP*alpha*(qValue/aValue)*qDeriv - &
                    & (alpha*((qValue/aValue)**2.0_DP)*aDeriv)+(beta/rhoParam)* &          !Convective
                    & ((SQRT(aValue)/2.0_DP)*aDeriv+ &                                     !A  gradient
                    & (aValue/(2.0_DP*SQRT(a0Param))-(aValue**1.5_DP)/a0Param)*a0Deriv+ &  !A0 gradient
                    & (aValue*(SQRT(aValue)))*(hDeriv/hParam) + &                          !H  gradient (nonlinear part)
                    & (aValue*(SQRT(aValue)))*(eDeriv/eParam)))* &                         !E  gradient (nonlinear part)
                    & dXidX(1,1)+kappa*(qValue/aValue))*rowPhi                             !Viscosity
                  residualVector%elementResidual%vector(rowElementDOFIdx)= &
                    & residualVector%elementResidual%vector(rowElementDOFIdx)+sum*jacobianGaussWeight
                ENDIF
              ENDIF
              
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
        ENDIF
      ENDDO !gaussPointIdx

      IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
        IF(updateResidual) THEN
          CALL Basis_NumberOfLocalNodesGet(dependentBasis,numberOfElementNodes,err,error,*999)
          CALL DomainElements_MaxElementParametersGet(dependentDomainElements,numberOfParameters,err,error,*999)
          CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,1,firstNode,err,error,*999)
          CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,numberOfParameters,lastNode,err,error,*999)
          !Get material constants
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,rhoParam,err,error,*999)
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,4,pExternal,err,error,*999)
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,5,Lref,err,error,*999)
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,6,Tref,err,error,*999)
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,7,Mref,err,error,*999)
          !
          !-- P R E S S U R E    C A L C U L A T I O N --
          !
          !Loop over the element nodes and versions
          NULLIFY(dependentDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
          NULLIFY(dependentDecompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
            & err,error,*999)
          CALL DecompositionElements_NumberOfElementsGet(dependentDecompositionElements,numberOfElements,err,error,*999)
          DO nodeIdx=1,numberOfElementNodes
            CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,nodeIdx,nodeNumber,err,error,*999)
            derivativeIdx = 1
            CALL DomainElements_ElementVersionGet(dependentDomainElements,derivativeIdx,nodeIdx,elementNumber,versionIdx, &
              & err,error,*999) 
            !Get current Area values
            CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,area,err,error,*999)
            IF(area < a0Param*0.001_DP) area = a0Param*0.001_DP
            !Get material parameters
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,1,a0Param,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,2,eParam,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeNumber,3,hParam,err,error,*999)
            beta = (4.0_DP*(SQRT(PI))*eParam*hParam)/(3.0_DP*a0Param)  !(kg/m2/s2)
            !Pressure equation in mmHg
            pressure=(pExternal+beta*(SQRT(area)-SQRT(a0Param)))!/(Mref/(Lref*Tref**2.0))!*0.0075_DP
            !Update the dependent field
            IF(elementNumber<=numberOfElements) THEN
              CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & versionIdx,1,nodeNumber,1,pressure,err,error,*999)
            ENDIF
          ENDDO !nodeIdx
          !
          !-- B R A N C H   F L U X   U P W I N D I N G --
          !
          !----------------------------------------------------
          ! In order to enforce conservation of mass and momentum across discontinuous
          ! branching topologies, flux is upwinded against the conservative branch values
          ! established by the characteristic solver.
          DO nodeIdx=1,numberOfElementNodes
            CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,nodeIdx,nodeNumber,err,error,*999)
            derivativeIdx = 1
            CALL DomainNodes_DerivativeNumberOfVersionsGet(dependentDomainNodes,derivativeIdx,nodeNumber,numberOfVersions, &
              & err,error,*999)
            
            ! Find the branch node on this element
            IF(numberOfVersions>1) THEN
              elementVersionNumber=dependentDomainElements%elements(elementNumber)%elementVersions(derivativeIdx,nodeIdx)
              CALL DomainElements_ElementVersionGet(dependentDomainElements,derivativeIdx,nodeIdx,elementNumber, &
                & elementVersionNumber,err,error,*999) 
              
              ! Find the wave direction - incoming or outgoing
              DO componentIdx = 1,2
                CALL FieldVariable_ParameterSetGetLocalNode(uIndependentVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                  & derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
                IF(ABS(normalWave) > ZERO_TOLERANCE) normal = normalWave
              ENDDO !componentIdx
              
              ! Get materials parameters for node on this element
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,1,a0Param,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,2,eParam,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,3,hParam,err,error,*999)
              beta = (4.0_DP*(SQRT(PI))*eParam*hParam)/(3.0_DP*a0Param)
              
              !Get current Q & A values for node on this element
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,1,qValue,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,2,aValue,err,error,*999)
              
              !Get upwind Q & A values based on the branch (characteristics) solver
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,1,QUpwind,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,2,AUpwind,err,error,*999)
              
              ! If A goes negative during nonlinear iteration, set to positive value to avoid segfault
              IF(aValue < a0Param*0.001_DP) aValue = a0Param*0.001_DP
              
              !Momentum Equation: F_upwind - F_Current
              momentum = ((alpha*(QUpwind**2.0_DP)/AUpwind+(AUpwind**1.5_DP-a0Param**1.5_DP)*(beta/(3.0_DP*rhoParam))) &
                & - (alpha*(qValue**2.0_DP)/aValue+(aValue**1.5_DP-a0Param**1.5_DP)*(beta/(3.0_DP*rhoParam))))*normal
              
              !Continuity Equation
              mass = (QUpwind-qValue)*normal
              
              !Add momentum/mass contributions to first/last node accordingly
              IF(nodeNumber==firstNode) THEN
                residualVector%elementResidual%vector(1)= &
                  & residualVector%elementResidual%vector(1)+momentum
                residualVector%elementResidual%vector(numberOfParameters+1)= &
                  & residualVector%elementResidual%vector(numberOfParameters+1)+mass
              ELSE IF(nodeNumber==lastNode) THEN
                residualVector%elementResidual%vector(numberOfParameters)= &
                  & residualVector%elementResidual%vector(numberOfParameters)+momentum
                residualVector%elementResidual%vector(numberOfParameters*2)= &
                  & residualVector%elementResidual%vector(numberOfParameters*2)+mass
              ENDIF
            ENDIF !version>1
          ENDDO !loop nodes
          ! Update any distributed pressure field values
          CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U2_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
        ENDIF
      ENDIF
      
      ! F a c e   I n t e g r a t i o n
      IF(updateResidual) THEN
        !If specified, also perform boundary line (2D) or face (3D) integration for neumann boundary conditions
        CALL NavierStokes_FiniteElementBoundaryIntegrate(equationsSet,elementNumber,residualVariable,.FALSE.,err,error,*999)
      ENDIF

      !
      !--   A S S E M B L E   M A T R I C E S  &  V E C T O R S   --
      !
      minRowElementDOFIdx=rowElementDOFIdx
      maxRowElementDOFIdx=columnElementDOFIdx
      minColumnElementDOFIdx=rowElementDOFIdx
      maxColumnElementDOFIdx=columnElementDOFIdx
      IF(esSpecification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
        & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        IF(firstAssemblyStiffness) THEN
          IF(updateStiffness) THEN
            DO rowElementDOFIdx=minRowElementDOFIdx+1,maxRowElementDOFIdx
              DO columnElementDOFIdx=1,minColumnElementDOFIdx
                !Transpose pressure type entries for mass equation
                stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  & -stiffnessMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx)
              ENDDO !columnElementDOFIdx
            ENDDO !rowElementDOFIdx
          ENDIF
        ENDIF
      ENDIF

!!TODO: WHERE ARE THE SCALE FACTOR ADJUSTMENTS?
      
    ENDIF !update
 
    EXITS("NavierStokes_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementResidualEvaluate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices and RHS for a Navier-Stokes equation finite element equations set.
  SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnMeshComponent,componentIdx, &
      & coordinateIdx,derivativeIdx,elementVersionNumber,esSpecification(3),firstNode,gaussPointIdx,lastNode,nodeIdx,nodeNumber, &
      & numberOfColumnElementParameters,numberOfDimensions,numberOfElementNodes,numberOfGauss,numberOfParameters, &
      & numberOfResidualComponents,numberOfRowElementParameters,numberOfVersions,numberOfXi,residualVariableType, &
      & rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowMeshComponent,solverGlobalNumber,solveType,versionIdx,xiIdx
    REAL(DP) :: a0Deriv,a0Param,aDeriv,alpha,aValue,beta,columnPhi,dColumnPhidXi(3),dRowPhidXi(3),dXidX(3,3),eDeriv,eParam, &
      & gaussWeight,hDeriv,hParam,jacobian,jacobianGaussWeight,kappa,mass,momentum1,momentum2,muParam,muScale,normal,normalWave, &
      & qDeriv,qValue,rhoParam,rowPhi,sum,uValue(3),uDeriv(3,3),wValue(3)
    LOGICAL  :: updateJacobian,updateResidual,updateRHS
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,independentBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition,independentDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,independentDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements, &
      & independentDomainElements,rowDomainElements
    TYPE(DomainNodesType), POINTER :: dependentDomainNodes
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology, &
      & independentDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,geometricField,materialsField,independentField
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,independentInterpPoint, &
      & uMaterialsInterpPoint,vMaterialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters, &
      & independentInterpParameters,uMaterialsInterpParameters,vMaterialsInterpParameters
    TYPE(FieldVariableType), POINTER :: geometricVariable,independentVariable,residualVariable,uMaterialsVariable, &
      & vMaterialsVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,quadratureScheme,rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_FiniteElementJacobianEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
      
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,1,jacobianMatrix,err,error,*999)
    CALL JacobianMatrix_UpdateMatrixGet(jacobianMatrix,updateJacobian,err,error,*999)

    IF(updateJacobian) THEN

      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(residualVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,residualVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(residualVariable,residualVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(residualVariable,numberOfResidualComponents,err,error,*999)

      dXidX=0.0_DP
      uValue=0.0_DP
      uDeriv=0.0_DP
      wValue=0.0_DP

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)

      NULLIFY(dependentField)
      CALL EquationsInterpolation_DependentFieldGet(equationsInterpolation,dependentField,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentDomainNodes)
      CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,residualVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,residualVariableType,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters, &
        & err,error,*999)

      NULLIFY(geometricField)
      CALL EquationsInterpolation_GeometricFieldGet(equationsInterpolation,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
      NULLIFY(geometricDomain)
      CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
      NULLIFY(geometricDomainTopology)
      CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
      NULLIFY(geometricDomainElements)
      CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
      NULLIFY(geometricBasis)
      CALL DomainElements_ElementBasisGet(geometricDomainElements,elementNumber,geometricBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfXi,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
        & err,error,*999)

      NULLIFY(materialsField)
      CALL EquationsInterpolation_MaterialsFieldGet(equationsInterpolation,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      NULLIFY(uMaterialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & uMaterialsInterpParameters,err,error,*999)
      NULLIFY(uMaterialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uMaterialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,uMaterialsInterpParameters, &
        & err,error,*999)        

      NULLIFY(vMaterialsVariable)
      NULLIFY(vMaterialsInterpParameters)
      NULLIFY(vMaterialsInterpPoint)
      NULLIFY(equationsSetField)
      NULLIFY(independentField)
      NULLIFY(independentVariable)
      NULLIFY(independentDecomposition)
      NULLIFY(independentDomain)
      NULLIFY(independentDomainTopology)
      NULLIFY(independentDomainElements)
      NULLIFY(independentBasis)
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)

      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
        !Do nothing
      CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE)
      CASE(EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE)
      CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE,EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & vMaterialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vMaterialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vMaterialsInterpParameters, &
          & err,error,*999)        
      CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)
        CALL EquationsInterpolation_IndependentFieldGet(equationsInterpolation,independentField,err,error,*999)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
        CALL Field_DecompositionGet(independentField,independentDecomposition,err,error,*999)
        CALL Decomposition_DomainGet(independentDecomposition,0,independentDomain,err,error,*999)
        CALL Domain_DomainTopologyGet(independentDomain,independentDomainTopology,err,error,*999)
        CALL DomainTopology_DomainElementsGet(independentDomainTopology,independentDomainElements,err,error,*999)
        CALL DomainElements_ElementBasisGet(independentDomainElements,elementNumber,independentBasis,err,error,*999)
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,independentInterpPoint, &
          & err,error,*999)        
        CALL Field_InterpolationParametersElementGet(FIELD_MESH_VELOCITY_SET_TYPE,elementNumber,independentInterpParameters, &
          & err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes fluid type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT

      !Loop over all Gauss points
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      DO gaussPointIdx=1,numberOfGauss

        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%numberOfXi,geometricInterpPointMetrics,err,error,*999)
        
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,uMaterialsInterpPoint, &
          & err,error,*999)

        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        DO xiIdx=1,numberOfXi
          DO coordinateIdx=1,numberOfDimensions
            dXidX(xiIdx,coordinateIdx)=geometricInterpPointMetrics%dXidX(xiIdx,coordinateIdx)
          ENDDO !coordinateIdx
        ENDDO !xiIdx

        !Get ALE velocity, w 
        IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
            & err,error,*999)
          wValue(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ELSE
          wValue(1:numberOfDimensions)=0.0_DP
        ENDIF

        !Get the constitutive law (non-Newtonian) viscosity based on shear rate
        IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
          !Note the constant from the U_VARIABLE is a scale factor
          muScale = uMaterialsInterpPoint%values(1,NO_PART_DERIV)
          !Get the gauss point based value returned from the CellML solver
          CALL FieldVariable_ParameterSetGetLocalGaussPoint(vMaterialsVariable,FIELD_VALUES_SET_TYPE,gaussPointIdx, &
            & elementNumber,1,muParam,err,error,*999)
          muParam=muParam*muScale
        ELSE
          muParam = uMaterialsInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        rhoParam = uMaterialsInterpPoint%values(2,NO_PART_DERIV)

        IF(esSpecification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
          & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN

          uValue(1:numberOfDimensions)=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
          DO xiIdx=1,numberOfXi
            uDeriv(1:numberOfDimensions,xiIdx)= &
              & dependentInterpPoint%values(1:numberOfDimensions,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
          ENDDO !xiIdx
          !Start with calculation of partial matrices
          !Here wValues must be ZERO if ALE part of linear matrix
          wValue=0.0_DP
        ENDIF

        IF(esSpecification(3)==EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE.OR.  &
          & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          !Loop over field components
          rowElementDOFIdx=0
          DO rowComponentIdx=1,numberOfDimensions
            CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
            NULLIFY(rowDomain)
            CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
            NULLIFY(rowDomainTopology)
            CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
            NULLIFY(rowDomainElements)
            CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
            NULLIFY(rowBasis)
            CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
            CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
            NULLIFY(rowQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
            jacobianGaussWeight=jacobian*gaussWeight
            DO rowElementParameterIdx=1,numberOfRowElementParameters
              rowElementDOFIdx=rowElementDOFIdx+1
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussPointIdx,rowPhi,err,error,*999)
              DO xiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                  & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,dRowPhidXi(xiIdx),err,error,*999)
              ENDDO !xiIdx
              columnElementDOFIdx=0
              !Loop over element columns
              DO columnComponentIdx=1,numberOfResidualComponents
                CALL FieldVariable_ComponentMeshComponentGet(residualVariable,columnComponentIdx,columnMeshComponent, &
                  & err,error,*999)
                NULLIFY(columnDomain)
                CALL Decomposition_DomainGet(dependentDecomposition,columnMeshComponent,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                NULLIFY(columnQuadratureScheme)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                  & err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  !Calculate some general values needed below
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                    & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  sum=0.0_DP
                  DO xiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,dColumnPhidXi(xiIdx),err,error,*999)
                    !Calculate J1 only
                    sum=sum+(columnPhi*uDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,columnComponentIdx)*rowPhi*rhoParam)
                  ENDDO !xiIdx
                  !Calculate matrix
                  jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  IF(columnComponentIdx==rowComponentIdx) THEN
                    !Calculate J2 only
                    sum=0.0_DP
                    !Calculate sum
                    DO componentIdx=1,numberOfDimensions
                      DO xiIdx=1,numberOfXi
                        sum=sum+rhoParam*(uValue(componentIdx)-wValue(componentIdx))*dColumnPhidXi(xiIdx)* &
                          & dXidX(xiIdx,componentIdx)*rowPhi
                      ENDDO !xiIdx
                    ENDDO !componentIdx
                    !Calculate matrix
                    jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  ENDIF !column component = row component
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
        ENDIF

        !------------------------------------------------------------------
        ! R e s i d u a l - b a s e d    S t a b i l i s a t i o n
        !------------------------------------------------------------------
        IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE) THEN
          CALL NavierStokes_ResidualBasedStabilisation(equationsSet,elementNumber,gaussPointIdx,muParam,rhoParam,.TRUE., &
            & err,error,*999)
        ENDIF
        
        !------------------------------------------------------------------
        ! 1 D  T r a n s i e n t
        !------------------------------------------------------------------
        
        !Start with Matrix Calculations
        IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
          & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
          qValue=dependentInterpPoint%values(1,NO_PART_DERIV)
          qDeriv=dependentInterpPoint%values(1,FIRST_PART_DERIV)
          aValue=dependentInterpPoint%values(2,NO_PART_DERIV)
          aDeriv=dependentInterpPoint%values(2,FIRST_PART_DERIV)
          CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,3,alpha,err,error,*999)
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,vMaterialsInterpPoint, &
            & err,error,*999)
          a0Param=vMaterialsInterpPoint%values(1,NO_PART_DERIV)
          a0Deriv=vMaterialsInterpPoint%values(1,FIRST_PART_DERIV)
          eParam=vMaterialsInterpPoint%values(2,NO_PART_DERIV)
          eDeriv=vMaterialsInterpPoint%values(2,FIRST_PART_DERIV)
          hParam=vMaterialsInterpPoint%values(3,NO_PART_DERIV)
          hDeriv=vMaterialsInterpPoint%values(3,FIRST_PART_DERIV)
          beta = (4.0_DP*(SQRT(PI))*eParam*hParam)/(3.0_DP*a0Param)  !(kg/m2/s2)
          kappa = 8.0_DP*PI*muParam/rhoParam ! viscous resistance operator
          
          ! If A goes negative during nonlinear iteration, give ZERO_TOLERANCE value to avoid segfault
          IF(aValue < a0Param*0.001_DP) aValue = a0Param*0.001_DP
          
          rowElementDOFIdx=0
          !Loop Over Element Rows
          DO rowComponentIdx=1,numberOfResidualComponents
            CALL FieldVariable_ComponentMeshComponentGet(residualVariable,rowComponentIdx,rowMeshComponent,err,error,*999)
            NULLIFY(rowDomain)
            CALL Decomposition_DomainGet(dependentDecomposition,rowMeshComponent,rowDomain,err,error,*999)
            NULLIFY(rowDomainTopology)
            CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
            NULLIFY(rowDomainElements)
            CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
            NULLIFY(rowBasis)
            CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
            CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
            NULLIFY(rowQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
            dXidX(1,1)=0.0_DP
            !Calculate dxi_dx in 3D
            DO xiIdx=1,numberOfXi
              DO coordinateIdx=1,numberOfDimensions
                dXidX(1,1)=dXidX(1,1)+(geometricInterpPointMetrics%dXidX(xiIdx,coordinateIdx))**2.0_DP
              ENDDO !coordinateIdx
            ENDDO !xiIdx
            dXidX(1,1)=SQRT(dXidX(1,1))
            !Loop Over Element rows
            DO rowElementParameterIdx=1,numberOfRowElementParameters
              rowElementDOFIdx=rowElementDOFIdx+1
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussPointIdx,rowPhi,err,error,*999)
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & FIRST_PART_DERIV,gaussPointIdx,dRowPhidXi(1),err,error,*999)
              columnElementDOFIdx=0
              !Loop Over Element Columns
              DO columnComponentIdx=1,numberOfResidualComponents
                CALL FieldVariable_ComponentMeshComponentGet(residualVariable,columnComponentIdx,columnMeshComponent, &
                  & err,error,*999)
                NULLIFY(columnDomain)
                CALL Decomposition_DomainGet(dependentDecomposition,columnMeshComponent,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                NULLIFY(columnQuadratureScheme)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                  & err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                    & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                    & FIRST_PART_DERIV,gaussPointIdx,dColumnPhidXi(1),err,error,*999)
                  
                  !Momentum Equation (dF/dQ)
                  IF(rowComponentIdx==1.AND.columnComponentIdx==1) THEN
                    sum=((alpha*2.0_DP*columnPhi*qDeriv/aValue +  &
                      & alpha*2.0_DP*qValue*dColumnPhidXi(1)/aValue+ &
                      & (-2.0_DP)*alpha*qValue*columnPhi*aDeriv/(aValue**2.0_DP))*dXidX(1,1)+ &   !Convective
                      & ((columnPhi*kappa/aValue)))*rowPhi                                        !Viscosity
                    jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  ENDIF
                  
                  !Momentum Equation (dF/dA)
                  IF(rowComponentIdx==1.AND.columnComponentIdx==2) THEN
                    sum=((((-2.0_DP*alpha*qValue*columnPhi*qDeriv)/(aValue**2.0_DP))+ &
                      & ((2.0_DP*alpha*columnPhi*(qValue**2.0_DP)*aDeriv)/(aValue**3.0_DP))+ &
                      & (-alpha*((qValue/aValue)**2.0_DP)*dColumnPhidXi(1))+ &                              !Convective
                      & ((0.5_DP*columnPhi*(1.0_DP/SQRT(aValue))*aDeriv+SQRT(aValue)*dColumnPhidXi(1))+ &   !Area Gradient
                      & ((1.0_DP/SQRT(a0Param))-((3.0_DP/(a0Param))*SQRT(aValue)))*(a0Deriv) + &            !Ref Area Gradient
                      & (2.0_DP*columnPhi*1.5_DP*SQRT(aValue))*hDeriv/hParam+ &                             !Thickness Gradient
                      & (2.0_DP*columnPhi*1.5_DP*SQRT(aValue))*eDeriv/eParam) &                             !Elasticity Gradient
                      & *beta/(2.0_DP*rhoParam))*dXidX(1,1)+(-columnPhi*kappa*qValue/aValue**2.0_DP))*rowPhi!Viscosity
                    jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  ENDIF
                  
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDDO !rowElementParameterIdx
          ENDDO !rowComponentIdx
        ENDIF
      ENDDO !gaussPointIdx

      !B o u n d a r y   I n t e g r a t i o n
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
        & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
        ! Calculate the Jacobian of the nonlinear boundary stabilisation term if beta > 0
        CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,beta,err,error,*999)
        IF(beta > ZERO_TOLERANCE) &
          & CALL NavierStokes_FiniteElementBoundaryIntegrate(equationsSet,elementNumber,residualVariable,.TRUE.,err,error,*999)
      CASE DEFAULT
        ! Do nothing for other equation set subtypes
      END SELECT
      
      IF(esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE .OR. &
        & esSpecification(3)==EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE) THEN
        CALL Basis_NumberOfLocalNodesGet(dependentBasis,numberOfElementNodes,err,error,*999)
        CALL DomainElements_MaxElementParametersGet(dependentDomainElements,numberOfParameters,err,error,*999)
        CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,1,firstNode,err,error,*999)
        CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,numberOfParameters,lastNode,err,error,*999)
        !Get material constants
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,rhoParam,err,error,*999)
        !
        !-- B R A N C H   F L U X   U P W I N D I N G --
        !
        !----------------------------------------------------
        ! In order to enforce conservation of mass and momentum across discontinuous
        ! branching topologies, flux is upwinded against the conservative branch values
        ! established by the characteristic solver.
        DO nodeIdx=1,numberOfElementNodes
          CALL DomainElements_ElementNodeGet(dependentDomainElements,elementNumber,1,firstNode,err,error,*999)
          derivativeIdx = 1
          CALL DomainNodes_DerivativeNumberOfVersionsGet(dependentDomainNodes,derivativeIdx,nodeIdx,numberOfVersions, &
            & err,error,*999) 
          !Find the branch node on this element
          IF(numberOfVersions>1) THEN
            CALL DomainElements_ElementVersionGet(dependentDomainElements,derivativeIdx,nodeIdx,elementNumber,versionIdx, &
              & err,error,*999) 
            
            !Find the wave direction - incoming or outgoing
            DO componentIdx = 1,2
              CALL FieldVariable_ParameterSetGetLocalNode(independentVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
                & derivativeIdx,nodeNumber,componentIdx,normalWave,err,error,*999)
              IF(ABS(normalWave) > ZERO_TOLERANCE) normal = normalWave
            ENDDO !componentIdx
            
            !Get materials parameters for node on this element
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
              & derivativeIdx,nodeNumber,1,a0Param,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
              & derivativeIdx,nodeNumber,2,eParam,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
              & derivativeIdx,nodeNumber,3,hParam,err,error,*999)
            beta = (4.0_DP*(SQRT(PI))*eParam*hParam)/(3.0_DP*a0Param)
            
            !Get current Q & A values
            CALL FieldVariable_ParameterSetGetLocalNode(residualVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
              & derivativeIdx,nodeNumber,1,qValue,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(residualVariable,FIELD_VALUES_SET_TYPE,elementVersionNumber, &
              & derivativeIdx,nodeNumber,2,aValue,err,error,*999)
            
            !Momentum Equation, d/dQ
            momentum1 = (-alpha*2.0_DP*qValue/aValue)*normal
            
            !Momentum Equation, d/dA
            momentum2 = (alpha*(qValue/aValue)**2.0_DP-1.5_DP*(aValue**0.5_DP)*(beta/(3.0_DP*rhoParam)))*normal
            
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
            ENDIF
          ENDIF !version>1
        ENDDO !loop nodes

      ENDIF
      
    ENDIF !update Jacobian
    
    EXITS("NavierStokes_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_FiniteElementJacobianEvaluate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the Navier-Stokes problem post solve.
  SUBROUTINE NavierStokes_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<The solver for the post solve
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,equationsSetNumber,esSpecification(3),inputIteration,iterationNumber,loopLevel,loopType, &
      & maxIterations,numberOfEquationsSets,outputIteration,pSpecification(3),solverGlobalNumber,solveType,solveType2, &
      & subLoopIndex,timestep
    REAL(DP) :: absoluteTolerance,currentTime,relativeTolerance,startTime,stopTime,timeIncrement
    LOGICAL :: continueLoop,convergedFlag,fluidEquationsSetFound
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DynamicSolverType), POINTER :: dynamicSolver
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldParameterSetType), POINTER :: upwindParameterSet
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_PostSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
      CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
      CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==2) CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE)
      CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
      SELECT CASE(solveType)
      CASE(SOLVER_NONLINEAR_TYPE)
        ! Characteristic solver- copy branch Q,A values to new parameter set
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
        CALL ControlLoop_CurrentWhileInformationGet(controlLoop,iterationNumber,maxIterations,absoluteTolerance, &
          & relativeTolerance,continueLoop,err,error,*999)
        IF(iterationNumber == 1) &
          & CALL FieldVariable_ParameterSetsCopy(dependentVariable,FIELD_VALUES_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP, &
          & err,error,*999)
      CASE(SOLVER_DYNAMIC_TYPE)
        !Navier-Stokes solver: do nothing
      CASE DEFAULT
        localError="The solver type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
          & " is invalid for a 1D Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        SELECT CASE(solverGlobalNumber)
        CASE(1)
          !Characteristic solver- copy branch Q,A values to new parameter set
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetsCopy(dependentVariable,FIELD_VALUES_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP, &
            & err,error,*999)          
        CASE(2)
          !1D Navier-Stokes solver
          CALL ControlLoop_LoopLevelGet(controlLoop,loopLevel,err,error,*999)
          IF(loopLevel==3) THEN
            ! check characteristic/ N-S convergence at branches
            !CALL NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
            & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_SIMPLE_TYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        IF(solverGlobalNumber==1) THEN
          ! DAE solver- do nothing
        ELSE
          localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
            & " is invalid for the CellML DAE simple loop of a 1D-0D coupled Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The solver control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for the a 1D-0D coupled Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      SELECT CASE(solverGlobalNumber)
      CASE(1)
        !Characteristic solver- copy branch Q,A values to new parameter set
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)        
        CALL FieldVariable_ParameterSetsCopy(dependentVariable,FIELD_VALUES_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP, &
          & err,error,*999)
      CASE(2)
        !Check characteristic/ N-S convergence at branches
        !CALL NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*999)
      CASE(3)
        !Advection solver output data if necessary
        CALL ControlLoop_CurrentWhileInformationGet(controlLoop,iterationNumber,maxIterations,absoluteTolerance, &
          & relativeTolerance,continueLoop,err,error,*999)
        IF(.NOT.continueLoop) THEN
          !1D NSE solver output data if N-S/Chars converged
          CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
          & " is invalid for a 1D Navier-Stokes and Advection problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        SELECT CASE(solverGlobalNumber)
        CASE(1)
          !Characteristic solver- copy branch Q,A values to new parameter set
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetsCopy(dependentVariable,FIELD_VALUES_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP, &
            & err,error,*999)
        CASE(2)
          !1D Navier-Stokes solver
          CALL ControlLoop_LoopLevelGet(controlLoop,loopLevel,err,error,*999)
          IF(controlLoop%controlLoopLevel==3) THEN
            !Check characteristic/ N-S convergence at branches
            !CALL NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
            & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_SIMPLE_TYPE)
        !DAE and advection solvers - output data if post advection solve
        CALL ControlLoop_SubLoopIndexGet(controlLoop,subLoopIndex,err,error,*999)
        IF(subLoopIndex==3) CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
      CASE DEFAULT
        localError="The solver control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for the a 1D-0D coupled Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==2) THEN
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
          & timestep,outputIteration,inputIteration,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsExists(solver,solverEquations,err,error,*999)
        IF(ASSOCIATED(solverEquations)) THEN
          convergedFlag = .FALSE.
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)          
          CALL NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*999)
        ENDIF
        CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
      ENDIF
    CASE(PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==2) THEN
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          !If this is a coupled constitutive (non-Newtonian) viscosity problem, update shear rate values
          !to be passed to the CellML solver at beginning of next timestep
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) &
            & CALL NavierStokes_ShearRateCalculate(equationsSet,err,error,*999)
        ENDDO !equationsSetIdx
      ENDIF
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        SELECT CASE(solverGlobalNumber)
        CASE(1)
          !Characteristic solver- copy branch Q,A values to new parameter set
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetsCopy(dependentVariable,FIELD_VALUES_SET_TYPE,FIELD_UPWIND_VALUES_SET_TYPE,1.0_DP, &
            & err,error,*999)          
        CASE(2)
          !1D Navier-Stokes solver
          CALL ControlLoop_LoopLevelGet(controlLoop,loopLevel,err,error,*999)
         IF(loopLevel==3) THEN
            ! check characteristic/ N-S convergence at branches
            !CALL NavierStokes_CoupleCharacteristics(controlLoop,solver,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
            & " is invalid for the iterative 1D-0D coupled Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(CONTROL_SIMPLE_TYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        IF(solverGlobalNumber == 1) THEN
          !DAE solver- do nothing
        ELSE
          localError="The solver global number of "//TRIM(NumberToVString(solverGlobalNumber,"*",err,error))// &
            & " is invalid for the CellML DAE simple loop of a 1D-0D coupled Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The solver control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for the a 1D-0D coupled Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
      CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
      SELECT CASE(solveType)
      CASE(SOLVER_LINEAR_TYPE)
        !Post solve for the linear solver
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solver2)
        CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
        CALL Solver_SolverTypeGet(solver2,solveType2,err,error,*999)
        IF(solveType2/=SOLVER_DYNAMIC_TYPE) CALL FlagError("Dynamic solver is not associated for ALE problem.",err,error,*999)
        NULLIFY(dynamicSolver)
        CALL Solver_DynamicSolverGet(solver2,dynamicSolver,err,error,*999)
        dynamicSolver%ale=.TRUE.
      CASE(SOLVER_DYNAMIC_TYPE)
        !Post solve for the dynamic solver
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        equationsSetIdx=1
        fluidEquationsSetFound=.FALSE.
        DO WHILE(.NOT.fluidEquationsSetFound.AND.equationsSetIdx<=numberOfEquationsSets)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
          IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
            & .AND.esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
            & .AND.(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
            & .OR.esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
            fluidEquationsSetFound=.TRUE.
          ELSE
            equationsSetIdx=equationsSetIdx+1
          ENDIF
        ENDDO !equations set search
        IF(.NOT.fluidEquationsSetFound) CALL FlagError("Fluid equations set not found.",err,error,*999)
        !CALL NavierStokes_WallShearStressCalculate(equationsSet,err,error,*999)
        CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
      CASE DEFAULT
        !Do nothing
      END SELECT
    CASE DEFAULT
      localError="The third problem specification of  "// &
        & TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_PostSolve")
    RETURN
999 ERRORSEXITS("NavierStokes_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_PostSolve

  !
  !================================================================================================================================
  !

  !>Update boundary conditions for Navier-Stokes flow pre solve
  SUBROUTINE NavierStokes_PreSolveUpdateBoundaryConditions(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<The solver to update the boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionCheckVariable,componentIdx,componentNumberVelocity,coordinateIdx, &
      & computationNode,coupledNodeNumber,currentTimeLoopIteration,dependentDof,dependentVariableType,derivativeIdx, &
      & dimensionIdx,elementIdx,elementNodeIdx,equationsSetIdx,equationsSetIndex,esSpecification(3),fluidNodeNumber, &
      & fluidVariableType,globalDerivIndex,globalDof,I,independentDof,independentVariableType,inputIteration,interpolationType, &
      & J,K,localDOF,localNodeNumber,loopOutputType,maxNumberOfElementParameters,nodeIdx,nodeNumber,numberOfDependentComponents, &
      & numberOfDimensions,numberOfEquationsSets,numberOfFittedNodes,numberOfGlobalNodes,numberOfInterfaceNodes,numberOfNodes, &
      & numberOfNodeDerivatives,numberOfNodesXIC(3),numberOfSourceTimesteps,numberOfVariables,numberOfVersions, &
      & outputIterationNumber,outputType,pSpecification(3),searchIdx,solidNodeNumber,solidVariableType,solverGlobalNumber, &
      & solveType,timeIdx,userNodeNumber,versionIdx,variableIdx,velocityInterpolationType
    REAL(DP) :: componentValues(3),currentTime,displacementValue,fluidGFValue,Lref,Mref,muParam,newLaplaceBoundaryValue,QP,QPP, &
      & rhoParam,solidDFValue,startTime,stopTime,tCoordinates(20,3),timeData,timeIncrement,Tref,VALUE,x(3),xiCoordinates(3)
    REAL(DP), POINTER :: analyticParameters(:),boundaryValues(:),geometricParameters(:),materialsParameters(:), &
      & meshVelocityValues(:),nodeData(:,:),qSpline(:),qValues(:),tValues(:)
    REAL(DP), ALLOCATABLE :: fittedNodes(:)
    LOGICAL :: fluidEquationsSetFound=.FALSE.,ghostNode,importDataFromFile,nodeExists,parameterSetCreated, &
      & solidEquationsSetFound=.FALSE.,solidNodeFound=.FALSE.
    CHARACTER(70) :: inputFile,tempString
    TYPE(BasisType), POINTER :: basis
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsSetType), POINTER :: equationsSet,solidEquationsSet,fluidEquationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField,fluidDependentField
    TYPE(FieldType), POINTER :: independentField,solidDependentField,fluidGeometricField
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,fluidDependentVariable,geometricVariable, &
      & independentVariable,materialsVariable,uMaterialsVariable,solidDependentVariable
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations,solidSolverEquations,fluidSolverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping,solidSolverMapping,fluidSolverMapping
    TYPE(SolverType), POINTER :: solver2
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup 

    ENTERS("NavierStokes_PreSolveUpdateBoundaryConditions",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
      & currentTimeLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(1))      
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      SELECT CASE(pSpecification(3))
      CASE(PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
        ! do nothing ???
      CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
        CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
        IF(solverGlobalNumber==2) THEN
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)          
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          !Fitting boundary condition- get values from file
          !TODO: this should be generalised with input filenames specified from the example file when IO is improved
          NULLIFY(equationsIndependent)
          CALL EquationsSet_IndependentExists(equationsSet,equationsIndependent,err,error,*999)
          IF(ASSOCIATED(equationsIndependent)) THEN
            !Read in field values to independent field
            NULLIFY(independentField)
            CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
            NULLIFY(independentVariable)
            CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
            NULLIFY(boundaryConditionsVariable)
            CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable, err,error,*999)
            !Read in field data from file
            !Loop over nodes and update independent field values - if a fixed fitted boundary, also update dependent
            componentNumberVelocity = 1
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
            numberOfDimensions = numberOfDependentComponents - 1
            !Get the nodes on this computation domain
            CALL FieldVariable_ComponentInterpolationGet(independentVariable,componentNumberVelocity, &
              & velocityInterpolationType,err,error,*999)
            IF(velocityInterpolationType==FIELD_NODE_BASED_INTERPOLATION) THEN
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(independentVariable,componentNumberVelocity,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              CALL DomainNodes_NumberOfGlobalNodesGet(domainNodes,numberOfGlobalNodes,err,error,*999)

              !Construct the filename based on the computation node and time step
              inputFile = './../interpolatedData/fitData'//TRIM(NumberToVString(currentTimeLoopIteration,"*",err,error))//'.dat'

              INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
              IF(importDataFromFile) THEN
                !Read fitted data from input file (if exists)
                CALL ControlLoop_OutputTypeGet(controlLoop,loopOutputType,err,error,*999)
                IF(loopOutputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                  & CALL WriteString(GENERAL_OUTPUT_TYPE,"Updating independent field and boundary nodes from "//inputFile, &
                  & err,error,*999)
                OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                READ(10,*) numberOfFittedNodes
                ALLOCATE(fittedNodes(numberOfFittedNodes))
                READ(10,*) fittedNodes
                DO nodeIdx=1, numberOfFittedNodes
                  userNodeNumber=INT(fittedNodes(nodeIdx),INTG)
                  CALL DomainNodes_NodeCheckExists(domainNodes,userNodeNumber,nodeExists,localNodeNumber,ghostNode,err,error,*999)
                  IF(nodeExists.AND..NOT.ghostNode) THEN
                    ! Node found on this computation node
                    READ(10,*) (componentValues(componentIdx), componentIdx=1,numberOfDimensions)
                    DO componentIdx=1,numberOfDimensions
                      CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,1,localNodeNumber,componentIdx,dependentDOF, &
                        & err,error,*999)
                      CALL FieldVariable_LocalNodeDOFGet(independentVariable,1,1,localNodeNumber,componentIdx,independentDOF, &
                        & err,error,*999)
                      VALUE = componentValues(componentIdx)
                      CALL FieldVariable_ParameterSetUpdateLocalDOF(independentVariable,FIELD_VALUES_SET_TYPE,independentDof, &
                        & VALUE,err,error,*999)
                      CALL FieldVariable_ComponentDOFGetUserNode(dependentVariable,1,1,userNodeNumber,componentIdx,localDof, &
                        & globalDof,err,error,*999)
                      CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDOF, &
                        & boundaryConditionCheckVariable,err,error,*999)
                      IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_FITTED) &
                        & CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDof, &
                        & VALUE,err,error,*999)
                    ENDDO !componentIdx
                  ELSE
                    !Dummy read if this node not on this computation node
                    READ(10,*)
                  ENDIF
                ENDDO !nodeIdx
                DEALLOCATE(fittedNodes)
                CLOSE(UNIT=10)
                ! Update any distributed field values
                CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateStart(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateFinish(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
              ENDIF !check import file exists
            ENDIF
          ENDIF !Equations set independent

          !Analytic equations
          NULLIFY(equationsAnalytic)
          CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
          IF(ASSOCIATED(equationsAnalytic)) THEN
            CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
            !Standard analytic functions
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
              !Update analytic time value with current time
              CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
              !Calculate analytic values
              IF(ASSOCIATED(boundaryConditions)) &
                & CALL NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
            CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN, &
              & EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
              !Geometric parameters
              NULLIFY(geometricParameters)
              CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
              ! Analytic parameters
              NULLIFY(analyticField)
              CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
              NULLIFY(analyticVariable)
              NULLIFY(analyticParameters)
              IF(ASSOCIATED(analyticField)) THEN
                CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
                CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
              ENDIF
              !Materials parameters
              NULLIFY(materialsField)
              CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
              NULLIFY(materialsVariable)
              NULLIFY(materialsParameters)
              NULLIFY(equationsMaterials)
              IF(ASSOCIATED(materialsField)) THEN
                CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)                
                CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
                CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
              ENDIF
              NULLIFY(equations)
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              NULLIFY(equationsInterpolation)
              CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
              NULLIFY(dependentInterpParameters)
              CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
                & dependentInterpParameters,err,error,*999)
              NULLIFY(dependentInterpPoint)
              CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
                & err,error,*999)
              CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
              DO variableIdx=1,numberOfVariables
                NULLIFY(dependentVariable)
                CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
                CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
                DO componentIdx=1,numberOfDependentComponents
                  CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,interpolationType,err,error,*999)
                  IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
                    & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                  NULLIFY(domain)
                  CALL FieldVariable_ComponentDomainGet(independentVariable,componentIdx,domain,err,error,*999)
                  NULLIFY(domainTopology)
                  CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
                  NULLIFY(domainNodes)
                  CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
                  CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
                  NULLIFY(domainElements)
                  CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
                  CALL DomainElements_MaxElementParametersGet(domainElements,maxNumberOfElementParameters,err,error,*999)
                  !Should be replaced by boundary node flag
                  DO nodeIdx=1,numberOfNodes
                    CALL DomainNodes_NodeSurroundingElementGet(domainNodes,1,nodeIdx,elementIdx,err,error,*999)
                    NULLIFY(basis)
                    CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
                    CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,dependentInterpParameters, &
                      & err,error,*999)
                    elementNodeIdx=0
                    xiCoordinates=0.0_DP
                    numberOfNodesXiC=1
                    CALL Basis_NumberOfNodesXiCGet(basis,numberOfNodesXIC,err,error,*999)
                    !\todo: change definitions as soon as adjacent elements / boundary elements calculation works for simplex
                    IF(maxNumberOfElementParameters==4 .OR. &
                      & maxNumberOfElementParameters==9 .OR. &
                      & maxNumberOfElementParameters==16 .OR. &
                      & maxNumberOfElementParameters==8 .OR. &
                      & maxNumberOfElementParameters==27 .OR. &
                      & maxNumberOfElementParameters==64) THEN
                      DO K=1,numberOfNodesXIC(3)
                        DO J=1,numberOfNodesXIC(2)
                          DO I=1,numberOfNodesXIC(1)
                            elementNodeIdx=elementNodeIdx+1
                            CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,nodeNumber, &
                              & err,error,*999)
                            IF(nodeNumber==nodeIdx) EXIT
                            xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXIC(1)-1))
                          END DO !I
                          CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,nodeNumber,err,error,*999)
                          IF(nodeNumber==nodeIdx) EXIT
                          xiCoordinates(1)=0.0_DP
                          xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXIC(2)-1))
                        END DO !J
                        CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,nodeNumber,err,error,*999)
                        IF(nodeNumber==nodeIdx) EXIT
                        xiCoordinates(1)=0.0_DP
                        xiCoordinates(2)=0.0_DP
                        IF(numberOfNodesXIC(3)/=1) xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXIC(3)-1))
                      END DO !K
                      CALL Field_InterpolateXi(NO_PART_DERIV,xiCoordinates,dependentInterpPoint,err,error,*999)
                    ELSE
                      !\todo: Use boundary flag
                      IF(maxNumberOfElementParameters==3) THEN
                        tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                        tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                        tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                      ELSE IF(maxNumberOfElementParameters==6) THEN
                        tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                        tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                        tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                        tCoordinates(4,1:2)=[0.5_DP,0.5_DP]
                        tCoordinates(5,1:2)=[1.0_DP,0.5_DP]
                        tCoordinates(6,1:2)=[0.5_DP,1.0_DP]
                      ELSE IF(maxNumberOfElementParameters==10.AND.numberOfDimensions==2) THEN
                        tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
                        tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
                        tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
                        tCoordinates(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                        tCoordinates(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
                        tCoordinates(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                      ELSE IF(maxNumberOfElementParameters==4) THEN
                        tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                        tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                        tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                      ELSE IF(maxNumberOfElementParameters==10.AND.numberOfDimensions==3) THEN
                        tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                        tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                        tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
                        tCoordinates(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
                        tCoordinates(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
                        tCoordinates(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
                        tCoordinates(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
                        tCoordinates(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
                      ELSE IF(maxNumberOfElementParameters==20) THEN
                        tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
                        tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
                        tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
                        tCoordinates(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
                        tCoordinates(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
                        tCoordinates(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
                        tCoordinates(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
                        tCoordinates(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
                        tCoordinates(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
                      ENDIF
                      DO K=1,maxNumberOfElementParameters
                        CALL DomainElements_ElementNodeGet(domainElements,K,elementIdx,nodeNumber,err,error,*999)
                        IF(nodeNumber==nodeIdx) EXIT
                      ENDDO !K
                      IF(numberOfDimensions==2) THEN
                        CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(K,1:2),dependentInterpPoint,err,error,*999)
                      ELSE IF(numberOfDimensions==3) THEN
                        CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(K,1:3),dependentInterpPoint,err,error,*999)
                      ENDIF
                    ENDIF
                    X=0.0_DP
                    X(1:numberOfDimensions)=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
                    !Loop over the derivatives
                    CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                    DO derivativeIdx=1,numberOfNodeDerivatives
                      CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivIndex, &
                        & err,error,*999)
                      !Define muParam, density=1
                      muParam=materialsParameters(1)
                      !Define rhoParam, density=2
                      rhoParam=materialsParameters(2)
                      CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,X, &
                        & currentTime,dependentVariableType,globalDerivIndex,componentIdx, &
                        & numberOfDimensions,dependentVariable%numberOfComponents, &
                        & analyticParameters,materialsParameters,VALUE,err,error,*999)
                      !Default to version 1 of each node derivative
                      CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                        & err,error,*999)
                      CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOF, &
                        & VALUE,err,error,*999)
                      NULLIFY(boundaryConditionsVariable)
                      CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                        & err,error,*999)
                      CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                        & boundaryConditionCheckVariable,err,error,*999)
                      IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) THEN
                        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF,VALUE, &
                          & err,error,*999)
                      ENDIF
                    ENDDO !derivativeIdx
                  ENDDO !nodeIdx
                ENDDO !componentIdx
                CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
              ENDDO !variableIdx
              CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid or not implemented."
              CALL FlagError(localError,err,error,*999)
            END SELECT !Standard/unit analytic subtypes           
          ENDIF ! Analytic boundary conditions

          !TODO implement non-analytic time-varying boundary conditions (i.e. from file)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        ENDIF !solver number 2

      CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
        !If analytic flow waveform, calculate and update
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)          
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          SELECT CASE(analyticFunctionType)
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
            CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
            ! Calculate analytic values
            CALL NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
            ! Perform spline interpolation of values from a file
            CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
            CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
            DO variableIdx=1,numberOfVariables
              NULLIFY(dependentVariable)
              CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
              NULLIFY(boundaryConditionsVariable)
              CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                & err,error,*999)
              IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
                DO componentIdx=1,numberOfDependentComponents
                  CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,interpolationType, &
                    & err,error,*999)
                  IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
                    & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                  NULLIFY(domain)
                  CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
                  NULLIFY(domainTopology)
                  CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
                  NULLIFY(domainNodes)
                  CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
                  ! Create the analytic field values type on the dependent field if it does not exist
                  CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  !Loop over the local nodes excluding the ghosts.
                  CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
                  DO nodeIdx=1,numberOfNodes
                    CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
                    CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                    DO derivativeIdx=1,numberOfNodeDerivatives
                      CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions, &
                        & err,error,*999)
                      DO versionIdx=1,numberOfVersions
                        CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                          & dependentDOF,err,error,*999)
                        !Update dependent field value if this is a splint BC
                         CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dependentDOF, &
                          & boundaryConditionCheckVariable,err,error,*999)
                        IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                          !Update analytic field if file exists and dependent field if boundary condition set
                          inputFile = './input/interpolatedData/1D/'
                          IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) inputFile = TRIM(inputFile) // 'U/component'
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
                              READ(10,*) (nodeData(timeIdx,dimensionIdx), dimensionIdx=1,2)
                            ENDDO !timeIdx
                            CLOSE(UNIT=10)
                            tValues = nodeData(:,1)
                            qValues = nodeData(:,2)
                            CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues,2,0.0_DP,2,0.0_DP,qSpline, &
                              & err,error,*999)
                            CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues,qSpline,currentTime,VALUE,QP,QPP, &
                              & err,error,*999)
                            DEALLOCATE(nodeData)
                            DEALLOCATE(qSpline)
                            DEALLOCATE(qValues)
                            DEALLOCATE(tValues)
                            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,dependentDOF, &
                              & VALUE,err,error,*999)
                            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                              & dependentDOF,VALUE,err,error,*999)
                          ENDIF
                        ENDIF ! check if import data file exists
                      ENDDO !versionIdx
                    ENDDO !derivativeIdx
                  ENDDO !nodeIdx
                  !Update distributed field values
                  CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            ENDDO !variableIdx
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_HEART)
            !Using heart lumped parameter model for input
            CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
            NULLIFY(materialsField)
            CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
            NULLIFY(uMaterialsVariable)
            CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
            CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,5,Lref,err,error,*999)
            CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,6,Tref,err,error,*999)
            CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,7,Mref,err,error,*999)
            CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
            DO variableIdx=1,numberOfVariables
              NULLIFY(dependentVariable)
              CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
              NULLIFY(boundaryConditionsVariable)
              CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                & err,error,*999)
              IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
                DO componentIdx=1,numberOfDependentComponents
                  CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,interpolationType,err,error,*999)
                  IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
                    & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                  NULLIFY(domain)
                  CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
                  NULLIFY(domainTopology)
                  CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
                  NULLIFY(domainNodes)
                  CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
                  !Loop over the local nodes excluding the ghosts.
                  CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
                  DO nodeIdx=1,numberOfNodes
                    CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
                    CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                    DO derivativeIdx=1,numberOfNodeDerivatives
                      CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions, &
                        & err,error,*999)
                      DO versionIdx=1,numberOfVersions
                        CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                          & dependentDOF,err,error,*999)
                        !Update dependent field value if this is a splint BC
                        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dependentDOF, &
                          & boundaryConditionCheckVariable,err,error,*999)
                        IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                          CALL Field_ParameterSetGetLocalNode(dependentField,FIELD_U1_VARIABLE_TYPE, &
                            & FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx,userNodeNumber,1,VALUE, &
                            & err,error,*999)
                          !Convert Q from ml/s to non-dimensionalised form.
                          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,dependentDOF, &
                            & ((Lref**3.0)/Tref)*VALUE,err,error,*999)
                        ENDIF
                      ENDDO !versionIdx
                    ENDDO !derivativeIdx
                  ENDDO !nodeIdx
                  !Update distributed field values
                  CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            ENDDO !variableIdx
          CASE DEFAULT
            !Do nothing (might have another use for analytic equations)
          END SELECT
        ENDIF ! Check for analytic equations
        !Update any multiscale boundary values (coupled 0D or non-reflecting)
        CALL NavierStokes_UpdateMultiscaleBoundary(equationsSet,boundaryConditions,timeIncrement,err,error,*999)
      CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
        & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
        !TODO: this should be set up so it uses the same pre_solve steps as the individual 3D/1D equations sets
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        SELECT CASE(solveType)
        CASE(SOLVER_DYNAMIC_TYPE)
          ! --- D y n a m i c    S o l v e r s ---
          CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
            & currentTimeLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
          DO equationsSetIdx=1,numberOfEquationsSets
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
            CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
            SELECT CASE(esSpecification(3))
              ! --- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
            CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
              !Fitting boundary condition- get values from file
              !TODO: this should be generalised when IO is improved
              NULLIFY(equationsIndependent)
              CALL EquationsSet_IndependentExists(equationsSet,equationsIndependent,err,error,*999)
              IF(ASSOCIATED(equationsIndependent)) THEN
                !Read in field values to independent field
                NULLIFY(dependentField)
                CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                NULLIFY(dependentVariable)
                CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
                NULLIFY(independentField)
                CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
                NULLIFY(independentVariable)
                CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
                NULLIFY(boundaryConditionsVariable)
                CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                  & err,error,*999)
                !Read in field data from file
                !Loop over nodes and update independent field values. If a fixed fitted boundary, also update dependent
                componentNumberVelocity = 1
                !Get the nodes on this computation domain
                CALL FieldVariable_ComponentInterpolationGet(independentVariable,componentNumberVelocity,interpolationType, &
                  & err,error,*999)
                IF(interpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
                  & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                NULLIFY(domain)
                CALL FieldVariable_ComponentDomainGet(independentVariable,componentNumberVelocity,domain,err,error,*999)
                NULLIFY(domainTopology)
                CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
                NULLIFY(domainNodes)
                CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
                CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
                CALL DomainNodes_NumberOfGlobalNodesGet(domainNodes,numberOfGlobalNodes,err,error,*999)
                
                ! Construct the filename based on the computation node and time step
                inputFile='./../interpolatedData/fitData'//TRIM(NumberToVString(currentTimeLoopIteration,"*",err,error))//'.dat'

                INQUIRE(FILE=inputFile, EXIST=importDataFromFile)
                IF(importDataFromFile) THEN
                  !Read fitted data from input file (if exists)
                  NULLIFY(workGroup)
                  CALL Solver_WorkGroupGet(solver,workGroup,err,error,*999)
                  CALL WorkGroup_GroupNodeNumberGet(workGroup,computationNode,err,error,*999)
                  CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
                  IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                    IF(computationNode==0) THEN
                      CALL WriteString(GENERAL_OUTPUT_TYPE,"Updating independent field and boundary nodes from " &
                        & //inputFile,err,error,*999)
                    ENDIF
                  ENDIF
                  OPEN(UNIT=10, FILE=inputFile, STATUS='OLD')
                  READ(10,*) numberOfFittedNodes
                  ALLOCATE(fittedNodes(numberOfFittedNodes))
                  READ(10,*) fittedNodes
                  DO nodeIdx=1,numberOfFittedNodes
                    userNodeNumber=INT(fittedNodes(nodeIdx),INTG)
                    CALL DomainNodes_NodeCheckExists(domainNodes,userNodeNumber,nodeExists,localNodeNumber,ghostNode, &
                      & err,error,*999)
                    IF(nodeExists.AND..NOT.ghostNode) THEN
                      !Node found on this computation node
                      READ(10,*) (componentValues(componentIdx), componentIdx=1,numberOfDimensions)
                      DO componentIdx=1,numberOfDimensions
                        CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,1,localNodeNumber,componentIdx,dependentDOF, &
                          & err,error,*999)
                        CALL FieldVariable_LocalNodeDOFGet(independentVariable,1,1,localNodeNumber,componentIdx,independentDOF, &
                          & err,error,*999)
                        VALUE = componentValues(componentIdx)
                        CALL FieldVariable_ParameterSetUpdateLocalDOF(independentVariable,FIELD_VALUES_SET_TYPE,independentDOF, &
                          & VALUE,err,error,*999)
                        CALL FieldVariable_ComponentDOFGetUserNode(dependentVariable,1,1,userNodeNumber,componentIdx,localDOF, &
                          & globalDof,err,error,*999)
                        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDOF, &
                          & boundaryConditionCheckVariable,err,error,*999)
                        IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF,VALUE, &
                            & err,error,*999)
                        ENDIF
                      ENDDO !componentIdx
                    ELSE
                      !Dummy read if this node not on this computation node
                      READ(10,*)
                    ENDIF
                  ENDDO !nodeIdx
                  DEALLOCATE(fittedNodes)
                  CLOSE(UNIT=10)
                  !Update any distributed field values
                  CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateStart(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateFinish(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                ENDIF !check import file exists
              ENDIF !Equations set independent
              !Analytic equations
              NULLIFY(equationsAnalytic)
              CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
              IF(ASSOCIATED(equationsAnalytic)) THEN
                CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
                !Standard analytic functions
                IF(analyticFunctionType/=EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                  localError="Analytic equations type "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                    & " is not yet implemented for a 3D Navier-Stokes equations set for a multiscale problem."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                ! Update analytic time value with current time
                CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
                !Calculate analytic values
                CALL NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions, err,error,*999)
             ENDIF
              ! --- 1 D    N a v i e r - S t o k e s   E q u a t i o n s ---
            CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
              & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
              !If analytic flow waveform, calculate and update
              NULLIFY(equationsAnalytic)
              CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
              IF(ASSOCIATED(equationsAnalytic)) THEN
                CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
                SELECT CASE(analyticFunctionType)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
                  & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
                  CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
                  !Calculate analytic values
                  CALL NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
                CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SPLINT_FROM_FILE)
                  !Perform spline interpolation of values from a file
                  CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
                  NULLIFY(dependentField)
                  CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                  CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
                  DO variableIdx=1,numberOfVariables
                    NULLIFY(dependentVariable)
                    CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
                    NULLIFY(boundaryConditionsVariable)
                    CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                      & err,error,*999)
                    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
                      DO componentIdx=1,numberOfDependentComponents
                        CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,interpolationType, &
                          & err,error,*999)
                        IF(interpolationType==FIELD_NODE_BASED_INTERPOLATION) &
                          & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                        NULLIFY(domain)
                        CALL FieldVariable_ComponentDomainGet(independentVariable,componentNumberVelocity,domain,err,error,*999)
                        NULLIFY(domainTopology)
                        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
                        NULLIFY(domainNodes)
                        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
                        !Loop over the local nodes excluding the ghosts.
                        DO nodeIdx=1,numberOfNodes
                          CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
                          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                          DO derivativeIdx=1,numberOfNodeDerivatives
                            CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions, &
                              & err,error,*999)
                            DO versionIdx=1,numberOfVersions
                              !Update analytic field if file exists and !the dependent field if boundary condition set
                              inputFile = './input/interpolatedData/1D/'
                              IF(dependentVariableType == FIELD_U_VARIABLE_TYPE) inputFile = TRIM(inputFile) // 'U/component'
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
                                !Create the analytic field values type on the dependent field if it does not exist
                                CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE, &
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
                                  READ(10,*) (nodeData(timeIdx,dimensionIdx), dimensionIdx=1,2)
                                ENDDO !timeIdx
                                CLOSE(UNIT=10)
                                tValues = nodeData(:,1)
                                qValues = nodeData(:,2)
                                CALL spline_cubic_set(numberOfSourceTimesteps,tValues,qValues,2,0.0_DP,2,0.0_DP,qSpline, &
                                  & err,error,*999)
                                CALL spline_cubic_val(numberOfSourceTimesteps,tValues,qValues,qSpline,currentTime,VALUE,QP,QPP, &
                                  & err,error,*999)
                                
                                DEALLOCATE(nodeData)
                                DEALLOCATE(qSpline)
                                DEALLOCATE(qValues)
                                DEALLOCATE(tValues)
                                
                                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx, &
                                  & componentIdx,dependentDOF,err,error,*999)
                                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                                  & dependentDof,VALUE,err,error,*999)
                                !Update dependent field value if this is a splint BC
                                CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dependentDOF, &
                                  & boundaryConditionCheckVariable,err,error,*999)
                                IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_FITTED) THEN
                                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE, &
                                    & dependentDOF,VALUE,err,error,*999)
                                ENDIF
                              ENDIF ! check if import data file exists
                            ENDDO !versionIdx
                          ENDDO !derivativeIdx
                        ENDDO !nodeIdx
                        !Update distributed field values
                        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
                        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                     ENDDO !componentIdx
                    ENDIF
                  ENDDO !variableIdx
                CASE DEFAULT
                  !Do nothing (might have another use for analytic equations)
                END SELECT
              ENDIF ! Check for analytic equations
              !Update any multiscale boundary values (coupled 0D or non-reflecting)
              CALL NavierStokes_UpdateMultiscaleBoundary(equationsSet,boundaryConditions,timeIncrement,err,error,*999)
            CASE DEFAULT
              localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                & " is not valid for a multiscale dynamic Navier-Stokes solver."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !equationsSetIdx
        CASE DEFAULT
          localError="The solve type of "//TRIM(NumberToVString(solver%solveType,"*",err,error))// &
            & " is invalid for pre_solve_update_BC step of a multiscale Navier-Stokes problem type."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE(PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE)
        !Pre solve for the linear solver
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(boundaryConditions)
        CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfCOmponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
        NULLIFY(boundaryValues)
        CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_NONLINEAR_TYPE,boundaryValues, &
          & numberOfDimensions,BOUNDARY_CONDITION_FIXED_INLET,controlLoop%timeLoop%inputNumber, &
          & currentTimeLoopIteration,currentTime,1.0_DP,err,error,*999)
        CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
          DO componentIdx=1,numberOfDependentComponents
            NULLIFY(domain)
            CALL FieldVariable_ComponentDomainGet(independentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF,err,error,*999)
                CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                  & boundaryConditionCheckVariable,err,error,*999)
                IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) &
                  & CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                  & boundaryValues(localDOF),err,error,*999)
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !variableIdx
        !\todo: This part should be read in out of a file eventually
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)          
        !Pre solve for the linear solver
        IF(solveType==SOLVER_LINEAR_TYPE) THEN
           IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
            & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          CALL Field_NumberOfCOmponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          NULLIFY(independentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
          NULLIFY(independentVariable)
          CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
          NULLIFY(boundaryValues)
          CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,boundaryValues, &
            & numberOfDimensions,BOUNDARY_CONDITION_MOVED_WALL,controlLoop%timeLoop%inputNumber, &
            & currentTimeLoopIteration,currentTime,1.0_DP,err,error,*999)
          CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
            DO componentIdx=1,numberOfDependentComponents
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                    & err,error,*999)
                  CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                    & boundaryConditionCheckVariable,err,error,*999)
                  IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                      & boundaryValues(localDOF),err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
          ENDDO !variableIdx
          CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
          !\todo: This part should be read in out of a file eventually
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Pre solve for the dynamic solver
        ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
          IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
            & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          NULLIFY(independentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
          NULLIFY(independentVariable)
          CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
          NULLIFY(meshVelocityValues)
          CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
            & err,error,*999)
          NULLIFY(boundaryValues)
          CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
          CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS(SOLVER_LINEAR_TYPE,boundaryValues, &
            & numberOfDimensions,BOUNDARY_CONDITION_FIXED_INLET,controlLoop%timeLoop%inputNumber, &
            & currentTimeLoopIteration,currentTime,1.0_DP,err,error,*999)
          CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
            DO componentIdx=1,numberOfDependentComponents
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                    & err,error,*999)
                  displacementValue=0.0_DP
                  CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                    & boundaryConditionCheckVariable,err,error,*999)
                  IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                      & meshVelocityValues(localDOF),err,error,*999)
                  ELSE IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                      & boundaryValues(localDOF),err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
          ENDDO !variableIdx
          CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
            & err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
          NULLIFY(dependentVariable)
          CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        ENDIF
        ! do nothing ???
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
          CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
          CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
          !Pre solve for the linear solver
          IF(solveType==SOLVER_LINEAR_TYPE) THEN
            IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
              & CALL WriteString(GENERAL_OUTPUT_TYPE,"Mesh movement change boundary conditions... ",err,error,*999)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            NULLIFY(boundaryConditions)
            CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
            NULLIFY(solverMapping)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
            NULLIFY(geometricField)
            CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
            NULLIFY(geometricVariable)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
            CALL FieldVariable_NumberOfCOmponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
            NULLIFY(dependentField)
            CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            NULLIFY(dependentVariable)
            CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
            NULLIFY(independentField)
            CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
            NULLIFY(boundaryConditionsVariable)
            CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
            NULLIFY(solvers)
            CALL Solver_SolversGet(solver,solvers,err,error,*999)
            !Update moving wall nodes from solid/fluid gap (as we solve for displacements of the mesh
            !in Laplacian smoothing step).
            IF(pSpecification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
              & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
              NULLIFY(solver2)
              CALL Solvers_SolverGet(solvers,2,solver2,err,error,*999)
            ELSE
              CALL Solvers_SolverGet(solvers,4,solver2,err,error,*999)
            ENDIF
            !Find the FiniteElasticity equations set as there is a NavierStokes equations set too
            NULLIFY(solidSolverEquations)
            CALL Solver_SolverEquationsGet(solver2,solidSolverEquations,err,error,*999)
            NULLIFY(solidSolverMapping)
            CALL SolverEquations_SolverMappingGet(solidSolverEquations,solidSolverMapping,err,error,*999)
            CALL SolverMapping_NumberOfEquationsSetsGet(solidSolverMapping,numberOfEquationsSets,err,error,*999)
            equationsSetIndex=1
            solidEquationsSetFound=.FALSE.
            DO WHILE(.NOT.solidEquationsSetFound.AND.equationsSetIndex<=numberOfEquationsSets)
              NULLIFY(solidEquationsSet)
              CALL SolverMapping_EquationsSetGet(solidSolverMapping,equationsSetIndex,solidEquationsSet,err,error,*999)
              CALL EquationsSet_SpecificationGet(solidEquationsSet,3,esSpecification,err,error,*999)
              IF(esSpecification(1)==EQUATIONS_SET_ELASTICITY_CLASS &
                & .AND.esSpecification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE &
                & .AND.((esSpecification(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_DYNAMIC_ST_VENANT_KIRCHOFF_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_DYNAMIC_MOONEY_RIVLIN_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_DYNAMIC_COMP_ST_VENANT_KIRCHOFF_SUBTYPE) &
                & .OR.(esSpecification(3)==EQUATIONS_SET_DYNAMIC_COMP_MOONEY_RIVLIN_SUBTYPE))) THEN
                solidEquationsSetFound=.TRUE.
              ELSE
                equationsSetIndex=equationsSetIndex+1
              ENDIF
            ENDDO
            IF(.NOT.solidEquationsSetFound) &
              & CALL FlagError("Solid equations set not found when trying to update boundary conditions.",err,error,*999)
            NULLIFY(solidDependentField)
            CALL EquationsSet_DependentFieldGet(solidEquationsSet,solidDependentField,err,error,*999)
            NULLIFY(solidDependentVariable)
            CALL Field_VariableIndexGet(solidDependentField,1,solidDependentVariable,solidVariableType,err,error,*999)              
            !Find the NavierStokes equations set as there is a FiniteElasticity equations set too
            NULLIFY(fluidSolverEquations)
            CALL Solver_SolverEquationsGet(solver2,fluidSolverEquations,err,error,*999)
            NULLIFY(fluidSolverMapping)
            CALL SolverEquations_SolverMappingGet(fluidSolverEquations,fluidSolverMapping,err,error,*999)
            CALL SolverMapping_NumberOfEquationsSetsGet(fluidSolverMapping,numberOfEquationsSets,err,error,*999)
            equationsSetIndex=1
            fluidEquationsSetFound=.FALSE.
            DO WHILE (.NOT.fluidEquationsSetFound.AND.equationsSetIndex<=numberOfEquationsSets)
              NULLIFY(fluidEquationsSet)
              CALL SolverMapping_EquationsSetGet(fluidSolverMapping,equationsSetIndex,fluidEquationsSet,err,error,*999)
              CALL EquationsSet_SpecificationGet(fluidEquationsSet,3,esSpecification,err,error,*999)
              IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
                & .AND.esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
                & .AND.(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
                & .OR.esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
                fluidEquationsSetFound=.TRUE.
              ELSE
                equationsSetIndex=equationsSetIndex+1
              ENDIF
            END DO
            IF(.NOT.fluidEquationsSetFound) &
              & CALL FlagError("Fluid equations set not found when trying to update boundary conditions.",err,error,*999)
            NULLIFY(fluidGeometricField)
            CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
            NULLIFY(fluidDependentField)
            CALL EquationsSet_DependentFieldGet(fluidEquationsSet,fluidDependentField,err,error,*999)
            NULLIFY(fluidDependentVariable)
            CALL Field_VariableIndexGet(fluidDependentField,1,fluidDependentVariable,fluidVariableType,err,error,*999)
            NULLIFY(interfaceCondition)
            CALL SolverMapping_InterfaceConditionGet(fluidSolverMapping,1,interfaceCondition,err,error,*999)
            NULLIFY(INTERFACE)
            CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
            NULLIFY(meshConnectivity)
            CALL Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*999)
            CALL InterfaceMeshConnectivity_NumberOfInterfaceNodesGet(meshConnectivity,numberOfInterfaceNodes,err,error,*999)
            variableIdx=1
            CALL FieldVariable_NumberOfComponentsGet(fluidDependentVariable,numberOfDependentComponents,err,error,*999)
            DO componentIdx=1,numberOfDependentComponents
              NULLIFY(domain)
              CALL FieldVariable_ComponentDomainGet(fluidDependentVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                    & err,error,*999)
                  CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                    & boundaryConditionCheckVariable,err,error,*999)
                  !Update moved wall nodes only
                  IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_MOVED_WALL) THEN
                    !NOTE: assuming same mesh and mesh nodes for fluid domain and moving mesh domain
                    fluidNodeNumber=nodeIdx
                    DO searchIdx=1,numberOfInterfaceNodes
                      CALL InterfaceMeshConnectivity_CoupledNodeNumberGet(meshConnectivity,2,searchIdx,coupledNodeNumber, &
                        & err,error,*999)
                      IF(coupledNodeNumber==nodeIdx) THEN
                        CALL InterfaceMeshConnectivity_CoupledNodeNumberGet(meshConnectivity,1,searchIdx,solidNodeNumber, &
                          & err,error,*999)
                        solidNodeFound=.TRUE.
                      ENDIF
                    ENDDO !searchIdx
                    IF(.NOT.solidNodeFound.OR.fluidNodeNumber==0) CALL FlagError("Solid interface node not found.",err,error,*999)
                    !Default to version number 1
                    IF(variableIdx==1) THEN
                      CALL Field_ParameterSetGetNode(fluidGeometricField,fluidVariableType,FIELD_VALUES_SET_TYPE,&
                        & 1,derivativeIdx,fluidNodeNumber,componentIdx,fluidGFValue,Err,Error,*999)
                    ELSE
                      fluidGFValue=0.0_DP
                    ENDIF
                    CALL Field_ParameterSetGetNode(solidDependentField,solidVariableType,FIELD_VALUES_SET_TYPE, &
                      & 1,derivativeIdx,solidNodeNumber,componentIdx,solidDFValue,Err,Error,*999)
                    newLaplaceBoundaryValue=solidDFValue-fluidGFValue
                    CALL Field_ParameterSetUpdateLocalDOF(dependentField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,localDOF,newLaplaceBoundaryValue,err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
            CALL Field_ParameterSetUpdateStart(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & err,error,*999)
            CALL Field_ParameterSetUpdateFinish(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & err,error,*999)
            !Pre solve for the dynamic solver
          ELSE IF(solveType==SOLVER_DYNAMIC_TYPE) THEN
            !IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
            !  & CALL WriteString(GENERAL_OUTPUT_TYPE,"Velocity field change boundary conditions... ",err,error,*999)
            !NULLIFY(solverEquations)
            !CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !NULLIFY(solverMapping)
            !CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            !CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
            !!Find the NavierStokes equations set as there is a finite elasticity equations set too
            !equationsSetIndex=1
            !ALENavierStokesEquationsSetFound=.FALSE.
            !DO WHILE(.NOT.ALENavierStokesEquationsSetFound.AND.equationsSetIndex<=numberOfEquationsSets)
            !  NULLIFY(equationsSet)
            !  CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIndex,equationsSet,err,error,*999)
            !  CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
            !  IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
            !    & .AND.esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
            !    & .AND.esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
            !    & .AND.esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
            !    ALENavierStokesEquationsSetFound=.TRUE.
            !  ELSE
            !    equationsSetIndex=equationsSetIndex+1
            !  ENDIF
            !ENDDO !equations set search
            !IF(.NOT.ALENavierStokesEquationsSetFound) &
            !  & CALL FlagError("ALE NavierStokes equations set not found when trying to update boundary conditions.", &
            !  & err,error,*999)
            !!Get boundary conditions
            !NULLIFY(boundaryConditions)
            !CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
            !NULLIFY(geometricField)
            !CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
            !NULLIFY(geometricVariable)
            !CALL Field_VariableGet(geometric,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
            !CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
            !NULLIFY(dependentField)
            !CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
            !NULLIFY(dependentVariable)
            !CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
            !NULLIFY(independentField)
            !CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
            !NULLIFY(independentVariable)
            !CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            !NULLIFY(boundaryConditionsVariable)
            !CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
            !NULLIFY(meshVelocityValues)
            !CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
            !  & err,error,*999)
            !NULLIFY(boundaryValues)
            !CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
            !!Get update for time-dependent boundary conditions
            !CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
            !  & currentTimeLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
            !IF(inputIteration==1) THEN
            !  componentBC=1
            !  CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(geometricField,solveType,inletNodes, &
            !    & boundaryValues,BOUNDARY_CONDITION_FIXED_INLET,inputIteration,currentTime,stopTime,err,error,*999)
            !  DO nodeIdx=1,SIZE(inletNodes)
            !    CALL FieldVariable_ParameterSetUpdateNode(dependentVariable,FIELD_VALUES_SET_TYPE,1,1,inletNodes(nodeIdx), &
            !      & componentBC,boundaryValues(nodeIdx),err,error,*999)
            !  ENDDO !nodeIdx
            !ELSE
            !  !Figure out which component we're applying BC at
            !  IF(iterationNumber==2) THEN
            !    componentBC=1
            !  ELSE
            !    componentBC=2
            !  ENDIF
            !  !Get inlet nodes and the corresponding velocities
            !  CALL FluidMechanics_IO_UpdateBoundaryConditionUpdateNodes(geometricField,solveType,inletNodes,boundaryValues, &
            !    & BOUNDARY_CONDITION_FIXED_INLET,inputIteration,currentTime,stopTime,err,error,*999)
            !  DO nodeIdx=1,SIZE(inletNodes)
            !    CALL Field_ParameterSetUpdateNode(dependentVariable,FIELD_VALUES_SET_TYPE,1,1,inletNodes(nodeIdx), &
            !      & componentBC,boundaryValues(nodeIdx),err,error,*999)
            !  ENDDO !nodeIdx
            !ENDIF
            !CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_MESH_VELOCITY_SET_TYPE,meshVelocityValues, &
            !  & err,error,*999)
            !CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_BOUNDARY_SET_TYPE,boundaryValues,err,error,*999)
            !CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
            !CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          ENDIF
          ! do nothing ???
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is not valid for a FiniteElasticity-NavierStokes problem type of a multi physics problem class."
          CALL FlagError(localError,Err,Error,*999)
        END SELECT
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for NavierStokes_PreSolve of a multi physics problem class."
        CALL FlagError(localError,Err,Error,*999)
      END SELECT
    CASE DEFAULT
      localError="The first problem specification of "//TRIM(NumberToVString(pSpecification(1),"*",err,error))// &
        & " is not valid for NavierStokes_PreSolveUpdateBoundaryConditions."
      CALL FlagError(localError,Err,Error,*999)
    END SELECT

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
  SUBROUTINE NavierStokes_PreSolveALEUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: componentIdx,currentIteration,derivativeIdx,equationsSetIdx,esSpecification(3),fluidNode, &
      & fluidNumberOfDimensions,geometricMeshComponent,inputIteration,interfaceNumberOfDimensions,interpolationType, &
      & inputType,inputOption,laplaceNumberOfDimensions,localDOF,nodeIdx,numberOfComponents,numberOfDependentComponents, &
      & numberOfEquationsSets,numberOfNodes,numberOfNodeDerivatives,numberOfVariables,numberOfVersions,outputIteration, &
      & pSpecification(3),solidNode,solveType,totalNumberOfNodes,variableIdx,variableType,versionIdx
    REAL(DP) :: alpha,currentTime,previousSolidNodePosition,solidDelta,solidNodePosition,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: meshDisplacementValues(:)
    LOGICAL :: fluidEquationsSetFound=.FALSE.
    LOGICAL :: solidEquationsSetFound=.FALSE.
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes   
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: laplaceEquationsSet,fluidEquationsSet,solidEquationsSet
    TYPE(FieldType), POINTER :: laplaceDependentField,laplaceGeometricField,fluidIndependentField,fluidGeometricField, &
      & solidDependentField,interfaceGeometricField
    TYPE(FieldVariableType), POINTER :: fluidGeometricVariable,fluidIndependentVariable,interfaceGeometricVariable, &
      & laplaceDependentVariable,laplaceGeometricVariable,solidDependentVariable
    TYPE(InterfaceType), POINTER :: fsiInterface
    TYPE(InterfaceConditionType), POINTER :: fsiInterfaceCondition
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: laplaceSolverEquations,fluidSolverEquations,fsiSolverEquations
    TYPE(SolverMappingType), POINTER :: laplaceSolverMapping,fsiSolverMapping,fluidSolverMapping
    TYPE(SolverType), POINTER :: dynamicSolver,laplaceSolver
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_PreSolveALEUpdateMesh",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime, &
      & currentIteration,outputIteration,inputIteration,err,error,*999)    

    SELECT CASE(pSpecification(1))
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      SELECT CASE(pSpecification(3))
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
      CASE(PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
        !Update mesh within the dynamic solver
        IF(solveType/=SOLVER_DYNAMIC_TYPE) CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
        IF(.NOT.solver%dynamicSolver%ale) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
        !Get the dependent field for the three component Laplace problem
        NULLIFY(laplaceSolver)
        CALL Solvers_SolverGet(solvers,1,laplaceSolver,err,error,*999)
        NULLIFY(laplaceSolverEquations)
        CALL Solver_SolverEquationsGet(laplaceSolver,laplaceSolverEquations,err,error,*999)
        NULLIFY(laplaceSolverMapping)
        CALL SolverEquations_SolverMappingGet(laplaceSolverEquations,laplaceSolverMapping,err,error,*999)
        NULLIFY(laplaceEquationsSet)
        CALL SolverMapping_EquationsSetGet(laplaceSolverMapping,1,laplaceEquationsSet,err,error,*999)
        NULLIFY(laplaceDependentField)
        CALL EquationsSet_DependentFieldGet(laplaceEquationsSet,laplaceDependentField,err,error,*999)
        NULLIFY(laplaceGeometricField)
        CALL EquationsSet_GeometricFieldGet(laplaceEquationsSet,laplaceGeometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(laplaceGeometricField,FIELD_U_VARIABLE_TYPE,laplaceNumberOfDimensions,err,error,*999)
        !Get the independent field for the ALE Navier-Stokes problem
        NULLIFY(dynamicSolver)
        CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
        NULLIFY(fluidSolverEquations)
        CALL Solver_SolverEquationsGet(dynamicSolver,fluidSolverEquations,err,error,*999)
        NULLIFY(fluidSolverMapping)
        CALL SolverEquations_SolverMappingGet(fluidSolverEquations,fluidSolverMapping,err,error,*999)
        NULLIFY(fluidEquationsSet)
        CALL SolverMapping_EquationsSetGet(fluidSolverMapping,1,fluidEquationsSet,err,error,*999)
        NULLIFY(fluidGeometricField)
        CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
        NULLIFY(fluidGeometricVariable)
        CALL Field_VariableGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidGeometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(fluidGeometricVariable,fluidNumberOfDimensions,err,error,*999)
        NULLIFY(fluidIndependentField)
        CALL EquationsSet_IndependentFieldGet(fluidEquationsSet,fluidIndependentField,err,error,*999)
        NULLIFY(fluidIndependentVariable)
        CALL Field_VariableGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,fluidIndependentVariable,err,error,*999)
        !Copy result from Laplace mesh movement to Navier-Stokes' independent field
        IF(fluidNumberOfDimensions/=laplaceNumberOfDimensions) &
          & CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",err,error,*999)
        DO componentIdx=1,fluidNumberOfDimensions
          CALL Field_ParametersToFieldParametersCopy(laplaceDependentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & componentIdx,fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,componentIdx, &
            & err,error,*999)
        ENDDO !componentIdx
        !Use calculated values to update mesh
        CALL FieldVariable_ComponentMeshComponentGet(fluidGeometricVariable,1,geometricMeshComponent,err,error,*999)
        NULLIFY(meshDisplacementValues)
        CALL Field_ParameterSetDataGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & meshDisplacementValues,err,error,*999)
        DO componentIdx=1,fluidNumberOfDimensions
          NULLIFY(domain)
          CALL FieldVariable_DomainGet(fluidGeometricVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(fluidGeometricVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                & err,error,*999)
              CALL FieldVariable_ParameterSetAddLocalDOF(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                & meshDisplacementValues(localDOF),err,error,*999)
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !componentIdx
        CALL FieldVariable_ParameterSetDataRestore(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & meshDisplacementValues,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        !Now use displacement values to calculate velocity values
        alpha=1.0_DP/timeIncrement
        CALL FieldVariable_ParameterSetsCopy(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
          & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
          !Update mesh within the dynamic solver
          IF(solveType/=SOLVER_DYNAMIC_TYPE) CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
          IF(.NOT.solver%dynamicSolver%ale) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
          !Get the dependent field for the Laplace problem
          IF(pSpecification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            NULLIFY(laplaceSolver)
            CALL Solvers_SolverGet(solvers,3,laplaceSolver,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            NULLIFY(laplaceSolver)
            CALL Solvers_SolverGet(solvers,5,laplaceSolver,err,error,*999)
          ELSE
            NULLIFY(laplaceSolver)
            CALL Solvers_SolverGet(solvers,4,laplaceSolver,err,error,*999)
          ENDIF
          NULLIFY(laplaceSolverEquations)
          CALL Solver_SolverEquationsGet(laplaceSolver,laplaceSolverEquations,err,error,*999)
          NULLIFY(laplaceSolverMapping)
          CALL SolverEquations_SolverMappingGet(laplaceSolverEquations,laplaceSolverMapping,err,error,*999)
          NULLIFY(laplaceEquationsSet)
          CALL SolverMapping_EquationsSetGet(laplaceSolverMapping,1,laplaceEquationsSet,err,error,*999)
          NULLIFY(laplaceGeometricField)
          CALL EquationsSet_GeometricFieldGet(laplaceEquationsSet,laplaceGeometricField,err,error,*999)
          NULLIFY(laplaceGeometricVariable)
          CALL Field_VariableGet(laplaceGeometricField,FIELD_U_VARIABLE_TYPE,laplaceGeometricVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(laplaceGeometricVariable,laplaceNumberOfDimensions,err,error,*999)
          NULLIFY(laplaceDependentField)
          CALL EquationsSet_DependentFieldGet(laplaceEquationsSet,laplaceDependentField,err,error,*999)
          NULLIFY(laplaceDependentVariable)
          CALL Field_VariableGet(laplaceDependentField,FIELD_U_VARIABLE_TYPE,laplaceDependentVariable,err,error,*999)
          !Get the independent field for the ALE Navier-Stokes problem
          !Get the dynamic solver
          IF(pSpecification(3)==PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            NULLIFY(dynamicSolver)
            CALL Solvers_SolverGet(solvers,2,dynamicSolver,err,error,*999)
          ELSE IF(pSpecification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE.OR. &
            & pSpecification(3)==PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE) THEN
            NULLIFY(dynamicSolver)
            CALL Solvers_SolverGet(solvers,4,dynamicSolver,err,error,*999)
          ELSE
            NULLIFY(dynamicSolver)
            CALL Solvers_SolverGet(solvers,3,dynamicSolver,err,error,*999)
          ENDIF
          NULLIFY(fsiSolverEquations)
          CALL Solver_SolverEquationsGet(dynamicSolver,fsiSolverEquations,err,error,*999)
          NULLIFY(fsiSolverMapping)
          CALL SolverEquations_SolverMappingGet(fsiSolverEquations,fsiSolverMapping,err,error,*999)
          CALL SolverMapping_NumberOfEquationsSetsGet(fsiSolverMapping,numberOfEquationsSets,err,error,*999)
          !Find the Navier Stokes equations set
          equationsSetIdx=1
          fluidEquationsSetFound=.FALSE.
          DO WHILE(.NOT.fluidEquationsSetFound.AND.equationsSetIdx<=numberOfEquationsSets)
            NULLIFY(fluidEquationsSet)
            CALL SolverMapping_EquationsSetGet(fsiSolverMapping,equationsSetIdx,fluidEquationsSet,err,error,*999)
            CALL EquationsSet_SpecificationGet(fluidEquationsSet,3,esSpecification,err,error,*999)
            IF(esSpecification(1)==EQUATIONS_SET_FLUID_MECHANICS_CLASS &
              & .AND.esSpecification(2)==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TYPE &
              & .AND.(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE &
              & .OR.esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)) THEN
              fluidEquationsSetFound=.TRUE.
            ELSE
              equationsSetIdx=equationsSetIdx+1
            ENDIF
          ENDDO
          IF(.NOT.fluidEquationsSetFound) &
            & CALL FlagError("ALE NavierStokes equations set not found when trying to update ALE mesh.",err,error,*999)
          NULLIFY(fluidGeometricField)
          CALL EquationsSet_GeometricFieldGet(fluidEquationsSet,fluidGeometricField,err,error,*999)
          NULLIFY(fluidGeometricVariable)
          CALL Field_VariableGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidGeometricVariable,err,error,*999)
          NULLIFY(fluidIndependentField)
          CALL EquationsSet_IndependentFieldGet(fluidEquationsSet,fluidIndependentField,err,error,*999)
          NULLIFY(fluidIndependentVariable)
          CALL Field_VariableGet(fluidIndependentField,FIELD_U_VARIABLE_TYPE,fluidIndependentVariable,err,error,*999)
          CALL Field_NumberOfComponentsGet(fluidGeometricField,FIELD_U_VARIABLE_TYPE,fluidNumberOfDimensions, &
            & err,error,*999)
          !Copy result from Laplace mesh movement to Navier-Stokes' independent field
          IF(fluidNumberOfDimensions/=laplaceNumberOfDimensions) &
            & CALL FlagError("Dimension of Laplace and ALE Navier-Stokes equations set is not consistent.",err,error,*999)
          DO componentIdx=1,fluidNumberOfDimensions
            CALL FieldVariable_ParametersToFieldVariableParametersCopy(laplaceDependentVariable,FIELD_VALUES_SET_TYPE, &
              & componentIdx,fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,componentIdx,err,error,*999)
          ENDDO !componentIdx
          !Find the solid mechanics equations set
          equationsSetIdx=1
          solidEquationsSetFound=.FALSE.
          DO WHILE (.NOT.solidEquationsSetFound.AND.equationsSetIdx<=fsiSolverMapping%numberOfEquationsSets)
            NULLIFY(solidEquationsSet)
            CALL SolverMapping_EquationsSetGet(fsiSolverMapping,equationsSetIdx,solidEquationsSet,err,error,*999)
            IF(solidEquationsSet%specification(1)==EQUATIONS_SET_ELASTICITY_CLASS.AND. &
              & solidEquationsSet%specification(2)==EQUATIONS_SET_FINITE_ELASTICITY_TYPE.AND. &
              & (solidEquationsSet%specification(3)==EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE.OR. &
              & solidEquationsSet%specification(3)==EQUATIONS_SET_MR_AND_GROWTH_LAW_IN_CELLML_SUBTYPE.OR. &
              & solidEquationsSet%specification(3)==EQUATIONS_SET_COMPRESSIBLE_FINITE_ELASTICITY_SUBTYPE.OR. &
              & solidEquationsSet%specification(3)==EQUATIONS_SET_DYNAMIC_ST_VENANT_KIRCHOFF_SUBTYPE.OR. &
              & solidEquationsSet%specification(3)==EQUATIONS_SET_DYNAMIC_MOONEY_RIVLIN_SUBTYPE)) THEN
              solidEquationsSetFound=.TRUE.
            ELSE
              equationsSetIdx=equationsSetIdx+1
            ENDIF
          ENDDO
          IF(.NOT.solidEquationsSetFound) &
            & CALL FlagError("Solid equations set not found when trying to update ALE mesh.",err,error,*999)
          NULLIFY(solidDependentField)
          CALL EquationsSet_DependentFieldGet(solidEquationsSet,solidDependentField,err,error,*999)
          NULLIFY(solidDependentVariable)
          CALL Field_VariableGet(solidDependentField,FIELD_U_VARIABLE_TYPE,solidDependentVariable,err,error,*999)
          !Loop over the interface nodes and see if there has been any change in the solid position since the
          !beginning of the time loop. Add this change to the mesh displacements.
          NULLIFY(fsiInterfaceCondition)
          CALL SolverMapping_InterfaceConditionGet(fsiSolverMapping,1,fsiInterfaceCondition,err,error,*999)
          NULLIFY(fsiInterface)
          CALL InterfaceCondition_InterfaceGet(fsiInterfaceCondition,fsiInterface,err,error,*999)
          NULLIFY(meshConnectivity)
          CALL Interface_MeshConnectivityGet(fsiInterface,meshConnectivity,err,error,*999)
          NULLIFY(interfaceGeometricField)
          CALL InterfaceCondition_GeometricFieldGet(fsiInterfaceCondition,interfaceGeometricField,err,error,*999)
          CALL Field_NumberOfComponentsGet(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,interfaceNumberOfDimensions, &
            & err,error,*999)
          NULLIFY(interfaceGeometricVariable)
          CALL Field_VariableGet(interfaceGeometricField,FIELD_U_VARIABLE_TYPE,interfaceGeometricVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(interfaceGeometricVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            CALL FieldVariable_ComponentInterpolationGet(interfaceGeometricVariable,componentIdx,interpolationType, &
              & err,error,*999)
            SELECT CASE(interpolationType)
            CASE(FIELD_NODE_BASED_INTERPOLATION)
              NULLIFY(domain)
              CALL FieldVariable_DomainGet(interfaceGeometricVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              CALL DomainNodes_TotalNumberOfNodesGet(domainNodes,totalNumberOfNodes,err,error,*999)
              DO nodeIdx=1,totalNumberOfNodes
                !!TODO: the fluid/solid meshes should not be position dependent in the coupled nodes
                CALL InterfaceMeshConnectivity_CoupledNodeNumberGet(meshConnectivity,1,nodeIdx,solidNode,err,error,*999)
                CALL InterfaceMeshConnectivity_CoupledNodeNumberGet(meshConnectivity,2,nodeIdx,fluidNode,err,error,*999)
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
                  DO versionIdx=1,numberOfVersions
                    CALL FieldVariable_ParameterSetGetNode(solidDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
                      & derivativeIdx,solidNode,componentIdx,solidNodePosition,err,error,*999)
                    CALL FieldVariable_ParameterSetGetNode(solidDependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,versionIdx, &
                      & derivativeIdx,solidNode,componentIdx,previousSolidNodePosition,err,error,*999)
                    solidDelta=solidNodePosition-previousSolidNodePosition
                    IF(ABS(solidDelta)>ZERO_TOLERANCE) THEN
                      CALL FieldVariable_ParameterSetAddNode(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
                        & versionIdx,derivativeIdx,fluidNode,componentIdx,solidDelta,err,error,*999)
                    ENDIF
                  ENDDO !versionIdx
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            CASE DEFAULT
              CALL FlagError("Interface geometric component does not have node based interpolation.",err,error,*999)
            END SELECT
          ENDDO !componentIdx
          NULLIFY(meshDisplacementValues)
          CALL FieldVariable_ParameterSetDataGet(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & meshDisplacementValues,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(fluidGeometricVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            NULLIFY(domain)
            CALL FieldVariable_DomainGet(fluidGeometricVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
                DO versionIdx=1,numberOfVersions
                  CALL FieldVariable_LocalNodeDOFGet(fluidGeometricVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                    & localDOF,err,error,*999)
                  CALL FieldVariable_ParameterSetAddLocalDOF(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                    & meshDisplacementValues(localDOF),err,error,*999)
                ENDDO !versionIdx
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
          CALL FieldVariable_ParameterSetDataRestore(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & meshDisplacementValues,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(fluidGeometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Now use displacement values to calculate velocity values
          alpha=1.0_DP/timeIncrement
          CALL FieldVariable_ParameterSetsCopy(fluidIndependentVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
        CASE DEFAULT
          localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is not valid for a Finite Elasticity-Navier-Stokes type of a multi physics problem class."
        END SELECT
      CASE DEFAULT
        localError="Problem type "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for a multi physics problem class."
      END SELECT
    CASE DEFAULT
      localError="Problem class "//TRIM(NumberToVString(pSpecification(1),"*",err,error))//" is not valid."
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
  SUBROUTINE NavierStokes_PreSolveALEUpdateParameters(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,localDOF,nodeIdx,numberOfComponents,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfVariables,pSpecification(3),solveType,variableIdx,variableType
    REAL(DP) :: currentTime,timeIncrement
    REAL(DP), POINTER :: meshStiffValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: dependentField,independentField
    TYPE(FieldVariableType), POINTER :: dependentVariable,independentVariable
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("NavierStokes_PreSolveALEUpdateParameters",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)

    NULLIFY(solvers)
    CALL Solver_SolversGet(solver,solvers,err,error,*999)
    NULLIFY(controlLoop)
    CALL Solvers_ControlLoopGet(solvers,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(1))
    CASE(PROBLEM_FLUID_MECHANICS_CLASS)
      SELECT CASE(pSpecification(3))
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
        CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
        IF(solveType/=SOLVER_LINEAR_TYPE) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
        !Get the independent field for the ALE Navier-Stokes problem
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(independentField)
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
        NULLIFY(independentVariable)
        CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
        NULLIFY(meshStiffValues)
        CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_VALUES_SET_TYPE,meshStiffValues,err,error,*999)
        CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
        DO variableIdx=1,numberOfVariables
          NULLIFY(dependentVariable)
          CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            NULLIFY(domain)
            CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
            DO nodeIdx=1,numberOfNodes
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                  & err,error,*999)
                !Calculation of K values dependent on current mesh topology
                meshStiffValues(localDOF)=1.0_DP
                CALL FieldVariable_ParameterSetUpdateLocalDOF(independentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                  & meshStiffValues(localDOF),err,error,*999)
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
        ENDDO !variableIdx
        CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_VALUES_SET_TYPE,meshStiffValues,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
          & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_MULTI_PHYSICS_CLASS)
      SELECT CASE(pSpecification(2))
      CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_TYPE)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_GROWTH_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_NAVIER_STOKES_ALE_SUBTYPE, &
          & PROBLEM_DYNAMIC_FINITE_ELASTICITY_RBS_NAVIER_STOKES_ALE_SUBTYPE)
          CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
          IF(solveType/=SOLVER_LINEAR_TYPE) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
          !Get the independent field for the ALE Navier-Stokes problem
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          NULLIFY(solverMapping)
          CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(independentField)
          CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
          NULLIFY(independentVariable)
          CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
          NULLIFY(meshStiffValues)
          CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_VALUES_SET_TYPE,meshStiffValues,err,error,*999)
          CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
            DO componentIdx=1,numberOfComponents
              NULLIFY(domain)
              CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
              DO nodeIdx=1,numberOfNodes
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                    & err,error,*999)
                  !Calculation of K values dependent on current mesh topology
                  meshStiffValues(localDOF)=1.0_DP
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(independentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                    & meshStiffValues(localDOF),err,error,*999)
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
          ENDDO !variableIdx
          CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_VALUES_SET_TYPE,meshStiffValues,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is not valid for a FiniteElasticity-NavierStokes type of a multi physics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second problem specification of "//TRIM(NumberToVString(pSpecification(2),"*",err,error))// &
          & " is not valid for NavierStokes_PreSolveALEUpdateParameters of a multi physics problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The first problem specification of "//TRIM(NumberToVString(pSpecification(1),"*",err,error))// &
        & " is not valid for NavierStokes_PreSolveALEUpdateParameters."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_PreSolveALEUpdateParameters")
    RETURN
999 ERRORSEXITS("NavierStokes_PreSolveALEUpdateParameters",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PreSolveALEUpdateParameters

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE NavierStokes_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,currentLoopIteration,equationsSetIdx,esSpecification(3),fileNameLength, &
      & inputIteration,numberOfDimensions,numberOfEquationsSets,outputIterationNumber,outputType,pSpecification(3), &
      & solverGlobalNumber
    REAL(DP) :: currentTime,timeIncrement,startTime,stopTime
    CHARACTER(20) :: file,outputFile                                    
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldsType), POINTER :: fields
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: region
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError,method,vFileName,filename
   
    ENTERS("NavierStokes_PostSolveOutputData",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL SYSTEM('mkdir -p ./output')
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE,PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        filename="./output/"//"STATIC_SOLUTION"
        method="FORTRAN"
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
        ENDIF
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        NULLIFY(fields)
        CALL Region_FieldsGet(region,fields,err,error,*999)
        CALL FIELD_IO_NODES_EXPORT(Fields,filename,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(Fields,filename,method,err,error,*999)
      ENDDO !equationsSetIdx

    CASE(PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_ALE_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==2) THEN
        CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime, &
          & currentLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Make sure the equations sets are up to date
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          IF(outputIterationNumber/=0) THEN
            IF(currentTime<=stopTime) THEN
              WRITE(outputFile,'("TimeStep_",I0)') currentLoopIteration
              file=outputFile
              method="FORTRAN"
              IF(MOD(currentLoopIteration,outputIterationNumber)==0)  THEN
                !Use standard field IO routines (also only export nodes after first step as not a moving mesh case)
                fileNameLength = LEN_TRIM(outputFile)
                vFileName = outputFile(1:FileNameLength)
                IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) &
                  & CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                NULLIFY(region)
                CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
                NULLIFY(fields)
                CALL Region_FieldsGet(region,fields,err,error,*999)
                CALL FIELD_IO_NODES_EXPORT(Fields,VFileName,method,err,error,*999)
                IF(currentLoopIteration==0) THEN
                  IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) &
                    & CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export elements... ",err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,method,err,error,*999)
                ENDIF
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                ENDIF
              ENDIF
              NULLIFY(equationsAnalytic)
              CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
              IF(ASSOCIATED(equationsAnalytic)) THEN
                CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
                IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                  & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                  NULLIFY(dependentField)
                  CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                  CALL AnalyticAnalysis_Output(dependentField,file,err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx
      ENDIF

    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE,PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
        & currentLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
          IF(outputIterationNumber/=0) THEN
            IF(currentTime<=stopTime) THEN
              WRITE(outputFile,'("TimeStep3D_",I0)') currentLoopIteration
              file=outputFile
              method="FORTRAN"
              IF(MOD(currentLoopIteration,outputIterationNumber)==0)  THEN
                !Use standard field IO routines (also only export nodes after first step as not a moving mesh case)
                fileNameLength = LEN_TRIM(outputFile)
                vFileName = outputFile(1:FileNameLength)
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                NULLIFY(region)
                CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
                NULLIFY(fields)
                CALL Region_FieldsGet(region,fields,err,error,*999)
                CALL FIELD_IO_NODES_EXPORT(fields,vFileName,method,err,error,*999)
                IF(currentLoopIteration==0) THEN
                  IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) &
                    & CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export elements... ",err,error,*999)
                  CALL FIELD_IO_ELEMENTS_EXPORT(Fields,VFileName,method,err,error,*999)
                ENDIF
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDIF

        CASE(EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
          & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)
          IF(outputIterationNumber/=0) THEN
            IF(currentTime<=stopTime) THEN
              NULLIFY(region)
              CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
              NULLIFY(fields)
              CALL Region_FieldsGet(region,fields,err,error,*999)
              file=outputFile
              filename="TimeStep1D_"//TRIM(NumberToVString(currentLoopIteration,"*",err,error))
              method="FORTRAN"
              IF(MOD(currentLoopIteration,outputIterationNumber)==0)  THEN
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                ENDIF
                CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
                ! Only export elements on first iteration (non-moving mesh case)
                IF(currentLoopIteration==0) THEN
                  CALL FIELD_IO_ELEMENTS_EXPORT(region%FIELDS,filename,method,err,error,*999)
                ENDIF
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                ENDIF
!!TODO: What are these statements doing?
                NULLIFY(geometricField)
                CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
                CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is not valid for a dynamic solver output for a multiscale Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !equationsSetIdx

    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)

      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
        & currentLoopIteration,outputIterationNumber,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(outputIterationNumber/=0) THEN
          IF(currentTime<=stopTime) THEN
            IF(currentLoopIteration<10) THEN
              WRITE(outputFile,'("TIME_STEP_000",I0)') currentLoopIteration
            ELSE IF(currentLoopIteration<100) THEN
              WRITE(outputFile,'("TIME_STEP_00",I0)') currentLoopIteration
            ELSE IF(currentLoopIteration<1000) THEN
              WRITE(outputFile,'("TIME_STEP_0",I0)') currentLoopIteration
            ELSE IF(currentLoopIteration<10000) THEN
              WRITE(outputFile,'("TIME_STEP_",I0)') currentLoopIteration
            ENDIF
            NULLIFY(region)
            CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
            file=outputFile
            filename="./output/"//"MainTime_"//TRIM(NumberToVString(currentLoopIteration,"*",err,error))
            method="FORTRAN"
            IF(MOD(currentLoopIteration,outputIterationNumber)==0)  THEN
              IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
              ENDIF
              NULLIFY(fields)
              CALL Region_FieldsGet(region,fields,err,error,*999)
              CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
              CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
              IF(outputType >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
              ENDIF
!!TODO: What are the next three lines doing?
              NULLIFY(geometricField)
              CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
              CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            ENDIF
            NULLIFY(equationsAnalytic)
            CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
            IF(ASSOCIATED(equationsAnalytic)) THEN
              CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
              IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1.OR. &
                & analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1) THEN
                NULLIFY(dependentField)
                CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                CALL AnalyticAnalysis_Output(dependentField,file,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("NavierStokes_PostSolveOutputData",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Sets up analytic parameters and calls NavierStokes_AnalyticFunctionsEvaluate to evaluate solutions to analytic problems
  SUBROUTINE NavierStokes_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionsCheckVariable,boundaryCount,componentIdx,derivativeIdx,dimensionIdx, &
      & elementIdx,elementNodeIdx,elementNodeNumber,elementNumber,globalDerivativeIndex,globalDOF,I,interpolationType,J,K, &
      & localDOF,maximumNumberOfElementParameters,nodeIdx,nodeNumber,numberOfComponents,numberOfDimensions,numberOfNodes, &
      & numberOfNodeDerivatives,numberOfNodesXiCoord(3),numberOfParameters,numberOfVariables,numberOfVersions,numberOfXi, &
      & parameterIdx,userNodeNumber,variableIdx,variableType,versionIdx
    REAL(DP) :: initialValue,nodeAnalyticParameters(10),time,tCoordinates(20,3),VALUE,x(3),xiCoordinates(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode
    TYPE(BasisType), POINTER :: basis
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldInterpolatedPointType), POINTER :: interpolatedPoint
    TYPE(FieldInterpolationParametersType), POINTER :: interpolationParameters
    TYPE(FieldParameterSetType), POINTER :: analyticParameterSet,pressureParameterSet
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,analyticVariable,materialsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_BoundaryConditionsAnalyticCalculate",err,error,*999)

    boundaryCount=0
    xiCoordinates(3)=0.0_DP

    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    ! Geometric parameters
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    NULLIFY(geometricParameters)
    CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(interpolationParameters)
    CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,interpolationParameters,err,error,*999)
    NULLIFY(interpolatedPoint)
    CALL Field_InterpolatedPointInitialise(interpolationParameters,interpolatedPoint,err,error,*999)
    !Analytic parameters
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    NULLIFY(analyticField)
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    NULLIFY(analyticVariable)
    NULLIFY(analyticParameters)
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    ENDIF
    !Materials parameters
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    NULLIFY(materialsVariable)
    NULLIFY(materialsParameters)
    IF(ASSOCIATED(materialsField)) THEN
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    ENDIF
    CALL EquationsSet_AnalyticTimeGet(equationsSet,time,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        boundaryCount=0
        CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,interpolationType,err,error,*999)
        IF(interpolationType==FIELD_NODE_BASED_INTERPOLATION) &
          & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
          CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          CALL DomainNodes_NodeSurroundingElementGet(domainNodes,1,nodeIdx,elementIdx,err,error,*999)
          NULLIFY(basis)
          CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,interpolationParameters,err,error,*999)
          elementNodeIdx=0
          xiCoordinates=0.0_DP
          CALL Basis_NumberOfXiGet(basis,numberOfXi,err,error,*999)
          numberOfNodesXiCoord=1
          CALL Basis_NumberOfNodesXiCGet(basis,numberOfNodesXiCoord,err,error,*999)

          ! --- Calculate analytic profile for validation ---
          SELECT CASE(analyticFunctionType)
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
            IF(variableIdx < 3) THEN
              !Get geometric position info for this node
              DO dimensionIdx=1,numberOfDimensions
                CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOF,err,error,*999)
                x(dimensionIdx)=geometricParameters(localDOF)
              ENDDO !dimensionIdx
              DO derivativeIdx=1,numberOfNodeDerivatives
                CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
                CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
                DO versionIdx=1,numberOfVersions
                  CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType,globalDerivativeIndex, &
                    & componentIdx,numberOfDimensions,numberOfComponents,analyticParameters,materialsParameters,VALUE, &
                    & err,error,*999)
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                    & localDOF,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOF, &
                    & VALUE,err,error,*999)
                ENDDO !versionIdx
              ENDDO !derivativeIdx
            ENDIF ! variableIdx < 3

            ! --- Set velocity boundary conditions with analytic value ---
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN, &
            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
            ! Get geometric position info for this node
            DO dimensionIdx=1,numberOfDimensions
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOF,err,error,*999)
              x(dimensionIdx)=geometricParameters(localDOF)
            ENDDO !dimensionIdx
            !Loop over the derivatives
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              IF(componentIdx<=numberOfXi.OR.analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
                DO versionIdx=1,numberOfVersions
                  !Get global and local DOF indices
                  CALL FieldVariable_ComponentDOFGetUserNode(dependentVariable,versionIdx,derivativeIdx,userNodeNumber, &
                    & componentIdx,localDof,globalDof,err,error,*999)
                  IF(analyticFunctionType==EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID) THEN
                    CALL FieldVariable_NumberOfComponentsGet(analyticVariable,numberOfParameters,err,error,*999)
                    DO parameterIdx=1,numberOfParameters
                      !populate nodeAnalyticParameters
                      CALL FieldVariable_ParameterSetGetLocalNode(analyticVariable,FIELD_VALUES_SET_TYPE, &
                        & versionIdx,derivativeIdx,nodeIdx,parameterIdx,nodeAnalyticParameters(parameterIdx), &
                        & err,error,*999)
                    ENDDO !parameterIdx
                    CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType,globalDerivativeIndex, &
                      & componentIdx,numberOfDimensions,numberOfComponents,nodeAnalyticParameters,materialsParameters,VALUE, &
                      & err,error,*999)
                  ELSE
                    CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType,globalDerivativeIndex, &
                      & componentIdx,numberOfDimensions,numberOfComponents,analyticParameters,materialsParameters,VALUE, &
                      & err,error,*999)
                  ENDIF
                  !update analytic field values
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOF, &
                    & VALUE,err,error,*999)
                  IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                    CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
                    IF(boundaryNode) THEN
                      NULLIFY(boundaryConditionsVariable)
                      CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                        & err,error,*999)
                      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                          & boundaryConditionsCheckVariable,err,error,*999)
                        !update dependent field values if fixed inlet or pressure BC
                        IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET .OR. &
                          & boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_PRESSURE) THEN
                          !Set DOF values
                          CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDof,VALUE, &
                            & err,error,*999)
                        ELSE IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_PRESSURE) THEN
                          !!Set neumann boundary pressure value on pressure nodes
                          !CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,1,1, &
                          !  & nodeIdx,componentIdx,VALUE,err,error,*999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !versionIdx
              ENDIF
            ENDDO !derivativeIdx
            
            ! --- Set Flow rate boundary conditions with analytic value ---
          CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA, &
            & EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
            !Get geometric position info for this node
            DO dimensionIdx=1,numberOfDimensions
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOF,err,error,*999)
              x(dimensionIdx)=geometricParameters(localDOF)
            ENDDO !dimensionIdx
            !Loop over the derivatives
            CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              IF(componentIdx==1.AND.variableType==FIELD_U_VARIABLE_TYPE) THEN
                CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
                DO versionIdx=1,numberOfVersions
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                    & err,error,*999)
                  IF(boundaryNode) THEN
                    NULLIFY(boundaryConditionsVariable)
                    CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                      & err,error,*999)
                    IF(ASSOCIATED(boundaryConditionsVariable)) THEN
                      CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,localDOF, &
                        & boundaryConditionsCheckVariable,err,error,*999)
                      IF(boundaryConditionsCheckVariable==BOUNDARY_CONDITION_FIXED_INLET) THEN
                        CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType, &
                          & globalDerivativeIndex,componentIdx,numberOfXi,numberOfComponents,analyticParameters, &
                          & materialsParameters,VALUE,err,error,*999)
                        !If we are a boundary node then set the analytic value on the boundary
                        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF,VALUE, &
                          & err,error,*999)
                      ELSE
                        CALL FieldVariable_ParameterSetGetLocalNode(dependentVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
                          & derivativeIdx,nodeIdx,componentIdx,VALUE,err,error,*999)
                        CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOF, &
                          & VALUE,err,error,*999)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !versionIdx
              ENDIF
            ENDDO !derivativeIdx
           
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
            CALL DomainElements_MaxElementParametersGet(domainElements,maximumNumberOfElementParameters,err,error,*999)
            IF(maximumNumberOfElementParameters==4.AND.numberOfDimensions==2.OR. &
              & maximumNumberOfElementParameters==9.OR. &
              & maximumNumberOfElementParameters==16.OR. &
              & maximumNumberOfElementParameters==8.OR. &
              & maximumNumberOfElementParameters==27.OR. &
              & maximumNumberOfElementParameters==64) THEN
              DO K=1,numberOfNodesXiCoord(3)
                DO J=1,numberOfNodesXiCoord(2)
                  DO I=1,numberOfNodesXiCoord(1)
                    elementNodeIdx=elementNodeIdx+1
                    CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,elementNodeNumber,err,error,*999)
                    IF(elementNodeNumber==nodeIdx) EXIT
                    xiCoordinates(1)=xiCoordinates(1)+(1.0_DP/(numberOfNodesXiCoord(1)-1))
                  ENDDO !I
                  CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,elementNodeNumber,err,error,*999)
                  IF(elementNodeNumber==nodeIdx) EXIT
                  xiCoordinates(1)=0.0_DP
                  xiCoordinates(2)=xiCoordinates(2)+(1.0_DP/(numberOfNodesXiCoord(2)-1))
                ENDDO !J
                CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,elementNodeNumber,err,error,*999)
                IF(elementNodeNumber==nodeIdx) EXIT
                xiCoordinates(1)=0.0_DP
                xiCoordinates(2)=0.0_DP
                IF(numberOfNodesXiCoord(3)/=1) xiCoordinates(3)=xiCoordinates(3)+(1.0_DP/(numberOfNodesXiCoord(3)-1))
              ENDDO !K
              CALL Field_InterpolateXi(NO_PART_DERIV,xiCoordinates,interpolatedPoint,err,error,*999)
              !Tri/Tet
              !\todo: Use boundary flag
            ELSE IF(maximumNumberOfElementParameters==3) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==6) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
              tCoordinates(4,1:2)=[0.5_DP,0.5_DP]
              tCoordinates(5,1:2)=[1.0_DP,0.5_DP]
              tCoordinates(6,1:2)=[0.5_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==2) THEN
              tCoordinates(1,1:2)=[0.0_DP,1.0_DP]
              tCoordinates(2,1:2)=[1.0_DP,0.0_DP]
              tCoordinates(3,1:2)=[1.0_DP,1.0_DP]
              tCoordinates(4,1:2)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(5,1:2)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP]
              tCoordinates(6,1:2)=[1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(7,1:2)=[1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(8,1:2)=[2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(9,1:2)=[1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(10,1:2)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP]
            ELSE IF(maximumNumberOfElementParameters==4) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==10.AND.numberOfDimensions==3) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
              tCoordinates(5,1:3)=[0.5_DP,0.5_DP,1.0_DP]
              tCoordinates(6,1:3)=[0.5_DP,1.0_DP,0.5_DP]
              tCoordinates(7,1:3)=[0.5_DP,1.0_DP,1.0_DP]
              tCoordinates(8,1:3)=[1.0_DP,0.5_DP,0.5_DP]
              tCoordinates(9,1:3)=[1.0_DP,1.0_DP,0.5_DP]
              tCoordinates(10,1:3)=[1.0_DP,0.5_DP,1.0_DP]
            ELSE IF(maximumNumberOfElementParameters==20) THEN
              tCoordinates(1,1:3)=[0.0_DP,1.0_DP,1.0_DP]
              tCoordinates(2,1:3)=[1.0_DP,0.0_DP,1.0_DP]
              tCoordinates(3,1:3)=[1.0_DP,1.0_DP,0.0_DP]
              tCoordinates(4,1:3)=[1.0_DP,1.0_DP,1.0_DP]
              tCoordinates(5,1:3)=[1.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(6,1:3)=[2.0_DP/3.0_DP,1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(7,1:3)=[1.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(8,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(9,1:3)=[1.0_DP/3.0_DP,1.0_DP,1.0_DP]
              tCoordinates(10,1:3)=[2.0_DP/3.0_DP,1.0_DP,1.0_DP]
              tCoordinates(11,1:3)=[1.0_DP,1.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(12,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP/3.0_DP]
              tCoordinates(13,1:3)=[1.0_DP,1.0_DP,1.0_DP/3.0_DP]
              tCoordinates(14,1:3)=[1.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(15,1:3)=[1.0_DP,1.0_DP/3.0_DP,1.0_DP]
              tCoordinates(16,1:3)=[1.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(17,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
              tCoordinates(18,1:3)=[2.0_DP/3.0_DP,2.0_DP/3.0_DP,1.0_DP]
              tCoordinates(19,1:3)=[2.0_DP/3.0_DP,1.0_DP,2.0_DP/3.0_DP]
              tCoordinates(20,1:3)=[1.0_DP,2.0_DP/3.0_DP,2.0_DP/3.0_DP]
            ENDIF
            DO K=1,maximumNumberOfElementParameters
              CALL DomainElements_ElementNodeGet(domainElements,k,elementIdx,elementNodeNumber,err,error,*999)
              IF(elementNumber==nodeIdx) EXIT
            ENDDO !K
            CALL Field_InterpolateXi(NO_PART_DERIV,tCoordinates(k,1:numberOfDimensions),interpolatedPoint,err,error,*999)
            x=0.0_DP
            x(1:numberOfDimensions)=interpolatedPoint%values(1:numberOfDimensions,NO_PART_DERIV)
           CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            !Loop over the derivatives
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              CALL NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType,globalDerivativeIndex, &
                & componentIdx,numberOfDimensions,numberOfComponents,analyticParameters,materialsParameters,VALUE,err,error,*999)
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx,localDOF, &
                  & err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOF,VALUE, &
                  & err,error,*999)
                IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                  IF(boundaryNode) THEN
                    !If we are a boundary node then set the analytic value on the boundary
                    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOF,BOUNDARY_CONDITION_FIXED, &
                      & VALUE,err,error,*999)
                    ! \todo: This is just a workaround for linear pressure fields in simplex element components
                    IF(componentIdx>numberOfDimensions) THEN
                      IF(maximumNumberOfElementParameters==3) THEN
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
                            CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOF, &
                              & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                            boundaryCount=boundaryCount+1
                          ENDIF
                        ENDIF
                      ELSE IF(maximumNumberOfElementParameters==4.AND.numberOfDimensions==3) THEN
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
                            CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOF, &
                              & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                            boundaryCount=boundaryCount+1
                          ENDIF
                        ENDIF
                        ! \todo: This is how it should be if adjacent elements would be working
                      ELSE IF(boundaryCount==0) THEN
                        CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOF, &
                          & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                        boundaryCount=boundaryCount+1
                      ENDIF
                    ENDIF
                  ELSE
                    !Set the initial condition.
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                      & initialValue,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !versionIdx
            ENDDO !derivativeIdx
            
          CASE DEFAULT
            localError="Analytic Function Type "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
              & " is not yet implemented for a Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        ENDDO !nodeIdx
      ENDDO !componentIdx
      !Update ghost and boundary values from local
      NULLIFY(pressureParameterSet)
      CALL FieldVariable_ParameterSetCheck(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,pressureParameterSet, &
        & err,error,*999)
      IF(ASSOCIATED(pressureParameterSet)) THEN
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      ENDIF
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    END DO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    CALL Field_InterpolatedPointFinalise(interpolatedPoint,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(interpolationParameters,err,error,*999)

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
  SUBROUTINE NavierStokes_AnalyticFunctionsEvaluate(analyticFunctionType,x,time,variableType,globalDerivIndex, &
    & componentNumber,numberOfDimensions,numberOfComponents,analyticParameters,materialsParameters,VALUE,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: x(:) !<x(dimension_idx). The geometric position to evaluate at (includes Y,Z for higher dim problems)
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivIndex !<The global derivative to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    INTEGER(INTG), INTENT(IN) :: numberOfDimensions !<The number of geometric dimensions
    INTEGER(INTG), INTENT(IN) :: numberOfComponents !<The number of components for the dependent field
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: i,j,n,m
    REAL(DP) :: lParam,hParam,uParam,pParam,muParam,nuParam,rhoParam,internalTime,currentTime,kParam
    REAL(DP) :: amplitude,yOffset,period,phaseShift,frequency,s,startTime,stopTime,tt,tmax,Qo
    REAL(DP) :: componentCoeff(4),delta(300),t(300),q(300)
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_AnalyticFunctionsEvaluate",err,error,*999)

    !\todo: Introduce user-defined or default values instead for density and viscosity
    internalTime=time
    currentTime=time

    SELECT CASE(analyticFunctionType)

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_POISEUILLE)
      !For fully developed 2D laminar flow through a channel, NSE should yield a parabolic profile,
      !U = Umax(1-y^2/H^2), Umax = (-dP/dx)*(H^2/(2*MU)), Umax = (3/2)*Umean
      !Note: assumes a flat inlet profile (uParam = Umean).
      !Nonlinear terms from NSE will effectively be 0 for Poiseuille flow
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          lParam = analyticParameters(1) ! channel length in x-direction
          hParam = analyticParameters(2) ! channel height in y-direction
          uParam = analyticParameters(3) ! mean (inlet) velocity
          pParam = analyticParameters(4) ! pressure value at outlet
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=(3.0_DP/2.0_DP)*uParam*(1.0_DP-((x(2)-hParam)**2)/(hParam**2))
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=0.0_DP
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE = (3.0_DP*muParam*uParam*(x(1)-lParam))/(hParam**2)+pParam
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE( NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_TAYLOR_GREEN)
      !Exact solution to 2D laminar, dynamic, nonlinear Taylor-Green vortex decay
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        nuParam = muParam/rhoParam ! kinematic viscosity
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          uParam = analyticParameters(1) ! characteristic velocity (initial amplitude)
          lParam = analyticParameters(2) ! length scale for square
          kParam = 2.0_DP*PI/lParam   ! scale factor for equations
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=-1.0_DP*uParam*COS(kParam*x(1))*SIN(kParam*x(2))*EXP(-2.0_DP*(kParam**2)*nuParam*currentTime)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=uParam*SIN(kParam*x(1))*COS(kParam*x(2))*EXP(-2.0_DP*(kParam**2)*nuParam*currentTime)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE =-1.0_DP*(uParam**2)*(rhoParam/4.0_DP)*(COS(2.0_DP*kParam*x(1))+ &
                & COS(2.0_DP*kParam*x(2)))*(EXP(-4.0_DP*(kParam**2)*nuParam*currentTime))
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE( NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_AORTA)
      SELECT CASE(numberOfDimensions)
      CASE(1)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !Input function
              period = 800
              tt=MOD(TIME,period)
              tmax=150.0_DP
              Qo=100000.0_DP
              VALUE=(Qo*tt/(tmax**2.0_DP))*EXP(-(tt**2.0_DP)/(2.0_DP*(tmax**2.0_DP)))
            ELSE
              CALL FlagError("Incorrect component specification for Aorta flow rate waveform ",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))// " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
          ! Do nothing
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Aorta flowrate waveform for "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
          & " dimension problem has not yet been implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_FLOWRATE_OLUFSEN)
      SELECT CASE(numberOfDimensions)
      CASE(1)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
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
                IF(t(j) <= (time/period)) THEN
                  m=j
                ENDIF
              END DO
              !Evaluate interpolant
              s=(time/period)-t(m)
              VALUE=(q(m)+s*delta(m))
            ELSE
              CALL FlagError("Incorrect component specification for Olufsen flow rate waveform.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE)
          ! Do nothing
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Olufsen flowrate waveform for "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
          & " dimension problem has not yet been implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_SINUSOID)
      !Returns a sinusoidal value for boundary nodes
      SELECT CASE(numberOfDimensions)
      CASE(2,3)
        componentCoeff(1) = analyticParameters(1)
        componentCoeff(2) = analyticParameters(2)
        componentCoeff(3) = analyticParameters(3)
        componentCoeff(4) = analyticParameters(4)
        amplitude = analyticParameters(5)
        yOffset = analyticParameters(6)
        frequency = analyticParameters(7)
        phaseShift = analyticParameters(8)
        startTime = analyticParameters(9)
        stopTime = analyticParameters(10)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(currentTime > startTime - ZERO_TOLERANCE .AND. &
              &  currentTime < stopTime + ZERO_TOLERANCE) THEN
              VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*currentTime+phaseShift))
            ELSE
              VALUE= componentCoeff(componentNumber)*(yOffset + amplitude*SIN(frequency*stopTime+phaseShift))
            ENDIF
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Sinusoidal analytic types for "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
          & " dimensional problems have not yet been implemented."
        CALL FlagError(localError,err,error,*999)
      END SELECT

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_ONE_DIM_1)
      IF(numberOfDimensions==1.AND.numberOfComponents==3) THEN
        !Polynomial function
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate Q
              VALUE=x(1)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate A
              VALUE=x(1)**2/10.0_DP**2
            ELSE IF(componentNumber==3) THEN
              !calculate P
              VALUE=x(1)**2/10.0_DP**2
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_1)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Polynomial function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=x(2)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=x(1)**2/10.0_DP**2
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=2.0_DP/3.0_DP*x(1)*(3.0_DP*muParam*10.0_DP**2-rhoParam*x(1)**2*x(2))/(10.0_DP ** 4)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE= 0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_2)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Exponential function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE= EXP((x(1)-x(2))/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE= EXP((x(1)-x(2))/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE= 2.0_DP*muParam/10.0_DP*EXP((x(1)-x(2))/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
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
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_3)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Sine and cosine function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(2.0_DP*PI*x(1)/10.0_DP)*SIN(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=COS(2.0_DP*PI*x(1)/10.0_DP)*COS(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=4.0_DP*muParam*PI/10.0_DP*SIN(2.0_DP*PI*x(2)/10.0_DP)*COS(2.0_DP*PI*x(1)/10.0_DP)+ &
                & 0.5_DP*rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)*COS(2.0_DP*PI*x(1)/10.0_DP)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=16.0_DP*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*x(2)/ 10.0_DP)*cos(2.0_DP*PI*x(1)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_TWO_DIM_5)
      IF(numberOfDimensions==2.AND.numberOfComponents==3) THEN
        !Taylor-Green vortex solution
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(x(1)/10.0_DP*2.0_DP*PI)*COS(x(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
              VALUE=SIN(x(1)/10.0_DP*PI)*COS(x(2)/10.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
              !VALUE=SIN(x(1))*COS(x(2))
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-COS(x(1)/10.0_DP*2.0_DP*PI)*SIN(x(2)/10.0_DP*2.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
              VALUE=-COS(x(1)/10.0_DP*PI)*SIN(x(2)/10.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
              !VALUE=-COS(x(1))*SIN(x(2))
            ELSE IF(componentNumber==3) THEN
              !calculate p
              VALUE=rhoParam/4.0_DP*(COS(2.0_DP*x(1)/10.0_DP*2.0_DP*PI)+COS(2.0_DP*x(2)/10.0_DP*2.0_DP*PI))* &
                & EXP(-4.0_DP*muParam/rhoParam*currentTime)
              VALUE=rhoParam/4.0_DP*(COS(2.0_DP*x(1)/10.0_DP*PI)+COS(2.0_DP*x(2)/10.0_DP*PI))* &
                & EXP(-4.0_DP*muParam/rhoParam*currentTime)
              !VALUE=rhoParam/4.0_DP*(COS(2.0_DP*x(1))+COS(2.0_DP*x(2)))
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
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
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_1)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Polynomial function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=x(2)**2/10.0_DP**2+x(3)**2/10.0_DP**2
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=x(1)**2/10.0_DP**2+x(3)**2/10.0_DP** 2
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=x(1)**2/10.0_DP**2+x(2)**2/10.0_DP** 2
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=2.0_DP/3.0_DP*x(1)*(6.0_DP*muParam*10.0_DP**2-rhoParam*x(2)*x(1)**2-3.0_DP* &
                & rhoParam*x(2)* &
                & x(3)**2-rhoParam*x(3)*x(1)**2-3.0_DP*rhoParam*x(3)*x(2)**2)/(10.0_DP**4)
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_2)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Exponential function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=EXP((x(1)-x(2))/10.0_DP)+EXP((x(3)-x(1))/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=EXP((x(1)-x(2))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=EXP((x(3)-x(1))/10.0_DP)+EXP((x(2)-x(3))/10.0_DP)
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=1.0_DP/10.0_DP*(2.0_DP*muParam*EXP((x(1)-x(2))/10.0_DP)- &
                & 2.0_DP*muParam*EXP((x(3)-x(1))/10.0_DP)+rhoParam*10.0_DP*EXP((x(1)-x(3))/10.0_DP)+ &
                & rhoParam*10.0_DP*EXP((x(2)-x(1))/10.0_DP))
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-2.0_DP*muParam*(2.0_DP*EXP(x(1)-x(2))+EXP(x(2)-x(3)))
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=-2.0_DP*muParam*(2.0_DP*EXP(x(3)-x(1))+EXP(x(2)-x(3)))
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_3)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Sine/cosine function
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=sin(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(2)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=2.0_DP*cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)*cos(2.0_DP*PI*x(2)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=-cos(2.0_DP*PI*x(1)/10.0_DP)*sin(2.0_DP*PI*x(2)/10.0_DP)*cos(2.0_DP*PI*x(3)/10.0_DP)
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=-COS(2.0_DP*PI*x(1)/10.0_DP)*(-12.0_DP*muParam*PI*SIN(2.0_DP*PI*x(2)/10.0_DP)* &
                & SIN(2.0_DP*PI*x(3)/10.0_DP)-rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)*10.0_DP+ &
                & 2.0_DP*rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*x(3)/10.0_DP)**2- &
                & rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)*10.0_DP*COS(2.0_DP*PI*x(2)/10.0_DP)**2)/10.0_DP/2.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=0.0_DP
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=36*muParam*PI**2/10.0_DP**2*cos(2.0_DP*PI*x(2)/10.0_DP)*sin(2.0_DP*PI*x(3)/10.0_DP)* &
                & cos(2.0_DP*PI*x(1)/10.0_DP)
            ELSE IF(componentNumber==3) THEN
              !calculate w
              VALUE=0.0_DP
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=0.0_DP
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF

    CASE(EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_4,EQUATIONS_SET_NAVIER_STOKES_EQUATION_THREE_DIM_5)
      IF(numberOfDimensions==3.AND.numberOfComponents==4) THEN
        !Taylor-Green vortex solution
        muParam = materialsParameters(1)
        rhoParam = materialsParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
          CASE(NO_GLOBAL_DERIV)
            IF(componentNumber==1) THEN
              !calculate u
              VALUE=SIN(x(1)/10.0_DP*PI)*COS(x(2)/10.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
            ELSE IF(componentNumber==2) THEN
              !calculate v
              VALUE=-COS(x(1)/10.0_DP*PI)*SIN(x(2)/10.0_DP*PI)*EXP(-2.0_DP*muParam/rhoParam*currentTime)
            ELSE IF(componentNumber==3) THEN
              !calculate v
              VALUE=0.0_DP
              !VALUE=-COS(x(1))*SIN(x(2))
            ELSE IF(componentNumber==4) THEN
              !calculate p
              VALUE=rhoParam/4.0_DP*(COS(2.0_DP*x(1)/10.0_DP*PI)+COS(2.0_DP*x(2)/10.0_DP*PI))* &
                & EXP(-4.0_DP*muParam/rhoParam*currentTime)
              !VALUE=rhoParam/4.0_DP*(COS(2.0_DP*x(1))+COS(2.0_DP*x(2)))
            ELSE
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivIndex)
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
              CALL FlagError("Not implemented.",err,error,*999)
            ENDIF
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivIndex,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The number of components does not correspond to the number of dimensions."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("NavierStokes_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("NavierStokes_AnalyticFunctionsEvaluate",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Update SUPG parameters for Navier-Stokes equation
  SUBROUTINE NavierStokes_ResidualBasedStabilisation(equationsSet,elementNumber,gaussNumber,mu,rho,jacobianFlag,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number
    INTEGER(INTG), INTENT(IN) :: gaussNumber !<The gauss point number
    REAL(DP), INTENT(IN) :: mu !<The dynamic viscosity
    REAL(DP), INTENT(IN) :: rho !<The fluid density
    LOGICAL, INTENT(IN) ::  jacobianFlag !<Flag indicating whether this was called from the jacobian or residual evaluation routine
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,crossStress,dimensionIdx,dimensionIdx2, &
      & esSpecification(3),i,j,k,l,numberOfDimensions,meshComponent1,meshComponent2,numberOfDependentComponents, &
      & numberOfElementParameters(4),numberOfXi,pressureIndex,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx, &
      & stabilisationType,variableType,velocityInterpolationOrder(4),xiIdx,xiIdx2
    REAL(DP) :: C1,columnPhi,columndPhi2dXi(3,3),doubleDotG,dPhidXPressure(27,3),dPhidXVelocity(27,3),dPhidXi,dXidX(3,3), &
      & elementInverse,jacobian,jacobianContinuity,jacobianGaussWeight,jacobianMomentum(3),LSIC,meshVelocity(3),momentumTerm, &
      & nuLSIC,pressure,pressureDeriv(3),pressureGaussWeight,PSPG,residualContinuity,residualMomentum(3),reynoldsStress,rowPhi, &
      & stabilisationValueDP,sum,sum2,SUPG,tauC,tauMp,tauMu,tauSUPS,timeIncrement,traceG,uDotGu,velocity(3),velocityDeriv(3,3), &
      & velocityGaussWeight,velocityPrevious(3),velocity2Deriv(3,3,3)
    LOGICAL :: linearElement
    TYPE(BasisType), POINTER :: basisVelocity,basisPressure
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: domain1,domain2
    TYPE(DomainElementsType), POINTER :: domainElements1,domainElements2
    TYPE(DomainTopologyType), POINTER :: domainTopology1,domainTopology2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,geometricField,independentField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters, &
      & independentInterpParameters,prevDependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,independentInterpPoint, &
      & prevDependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(QuadratureSchemeType), POINTER :: quadratureVelocity,quadraturePressure
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("NavierStokes_ResidualBasedStabilisation",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)      
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(equationsSetField)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      NULLIFY(independentField)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) &
        & CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(dependentVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
      NULLIFY(jacobianMatrix)
      CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,1,jacobianMatrix,err,error,*999)

      !Set general and specific pointers
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
      CALL FieldVariable_ComponentMeshComponentGet(dependentVariable,1,meshComponent1,err,error,*999)
      CALL FieldVariable_ComponentMeshComponentGet(dependentVariable,numberOfDependentComponents,meshComponent2,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(domain1)
      CALL FieldVariable_ComponentDomainGet(dependentVariable,1,domain1,err,error,*999)
      NULLIFY(domainTopology1)
      CALL Domain_DomainTopologyGet(domain1,domainTopology1,err,error,*999)
      NULLIFY(domainElements1)
      CALL DomainTopology_DomainElementsGet(domainTopology1,domainElements1,err,error,*999)
      NULLIFY(basisVelocity)
      CALL DomainElements_ElementBasisGet(domainElements1,elementNumber,basisVelocity,err,error,*999)
      CALL Basis_NumberOfXiGet(basisVelocity,numberOfXi,err,error,*999)
      NULLIFY(domain2)
      CALL FieldVariable_ComponentDomainGet(dependentVariable,numberOfDependentComponents,domain2,err,error,*999)
      NULLIFY(domainTopology2)
      CALL Domain_DomainTopologyGet(domain2,domainTopology2,err,error,*999)
      NULLIFY(domainElements2)
      CALL DomainTopology_DomainElementsGet(domainTopology2,domainElements2,err,error,*999)
      NULLIFY(basisPressure)
      CALL DomainElements_ElementBasisGet(domainElements2,elementNumber,basisPressure,err,error,*999)

      CALL Basis_InterpolationOrderGet(basisVelocity,velocityInterpolationOrder,err,error,*999)
      IF(velocityInterpolationOrder(1)<=BASIS_LINEAR_INTERPOLATION_ORDER) THEN
        linearElement = .TRUE.
      ELSE
        !Higher order element type- can calculate 2nd order terms
        linearElement = .FALSE.
      ENDIF

      NULLIFY(quadratureVelocity)
      CALL Basis_QuadratureSchemeGet(basisVelocity,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureVelocity,err,error,*999)
      NULLIFY(quadraturePressure)
      CALL Basis_QuadratureSchemeGet(basisPressure,BASIS_DEFAULT_QUADRATURE_SCHEME,quadraturePressure,err,error,*999)

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,variableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,dependentInterpPoint,err,error,*999)
      NULLIFY(prevDependentInterpParameters)
      NULLIFY(prevDependentInterpPoint)
      IF(esSpecification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_PreviousDependentParametersGet(equationsInterpolation,variableType, &
          & prevDependentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_PreviousDependentPointGet(equationsInterpolation,variableType, &
          & prevDependentInterpPoint,err,error,*999)
      ENDIF
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
       ENDIF

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
        IF(esSpecification(3)/=EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE.AND.timeIncrement < ZERO_TOLERANCE) THEN
          CALL FlagError("Please set the equations set field time increment to a value > 0.",err,error,*999)
        ENDIF
        ! Stabilisation type (default 1 for RBS)
        CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,4,stabilisationValueDP, &
          & err,error,*999)
        stabilisationType=NINT(stabilisationValueDP)
        !User specified or previously calculated C1
        CALL Field_ParameterSetGetLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,elementNumber,10, &
          & elementInverse,err,error,*999)

        !Get previous timestep values
        velocityPrevious=0.0_DP
        IF(esSpecification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber, &
            & prevDependentInterpParameters,err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber, &
            & prevDependentInterpPoint,err,error,*999)
          velocityPrevious(1:numberOfDimensions)=prevDependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        !Interpolate current solution velocity/pressure field values
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters, &
          & err,error,*999)
        IF(linearElement) THEN
          !Get 1st order derivatives for current timestep value
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,dependentInterpPoint, &
            & err,error,*999)
        ELSE
          !Get 2nd order derivatives for current timestep value
          CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,dependentInterpPoint, &
            & err,error,*999)
        ENDIF
        velocity=0.0_DP
        velocityDeriv=0.0_DP
        velocity2Deriv=0.0_DP
        pressure=0.0_DP
        pressureDeriv=0.0_DP
        DO dimensionIdx=1,numberOfDimensions
          velocity(dimensionIdx)=dependentInterpPoint%values(dimensionIdx,NO_PART_DERIV)
          velocityDeriv(dimensionIdx,1)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1)
          velocityDeriv(dimensionIdx,2)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S2)
          IF(.NOT.linearElement) THEN
            velocity2Deriv(dimensionIdx,1,1)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1_S1)
            velocity2Deriv(dimensionIdx,1,2)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1_S2)
            velocity2Deriv(dimensionIdx,2,1)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1_S2)
            velocity2Deriv(dimensionIdx,2,2)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S2_S2)
          ENDIF
          IF(numberOfDimensions > 2) THEN
            velocityDeriv(dimensionIdx,3)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S3)
            IF(.NOT. linearElement) THEN
              velocity2Deriv(dimensionIdx,1,3)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1_S3)
              velocity2Deriv(dimensionIdx,2,3)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S2_S3)
              velocity2Deriv(dimensionIdx,3,1)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1_S3)
              velocity2Deriv(dimensionIdx,3,2)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S2_S3)
              velocity2Deriv(dimensionIdx,3,3)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S3_S3)
            ENDIF
          ENDIF
        ENDDO !dimensionIdx
        pressureIndex = numberOfDimensions + 1
        pressure=dependentInterpPoint%values(pressureIndex,NO_PART_DERIV)
        pressureDeriv(1)=dependentInterpPoint%values(pressureIndex,PART_DERIV_S1)
        pressureDeriv(2)=dependentInterpPoint%values(pressureIndex,PART_DERIV_S2)
        IF(numberOfDimensions > 2) pressureDeriv(3)=dependentInterpPoint%values(pressureIndex,PART_DERIV_S3)
        dXidX=0.0_DP
        DO dimensionIdx=1,numberOfDimensions
          DO xiIdx=1,numberOfXi
            dXidX(xiIdx,dimensionIdx)=geometricInterpPointMetrics%dXidX(xiIdx,dimensionIdx)
          ENDDO !xiIdx
        ENDDO !dimensionIdx
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(quadratureVelocity,gaussNumber,velocityGaussWeight,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(quadraturePressure,gaussNumber,pressureGaussWeight,err,error,*999)

        meshVelocity=0.0_DP
        IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
            & err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,independentInterpPoint,err, &
            & error,*999)
          meshVelocity(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        !Get number of element parameters for each dependent component
        CALL Basis_NumberOfElementParametersGet(basisVelocity,numberOfElementParameters(1),err,error,*999)
        numberOfElementParameters(numberOfDimensions+1:numberOfDimensions)=numberOfElementParameters(1)
        CALL Basis_NumberOfElementParametersGet(basisPressure,numberOfElementParameters(numberOfDimensions+1),err,error,*999)
        !Calculate dPhi/dX
        dPhidXVelocity=0.0_DP
        dPhidXPressure=0.0_DP
        DO rowElementParameterIdx=1,numberOfElementParameters(1)
          DO dimensionIdx=1,numberOfDimensions
            dPhidXVelocity(rowElementParameterIdx,dimensionIdx)=0.0_DP
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussNumber,dPhidXi,err,error,*999)
              dPhidXVelocity(rowElementParameterIdx,dimensionIdx)=dPhidXVelocity(rowElementParameterIdx,dimensionIdx) + &
                & dPhidXi*dXidX(xiIdx,dimensionIdx)
            ENDDO !xiIdx
          ENDDO !dimensionIdx
        ENDDO !rowElementParameterIdx
        DO rowElementParameterIdx=1,numberOfElementParameters(numberOfDimensions+1)
          DO dimensionIdx=1,numberOfDimensions
            dPhidXPressure(rowElementParameterIdx,dimensionIdx)=0.0_DP
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadraturePressure,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussNumber,dPhidXi,err,error,*999)
              dPhidXPressure(rowElementParameterIdx,dimensionIdx)=dPhidXPressure(rowElementParameterIdx,dimensionIdx) + &
                & dPhidXi*dXidX(xiIdx,dimensionIdx)
            ENDDO !xiIdx
          ENDDO !dimensionIdx
        ENDDO !rowElementParameterIdx
        jacobianGaussWeight=jacobian*velocityGaussWeight

        !----------------------------------------------------------------------------------
        ! C a l c u l a t e   d i s c r e t e   r e s i d u a l s
        !----------------------------------------------------------------------------------
        sum = 0.0_DP
        residualMomentum = 0.0_DP
        residualContinuity = 0.0_DP
        ! Calculate momentum residual
        DO dimensionIdx=1,numberOfDimensions
          sum = 0.0_DP
          ! velocity time derivative
          IF(esSpecification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
!!TODO: Should interpolate previous mesh velocity so that the delta velocity should be
!!      ((velocity - meshVelocity) - (previousVelocity - previousMeshVelocity)) however
!!      mesh velocity doesn't really change that much and so previousMeshVelocity ~ meshVelocity
!!      and so it will cancel out.
            sum = rho*(velocity(dimensionIdx)-velocityPrevious(dimensionIdx))/timeIncrement
          ENDIF
!! TODO: j is a xi index for the pressure gradient term below but a dimension idx for the viscous stress term below???
          DO j=1,numberOfDimensions
            ! pressure gradient
            sum = sum + pressureDeriv(j)*dXidX(j,dimensionIdx)
            DO k=1,numberOfDimensions
              !Convective term
              sum = sum +rho*((velocity(j)-meshVelocity(j))*(velocityDeriv(dimensionIdx,k)*dXidX(k,j)))
              IF(.NOT. linearElement) THEN
                DO l=1,numberOfDimensions
                  ! viscous stress: only if quadratic or higher basis defined for laplacian
                  sum = sum - mu*(velocity2Deriv(dimensionIdx,k,l)*dXidX(k,j)*dXidX(l,j))
                ENDDO !l
              ENDIF
            ENDDO !k
          ENDDO !j
          residualMomentum(dimensionIdx) = sum
        ENDDO !dimensionIdx
        ! Calculate continuity residual
        sum = 0.0_DP
        DO dimensionIdx=1,numberOfDimensions
          DO xiIdx=1,numberOfXi
            sum= sum + velocityDeriv(dimensionIdx,xiIdx)*dXidX(xiIdx,dimensionIdx)
          ENDDO !xiIdx
        ENDDO !dimensonIdx
        residualContinuity = sum

        ! Constant of element inverse inequality
        IF(elementInverse > -ZERO_TOLERANCE) THEN
          ! Use user-defined value if specified (default -1)
          C1 = elementInverse
        ELSE IF(linearElement) THEN
          C1=3.0_DP
        ELSE
          IF(numberOfDimensions==2.AND.numberOfElementParameters(1)==9.AND. &
            & velocityInterpolationOrder(1)==BASIS_QUADRATIC_INTERPOLATION_ORDER) THEN
            C1=24.0_DP
          ELSE IF(numberOfDimensions==3.AND.numberOfElementParameters(1)==27.AND. &
            & velocityInterpolationOrder(1)==BASIS_QUADRATIC_INTERPOLATION_ORDER) THEN
            C1=12.0_DP
            !TODO: Expand C1 for more element types
          ELSE
            localError="Element inverse estimate undefined on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF
        ! Update element inverse value if calculated
        IF(ABS(C1-elementInverse) > ZERO_TOLERANCE) THEN
          CALL Field_ParameterSetUpdateLocalElement(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & elementNumber,10,C1,err,error,*999)
          ! Should probably move this field update it only happens one time for each element, when C1 undefined
!!TODO: CHANGE THIS. IT WILL CURRENTLY UPDATE FOR EVERY GAUSS POINT FOR EVERY ELEMENT.
          CALL Field_ParameterSetUpdateStart(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSetField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & err,error,*999)
        ENDIF

        !----------------------------------------------------------
        ! S t a b i l i z a t i o n    C o n s t a n t s    (Taus)
        !----------------------------------------------------------
        IF(stabilisationType == 1 .OR. stabilisationType == 2) THEN
          !Bazilevs method for calculating tau
          uDotGu = 0.0_DP
          DO i=1,numberOfDimensions
            DO j=1,numberOfDimensions
              uDotGu = uDotGu + (velocity(i)-meshVelocity(i))*geometricInterpPointMetrics%gu(i,j)*(velocity(j)-meshVelocity(j))
            ENDDO !j
          ENDDO !i
          doubleDotG = 0.0_DP
          DO i=1,numberOfDimensions
            DO j=1,numberOfDimensions
              doubleDotG = doubleDotG + geometricInterpPointMetrics%gu(i,j)*geometricInterpPointMetrics%gu(i,j)
            ENDDO !i
          ENDDO !j
          !Calculate tauSUPS (used for both PSPG and SUPG weights)
          IF(esSpecification(3) == EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
            tauSUPS = (uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
          ELSE
            tauSUPS = ((4.0_DP/(timeIncrement**2.0_DP)) + uDotGu + (C1*((mu/rho)**2.0_DP)*doubleDotG))**(-0.5_DP)
          ENDIF

          ! Calculate nu_LSIC (Least-squares incompressibility constraint)
          CALL Trace(geometricInterpPointMetrics%gu(1:numberOfDimensions,1:numberOfDimensions),traceG,err,error,*999)
          nuLSIC = 1.0_DP/(tauSUPS*traceG)

          tauMp = tauSUPS
          tauMu = tauSUPS
          tauC = nuLSIC

        ELSE
          localError="A tau factor has not been defined for the stabilisation type of "// &
            & TRIM(NumberToVString(stabilisationType,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF

        !-------------------------------------------------------------------------------------------------
        ! A d d   s t a b i l i z a t i o n   f a c t o r s   t o   e l e m e n t   m a t r i c e s
        !-------------------------------------------------------------------------------------------------
        jacobianMomentum = 0.0_DP
        jacobianContinuity = 0.0_DP
        rowElementDOFIdx = 0
        DO rowComponentIdx=1,numberOfDimensions+1
          DO rowElementParameterIdx=1,numberOfElementParameters(rowComponentIdx)
            rowElementDOFIdx = rowElementDOFIdx + 1
            IF(rowComponentIdx <= numberOfDimensions) THEN
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussNumber,rowPhi,err,error,*999)
              jacobianGaussWeight=jacobian*velocityGaussWeight
            ELSE
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadraturePressure,rowElementParameterIdx,NO_PART_DERIV, &
                & gaussNumber,rowPhi,err,error,*999)
              jacobianGaussWeight=jacobian*pressureGaussWeight
            ENDIF
            !------------------
            ! J A C O B I A N
            !------------------
            IF(jacobianFlag) THEN
              columnElementDOFIdx = 0
              DO columnComponentIdx=1,numberOfDimensions+1
                DO columnElementParameterIdx=1,numberOfElementParameters(columnComponentIdx)
                  columnElementDOFIdx=columnElementDOFIdx+1
                  ! Note that we still need to assemble the vector momentum jacobian for PSPG in the continuity row
                  IF(columnComponentIdx <= numberOfDimensions) THEN
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx,NO_PART_DERIV, &
                      & gaussNumber,columnPhi,err,error,*999)
                  ELSE
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadraturePressure,columnElementParameterIdx,NO_PART_DERIV, &
                      & gaussNumber,columnPhi,err,error,*999)
                  ENDIF

                  !Calculate jacobians of the discrete residual terms
                  jacobianMomentum = 0.0_DP
                  IF(columnComponentIdx == numberOfDimensions+1) THEN
                    ! d(Momentum(rowComponentIdx))/d(Pressure)
                    DO dimensionIdx=1,numberOfDimensions
                      jacobianMomentum(dimensionIdx) = dPhidXPressure(columnElementParameterIdx,dimensionIdx)
                    ENDDO !dimensionIdx
                    jacobianContinuity=0.0_DP
                  ELSE
                    columndPhi2dXi=0.0_DP
                    IF(.NOT. linearElement) THEN
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                        & PART_DERIV_S1_S1,gaussNumber,columndPhi2dXi(1,1),err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                        & PART_DERIV_S1_S2,gaussNumber,columndPhi2dXi(1,2),err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                        & PART_DERIV_S1_S2,gaussNumber,columndPhi2dXi(2,1),err,error,*999)
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                        & PART_DERIV_S2_S2,gaussNumber,columndPhi2dXi(2,2),err,error,*999)
                      IF(numberOfDimensions > 2) THEN
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                          & PART_DERIV_S1_S3,gaussNumber,columndPhi2dXi(1,3),err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                          & PART_DERIV_S2_S3,gaussNumber,columndPhi2dXi(2,3),err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                          & PART_DERIV_S1_S3,gaussNumber,columndPhi2dXi(3,1),err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                          & PART_DERIV_S2_S3,gaussNumber,columndPhi2dXi(3,2),err,error,*999)
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,columnElementParameterIdx, &
                          & PART_DERIV_S3_S3,gaussNumber,columndPhi2dXi(3,3),err,error,*999)
                      ENDIF
                    ENDIF
                    !d(Momentum)/d(Velocity(columnComponentIdx))
                    jacobianMomentum = 0.0_DP
                    DO dimensionIdx=1,numberOfDimensions
                      sum = 0.0_DP
                      !Note: Convective term split using product rule
                      !Convective term 1: applies to all velocity components
                      DO xiIdx=1,numberOfXi
                        sum = sum + rho*columnPhi*velocityDeriv(dimensionIdx,xiIdx)*dXidX(xiIdx,columnComponentIdx)
                      ENDDO !xiIdx
                      !Diagonal terms
                      IF(dimensionIdx==columnComponentIdx) THEN
                        !Transient
                        IF(esSpecification(3) /= EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE) THEN
                          sum = sum + rho*columnPhi/timeIncrement
                        ENDIF
                        !Convective 2: columnComponentIdx component only
                        DO dimensionIdx2=1,numberOfDimensions
                          sum = sum + rho*(velocity(dimensionIdx2)-meshVelocity(dimensionIdx2))* &
                            & dPhidXVelocity(columnElementParameterIdx,dimensionIdx2)
                        ENDDO !dimensionIdx2
                        IF(.NOT. linearElement) THEN
                          !Viscous laplacian term
                          DO dimensionIdx2=1,numberOfDimensions
                            DO xiIdx=1,numberOfXi
                              DO xiIdx2=1,numberOfXi
                                sum=sum-mu*columndPhi2dXi(xiIdx,xiIdx2)*dXidX(xiIdx,dimensionIdx2)*dXidX(xiIdx2,dimensionIdx2)
                              ENDDO !xiIdx2
                            ENDDO !xiIdx
                          ENDDO !dimensionIdx2
                        ENDIF
                      ENDIF
                      jacobianMomentum(dimensionIdx)=sum
                    ENDDO !dimensionIdx
                    ! Continuity/velocity
                    jacobianContinuity = dPhidXVelocity(columnElementParameterIdx,columnComponentIdx)
                  ENDIF
                  ! Calculate jacobian of discrete residual * RBS factors (apply product rule if neccesary)

                  IF(rowComponentIdx == numberOfDimensions+1) THEN
                    ! PSPG: Pressure stabilising Petrov-Galerkin
                    PSPG = 0.0_DP
                    sum = 0.0_DP
                    DO dimensionIdx=1,numberOfDimensions
                      sum = sum + dPhidXPressure(rowElementParameterIdx,dimensionIdx)*jacobianMomentum(dimensionIdx)
                    ENDDO !dimensionIdx
                    PSPG = tauMp*sum/rho*jacobianGaussWeight

                    jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+PSPG

                  ELSE
                    ! SUPG: Streamline upwind/Petrov-Galerkin
                    ! LSIC: Least-squares incompressibility constraint
                    SUPG=0.0_DP
                    LSIC=0.0_DP

                    sum=0.0_DP
                    IF(columnComponentIdx <= numberOfDimensions) THEN
                      SUPG= SUPG + columnPhi*dPhidXVelocity(rowElementParameterIdx,columnComponentIdx)* &
                        & residualMomentum(rowComponentIdx)
                    ENDIF
                    DO dimensionIdx=1,numberOfDimensions
                      sum = sum + (velocity(dimensionIdx)-meshVelocity(dimensionIdx))* &
                        & dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
                    ENDDO !dimensionIdx
                    SUPG = tauMu*(SUPG + sum*jacobianMomentum(rowComponentIdx))

                    sum=0.0_DP
                    DO dimensionIdx=1,numberOfDimensions
                      sum = sum + dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
                    ENDDO !dimensionIdx
                    LSIC = tauC*rho*dPhidXVelocity(rowElementParameterIdx,rowComponentIdx)*jacobianContinuity

                    momentumTerm = (SUPG + LSIC)*jacobianGaussWeight

                    IF(stabilisationType == 2) THEN
                      ! Additional terms for RBVM
                      crossStress=0.0_DP
                      reynoldsStress=0.0_DP
                      crossStress = 0.0_DP
                      IF(columnComponentIdx <= numberOfDimensions) THEN
                        IF(rowComponentIdx == columnComponentIdx) THEN
                          DO dimensionIdx=1,numberOfDimensions
                            crossStress= crossStress + dPhidXVelocity(columnElementParameterIdx,dimensionIdx)* &
                              & residualMomentum(dimensionIdx)
                          ENDDO !dimensionIdx
                        ENDIF
                      ENDIF
                      sum2=0.0_DP
                      DO dimensionIdx=1,numberOfDimensions
                        sum=0.0_DP
                        ! dU_rowComponentIdx/dX_i
                        DO xiIdx=1,numberOfXi
                          sum= sum + velocityDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,dimensionIdx)
                        ENDDO !xiIdx
                        ! Jm_i*dU_rowComponentIdx/dX_i
                        sum2 = sum2 + jacobianMomentum(dimensionIdx)*sum
                      ENDDO !dimensionIdx
                      crossStress = -tauMu*(crossStress + sum2)

                      reynoldsStress = 0.0_DP
                      sum = 0.0_DP
                      !Rm_rowComponentIdx.Rm_i.dPhi/dX_i
                      DO dimensionIdx=1,numberOfDimensions
                        sum = sum + jacobianMomentum(rowComponentIdx)*residualMomentum(dimensionIdx)* &
                          & dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
                        sum = sum + jacobianMomentum(dimensionIdx)*residualMomentum(rowComponentIdx)* &
                          & dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
                      ENDDO !dimensionIdx
                      reynoldsStress = -tauMu*tauMu*sum

                      momentumTerm = momentumTerm + (crossStress + reynoldsStress)*jacobianGaussWeight
                    ENDIF

                    ! Add stabilisation to element jacobian
                    jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+momentumTerm

                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx

              !-----------------
              ! R E S I D U A L
              !-----------------
            ELSE
              IF(rowComponentIdx == numberOfDimensions+1) THEN
                ! PSPG: Pressure stabilising Petrov-Galerkin
                sum = 0.0_DP
                DO dimensionIdx=1,numberOfDimensions
                  sum = sum + dPhidXPressure(rowElementParameterIdx,dimensionIdx)*residualMomentum(dimensionIdx)
                ENDDO !dimensionIdx
                PSPG = sum*(tauMp/rho)*jacobianGaussWeight
                residualVector%elementResidual%vector(rowElementDOFIdx)= &
                  & residualVector%elementResidual%vector(rowElementDOFIdx) + PSPG

              ELSE
                ! SUPG: Streamline upwind/Petrov-Galerkin
                ! LSIC: Least-squares incompressibility constraint
                SUPG=0.0_DP
                LSIC=0.0_DP

                ! u_i*Rm_mh*dv_mh/dx_i
                sum=0.0_DP
                DO dimensionIdx=1,numberOfDimensions
                  sum=sum+(velocity(dimensionIdx)-meshVelocity(dimensionIdx))*dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
                ENDDO !dimensionIdx
                SUPG = tauMu*sum*residualMomentum(rowComponentIdx)

                LSIC = tauC*rho*dPhidXVelocity(rowElementParameterIdx,rowComponentIdx)*residualContinuity
                momentumTerm = (SUPG + LSIC)*jacobianGaussWeight

                IF(stabilisationType ==2) THEN
                  ! Additional terms for RBVM
                  crossStress=0.0_DP
                  reynoldsStress=0.0_DP
                  sum2 = 0.0_DP
                  DO dimensionIdx=1,numberOfDimensions
                    sum = 0.0_DP
                    ! dU_mh/dX_i
                    DO xiIdx=1,numberOfXi
                      sum= sum + velocityDeriv(rowComponentIdx,xiIdx)*dXidX(xiIdx,dimensionIdx)
                    ENDDO !xiIdx
                    ! Rm_i.dU_mh/dX_i
                    sum2= sum2 + residualMomentum(dimensionIdx)*sum
                  ENDDO !dimensionIdx
                  crossStress= -tauMu*rowPhi*sum2

                  reynoldsStress = 0.0_DP
                  sum = 0.0_DP
                  !Rm_mh.Rm_i.dPhi/dX_i
                  DO dimensionIdx=1,numberOfDimensions
                    sum = sum + dPhidXVelocity(rowElementParameterIdx,dimensionIdx)*residualMomentum(dimensionIdx)* &
                      & residualMomentum(rowComponentIdx)
                  ENDDO !dimensionIdx
                  reynoldsStress = -sum*(tauMu*tauMu)/rho
                  momentumTerm = momentumTerm + (crossStress + reynoldsStress)*jacobianGaussWeight
                ENDIF

                ! Add stabilisation to element residual
                residualVector%elementResidual%vector(rowElementDOFIdx)= &
                  & residualVector%elementResidual%vector(rowElementDOFIdx) + momentumTerm
              ENDIF
            ENDIF ! jacobian/residual
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx

      ENDIF ! check stabilisation type
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dimensionIdx,esSpecification(3),gaussNumber,info,lWork,numberOfDependentComponents, &
      & numberOfElementParameters,numberOfDimensions,numberOfGauss,numberOfXi,rowComponentIdx,rowElementParameterIdx, &
      & variableType,xiIdx
    REAL(DP) :: avgVelocity(3),cellReynoldsNumber,cellCourantNumber,CMatrix(27,3),currentTime,dPhidXVelocity(27,3), &
      & dRowPhidXi,dXidX(3,3),gaussWeight,jacobian,jacobianGaussWeight,KMatrix(27,3),meshVelocity(3),MMatrix(27,3),mu,muScale, &
      & normCMatrix,normKMatrix,normMMatrix,rho,rowPhi,sum,sum2,svd(3),timeIncrement,U(27,27),velocity(3),velocityDeriv(3,3), &
      & velocityNorm,velocityPrevious(3),VT(3,3)
    REAL(DP), ALLOCATABLE :: work(:)
    TYPE(BasisType), POINTER :: basisVelocity
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: equationsSetField
    TYPE(QuadratureSchemeType), POINTER :: quadratureVelocity
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters, &
      & independentInterpParameters,prevDependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,independentInterpPoint, &
      & prevDependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable,uMaterialsVariable,vEquationsSetVariable,vMaterialsVariable
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("NavierStokes_CalculateElementMetrics",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
      NULLIFY(equationsSetField)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      NULLIFY(vEquationsSetVariable)
      CALL Field_VariableGet(equationsSetField,FIELD_V_VARIABLE_TYPE,vEquationsSetVariable,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      NULLIFY(vMaterialsVariable)
      IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) &
        & CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
      NULLIFY(independentField)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) &
        & CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      !Set general and specific pointers
      NULLIFY(dependentVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(domain)
      CALL FieldVariable_ComponentDomainGet(dependentVariable,1,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
      NULLIFY(basisVelocity)
      CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basisVelocity,err,error,*999)
      CALL Basis_NumberOfXiGet(basisVelocity,numberOfXi,err,error,*999)
      CALL Basis_NumberOfElementParametersGet(basisVelocity,numberOfElementParameters,err,error,*999)
      NULLIFY(quadratureVelocity)
      CALL Basis_QuadratureSchemeGet(basisVelocity,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureVelocity,err,error,*999)
 
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
       NULLIFY(geometricInterpPointMetrics)
       CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
         & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
        & err,error,*999)
      NULLIFY(prevDependentInterpParameters)
      CALL EquationsInterpolation_PreviousDependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & prevDependentInterpParameters,err,error,*999)
      NULLIFY(prevDependentInterpPoint)
      CALL EquationsInterpolation_PreviousDependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & prevDependentInterpPoint,err,error,*999)
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
         CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
      ENDIF
        
      ! Get time step size
      CALL EquationsSet_TimesGet(equationsSet,currentTime,timeIncrement,err,error,*999)
      !CALL Field_ParameterSetGetConstant(equationsSetField,FIELD_U1_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      ! & 3,timeIncrement,err,error,*999)

      ! Loop over gauss points
      CMatrix = 0.0_DP
      MMatrix = 0.0_DP
      KMatrix = 0.0_DP
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      ENDIF
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)

      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,prevDependentInterpParameters, &
        & err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) &
        & CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
        & err,error,*999)

      CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureVelocity,numberOfGauss,err,error,*999)
      
      avgVelocity = 0.0_DP      
      DO gaussNumber = 1,numberOfGauss

        !Get the constitutive law (non-Newtonian) viscosity based on shear rate
        IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
          ! Note the constant from the U_VARIABLE is a scale factor
          muScale = mu
          ! Get the gauss point based value returned from the CellML solver
          CALL FieldVariable_ParameterSetGetLocalGaussPoint(vMaterialsVariable,FIELD_VALUES_SET_TYPE,gaussNumber,elementNumber,1, &
            & mu,err,error,*999)
          mu=mu*muScale
        ENDIF

        !Get previous timestep values
        velocityPrevious=0.0_DP
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,prevDependentInterpPoint, &
          & err,error,*999)
        velocityPrevious(1:numberOfDimensions)=prevDependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)

        ! Interpolate current solution velocity and first deriv field values
        ! Get 1st order derivatives for current timestep value
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,dependentInterpPoint, &
          & err,error,*999)
        velocity=0.0_DP
        velocityDeriv=0.0_DP
        DO dimensionIdx=1,numberOfDimensions
          velocity(dimensionIdx)=dependentInterpPoint%values(dimensionIdx,NO_PART_DERIV)
          velocityDeriv(dimensionIdx,1)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S1)
          velocityDeriv(dimensionIdx,2)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S2)
          IF(numberOfDimensions > 2) velocityDeriv(dimensionIdx,3)=dependentInterpPoint%values(dimensionIdx,PART_DERIV_S3)
        ENDDO !dimensionIdx

        meshVelocity=0.0_DP
        IF(esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,independentInterpPoint, &
            & err,error,*999)
          meshVelocity(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        ! get dXi/dX deriv
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussNumber,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        dXidX=0.0_DP
        dXidX(1:numberOfDimensions,1:numberOfDimensions)= &
          & geometricInterpPointMetrics%dXidX(1:numberOfDimensions,1:numberOfDimensions)

        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(quadratureVelocity,gaussNumber,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        ! Calculate dPhi/dX
        dPhidXVelocity=0.0_DP
        DO rowElementParameterIdx=1,numberOfElementParameters
          DO dimensionIdx=1,numberOfDimensions
            dPhidXVelocity(rowElementParameterIdx,dimensionIdx)=0.0_DP
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussNumber,dRowPhidXi,err,error,*999)
              dPhidXVelocity(rowElementParameterIdx,dimensionIdx)=dPhidXVelocity(rowElementParameterIdx,dimensionIdx) + &
                & dRowPhidXi*dXidX(xiIdx,dimensionIdx)
            ENDDO !xiIdx
          ENDDO !dimensionIdx
        ENDDO !rowElementParameterIdx

        DO rowComponentIdx=1,numberOfDimensions
          DO rowElementParameterIdx=1,numberOfElementParameters
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureVelocity,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussNumber,rowPhi,err,error,*999)

            ! c_(a,i)
            sum=0.0_DP
            DO dimensionIdx=1,numberOfDimensions
              DO xiIdx=1,numberOfXi
                sum = sum + (velocity(dimensionIdx)-meshVelocity(dimensionIdx))*velocityDeriv(rowComponentIdx,xiIdx)* &
                  & dXidX(xiIdx,dimensionIdx)
              ENDDO !xiIdx
            ENDDO !dimensionIdx
            CMatrix(rowElementParameterIdx,rowComponentIdx)=CMatrix(rowElementParameterIdx,rowComponentIdx) + &
              & rho*rowPhi*sum*jacobianGaussWeight

            ! ~k_(a,i)
            sum=0.0_DP
            DO dimensionIdx=1,numberOfDimensions
              sum = sum + (velocity(dimensionIdx)-meshVelocity(dimensionIdx))*dPhidXVelocity(rowElementParameterIdx,dimensionIdx)
            ENDDO !dimensionIdx
            sum2=0.0_DP
            DO dimensionIdx=1,numberOfDimensions
              DO xiIdx=1,numberOfXi
                sum2 = sum2 + (velocity(dimensionIdx)-meshVelocity(dimensionIdx))*velocityDeriv(rowComponentIdx,xiIdx)* &
                  & dXidX(xiIdx,dimensionIdx)
              ENDDO !xiIdx
            ENDDO !dimensionIdx
            KMatrix(rowElementParameterIdx,rowComponentIdx)=KMatrix(rowElementParameterIdx,rowComponentIdx)+ &
              & rho*sum*sum2*jacobianGaussWeight

            ! m_(a,i)
!!TODO: Should interpolate previous mesh velocity so that the delta velocity should be
!!      ((velocity - meshVelocity) - (previousVelocity - previousMeshVelocity)) however
!!      mesh velocity doesn't really change that much and so previousMeshVelocity ~ meshVelocity
!!      and so it will cancel out.
            MMatrix(rowElementParameterIdx,rowComponentIdx)=MMatrix(rowElementParameterIdx,rowComponentIdx)+ &
              & rho*rowPhi*(velocity(rowComponentIdx)-velocityPrevious(rowComponentIdx))/timeIncrement*jacobianGaussWeight

          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx

        avgVelocity= avgVelocity + (velocity-meshVelocity)/quadratureVelocity%numberOfGauss
        
      ENDDO ! gauss loop

      lWork=MAX(1,3*MIN(numberOfElementParameters,numberOfDimensions)+ &
        & MAX(numberOfElementParameters,numberOfDimensions),5*MIN(numberOfElementParameters,numberOfDimensions))
      ALLOCATE(work(lWork))

      ! compute the singular value decomposition (SVD) using LAPACK
      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,CMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,work,lWork,info)
      normCMatrix=svd(1)
      IF(info /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,KMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,work,lWork,info)
      normKMatrix=svd(1)
      IF(info /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF

      CALL DGESVD('A','A',numberOfElementParameters,numberOfDimensions,MMatrix,numberOfElementParameters,svd, &
        & U,numberOfElementParameters,VT,numberOfDimensions,work,lWork,info)
      normMMatrix=svd(1)
      IF(info /= 0) THEN
        localError="Error calculating SVD on element "//TRIM(NumberToVString(elementNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DEALLOCATE(work)

      CALL L2Norm(avgVelocity,velocityNorm,err,error,*999)
      cellReynoldsNumber = 0.0_DP
      cellCourantNumber = 0.0_DP
      IF(velocityNorm > ZERO_TOLERANCE) THEN
        IF(normKMatrix > ZERO_TOLERANCE) THEN
          cellReynoldsNumber = velocityNorm**2.0_DP/(mu/rho)*normCMatrix/normKMatrix
        ENDIF
        IF(normMMatrix > ZERO_TOLERANCE) THEN
          cellCourantNumber = timeIncrement/2.0_DP*normCMatrix/normMMatrix
        ENDIF
      ENDIF
      CALL FieldVariable_ParameterSetUpdateLocalElement(vEquationsSetVariable,FIELD_VALUES_SET_TYPE,elementNumber,2, &
        & velocityNorm,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateLocalElement(vEquationsSetVariable,FIELD_VALUES_SET_TYPE,elementNumber,3, &
        & cellCourantNumber,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateLocalElement(vEquationsSetVariable,FIELD_VALUES_SET_TYPE,elementNumber,4, &
        & cellReynoldsNumber,err,error,*999)

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<The equations set to calculate the RHS term for
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate the RHS term for
    TYPE(FieldVariableType), POINTER :: dependentVariable
    LOGICAL, INTENT(IN) ::  jacobianFlag !<Flag indicating whether this was called from the jacobian or residual evaluation routine
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: boundaryIdx,boundaryNumber,boundaryType,columnComponentIdx,columnElementBaseDOFIdx,columnElementDOFIndex, &
      & columnElementParameterIdx,columnNodeIdx,columnNodeDerivativeIdx,columnParameterIdx,componentIdx,dependentBasisType1, &
      & dependentBasisType2,dependentVariableType,elementNodeIdx,elementNodeIdx2,esSpecification(3),gaussIdx, &
      & globalNodeDerivativeIdx,globalNodeDerivativeIdx2,lineXiDirection,numberOfDependentComponents,numberOfDimensions, &
      & numberOfElementBoundaries,numberOfElementParameters1,numberOfElementParameters2,numberOfGauss,numberOfNodes1, &
      & numberOfNodes2,numberOfNodeDerivatives,numberOfXi2,orientation,rowComponentIdx,rowElementBaseDOFIdx,rowElementDOFIndex, &
      & rowElementParameterIdx,rowNodeIdx,rowNodeDerivativeIdx,rowParameterIdx,xiDirection(4),xiIdx,xiNormalDirection
    REAL(DP) :: beta,boundaryNormal(3),boundaryValue,columnPhi,columndPhidXi(3),density,gaussWeight,jacobian, &
      & jacobianGaussWeights,meshVelocity(3),normalDifference,normalFlow,normalProjection(3),normalTolerance,pressure,rowPhi, &
      & stabilisationTerm,unitNormal(3),velocity(3)
    LOGICAL :: boundaryFace,boundaryLine,calculateFaces,calculateLines,integratedBoundary
    TYPE(BasisType), POINTER :: basis1,basis2,dependentBasis1,dependentBasis2
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionFaceType), POINTER :: face
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces
    TYPE(DecompositionLineType), POINTER :: line
    TYPE(DecompositionLinesType), POINTER :: decompositionLines
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain1,domain2
    TYPE(DomainElementsType), POINTER :: domainElements1,domainElements2
    TYPE(DomainFacesType), POINTER :: domainFaces1,domainFaces2
    TYPE(DomainLinesType), POINTER :: domainLines1,domainLines2
    TYPE(DomainTopologyType), POINTER :: domainTopology1,domainTopology2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: geometricField,equationsSetField,dependentField,independentField,materialsField,variableField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters, &
      & independentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,independentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: geometricVariable,u1EquationsSetVariable,uMaterialsVariable,vEquationsSetVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme1,quadratureScheme2
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_FiniteElementBoundaryIntegrate",err,error,*999)

    ! Get pointers and perform sanity checks
    IF(.NOT.ASSOCIATED(dependentVariable)) CALL FlagError("Dependent variable is not associated.",err,error,*999)
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)

      CALL FieldVariable_VariableTypeGet(dependentVariable,dependentVariableType,err,error,*999)
      NULLIFY(equationsSetField)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      NULLIFY(vEquationsSetVariable)
      CALL Field_VariableGet(equationsSetField,FIELD_V_VARIABLE_TYPE,vEquationsSetVariable,err,error,*999)
      NULLIFY(u1EquationsSetVariable)
      CALL Field_VariableGet(equationsSetField,FIELD_U1_VARIABLE_TYPE,u1EquationsSetVariable,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
      NULLIFY(jacobianMatrix)
      IF(jacobianFlag) CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,1,jacobianMatrix,err,error,*999)
      NULLIFY(independentField)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) &
        & CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
 
      ! Check whether this element contains an integrated boundary type
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable,FIELD_VALUES_SET_TYPE,elementNumber,9, &
        & boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      integratedBoundary = .FALSE.
      IF(boundaryType == BOUNDARY_CONDITION_PRESSURE) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_FIXED_PRESSURE) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_COUPLING_STRESS) integratedBoundary = .TRUE.
      IF(boundaryType == BOUNDARY_CONDITION_FIXED_CELLML) integratedBoundary = .TRUE.

      !Get the mesh decomposition and basis
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
      NULLIFY(domain1)
      CALL FieldVariable_ComponentDomainGet(dependentVariable,1,domain1,err,error,*999)
      NULLIFY(domainTopology1)
      CALL Domain_DomainTopologyGet(domain1,domainTopology1,err,error,*999)
      NULLIFY(domainElements1)
      CALL DomainTopology_DomainElementsGet(domainTopology1,domainElements1,err,error,*999)
      NULLIFY(domainFaces1)
      CALL DomainTopology_DomainFacesGet(domainTopology1,domainFaces1,err,error,*999)
      NULLIFY(domainLines1)
      CALL DomainTopology_DomainLinesGet(domainTopology1,domainLines1,err,error,*999)
      NULLIFY(dependentBasis1)
      CALL DomainElements_ElementBasisGet(domainElements1,elementNumber,dependentBasis1,err,error,*999)
      CALL Basis_TypeGet(dependentBasis1,dependentBasisType1,err,error,*999)
      CALL Basis_NumberOfElementParametersGet(dependentBasis1,numberOfElementParameters1,err,error,*999)
      NULLIFY(domain2)
      CALL FieldVariable_ComponentDomainGet(dependentVariable,numberOfDependentComponents,domain2,err,error,*999)
      NULLIFY(domainTopology2)
      CALL Domain_DomainTopologyGet(domain2,domainTopology2,err,error,*999)
      NULLIFY(domainElements2)
      CALL DomainTopology_DomainElementsGet(domainTopology2,domainElements2,err,error,*999)
      NULLIFY(domainFaces2)
      CALL DomainTopology_DomainFacesGet(domainTopology2,domainFaces2,err,error,*999)
      NULLIFY(domainLines2)
      CALL DomainTopology_DomainLinesGet(domainTopology2,domainLines2,err,error,*999)
      NULLIFY(dependentBasis2)
      CALL DomainElements_ElementBasisGet(domainElements2,elementNumber,dependentBasis2,err,error,*999)
      CALL Basis_TypeGet(dependentBasis2,dependentBasisType2,err,error,*999)
      CALL Basis_NumberOfElementParametersGet(dependentBasis2,numberOfElementParameters2,err,error,*999)
      NULLIFY(variableField)
      CALL FieldVariable_FieldGet(dependentVariable,variableField,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(variableField,decomposition,err,error,*999)
      CALL Decomposition_CalculateFacesGet(decomposition,calculateFaces,err,error,*999)
      CALL Decomposition_CalculateLinesGet(decomposition,calculateLines,err,error,*999)
      NULLIFY(decompositionTopology)
      CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
      NULLIFY(decompositionElements)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology,decompositionElements,err,error,*999)
      NULLIFY(decompositionFaces)
      CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology,decompositionFaces,err,error,*999)
      NULLIFY(decompositionLines)
      CALL DecompositionTopology_DecompositionLinesGet(decompositionTopology,decompositionLines,err,error,*999)
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,dependentVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,dependentVariableType,dependentInterpPoint, &
        & err,error,*999)
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
      ENDIF

      !Determine if this is a 2D or 3D problem with line/face parameters calculated
      IF(integratedBoundary) THEN
        IF(numberOfDimensions/=2.AND.numberOfDimensions/=3) THEN
          localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
            & " is invalid for a 2D or 3D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        IF(numberOfDimensions==3.AND.calculateFaces) THEN
          CALL Basis_NumberOfLocalFacesGet(dependentBasis1,numberOfElementBoundaries,err,error,*999)
        ELSE IF(numberOfDimensions==2.AND.calculateLines) THEN
          CALL Basis_NumberOfLocalLinesGet(dependentBasis1,numberOfElementBoundaries,err,error,*999)
        ELSE
          integratedBoundary = .FALSE.
        ENDIF
      ENDIF

      !Only apply to required element boundaries
      IF(integratedBoundary) THEN
        !Get the density
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,density,err,error,*999)
        ! Get the boundary element parameters
        CALL FieldVariable_ParameterSetGetConstant(u1EquationsSetVariable,FIELD_VALUES_SET_TYPE,1,beta,err,error,*999)
        boundaryNormal = 0.0_DP
        DO componentIdx=1,3
          CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable,FIELD_VALUES_SET_TYPE,elementNumber, &
            & componentIdx+4,boundaryNormal(componentIdx),err,error,*999)
        ENDDO !componentIdx

        ! Loop over the boundaries (lines or faces) for this element
        DO boundaryIdx=1,numberOfElementBoundaries

          IF(numberOfDimensions == 3) THEN
            !Get 3D face specific parameters
            CALL DecompositionElements_ElementFaceNumberGet(decompositionElements,boundaryIdx,elementNumber,boundaryNumber, &
              & err,error,*999)
            CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces,boundaryNumber,boundaryFace,err,error,*999)
            
            !This speeds things up but is also important, as non-boundary faces have an xi direction that might
            !correspond to the other element.
            IF(.NOT.(boundaryFace)) CYCLE
            
            SELECT CASE(dependentBasisType1)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces,boundaryNumber,xiNormalDirection,err,error,*999)
              xiDirection(3)=ABS(xiNormalDirection)
              xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
              xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
              orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*xiNormalDirection)
            CASE(BASIS_SIMPLEX_TYPE)
              orientation=1
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType1,"*",err,error))// &
                & " is not yet implemented for Navier-Stokes boundary integration."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            NULLIFY(basis1)
            CALL DomainFaces_FaceBasisGet(domainFaces1,boundaryNumber,basis1,err,error,*999)
            NULLIFY(basis2)
            CALL DomainFaces_FaceBasisGet(domainFaces2,boundaryNumber,basis2,err,error,*999)
            CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,boundaryNumber,geometricInterpParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,boundaryNumber,dependentInterpParameters, &
              & err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,boundaryNumber,independentInterpParameters, &
                & err,error,*999)
            ENDIF

          ELSE IF(numberOfDimensions==2) THEN
            !Get 2D line specific parameters
            CALL DecompositionElements_ElementLineNumberGet(decompositionElements,boundaryIdx,elementNumber,boundaryNumber, &
              & err,error,*999)
            CALL DecompositionLines_LineBoundaryLineGet(decompositionLines,boundaryNumber,boundaryLine,err,error,*999)
            
            !This speeds things up but is also important, as non-boundary lines have an xi direction that might
            !correspond to the other element.
            IF(.NOT.(boundaryLine)) CYCLE
            
            SELECT CASE(dependentBasisType1)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              CALL DecompositionLines_LineXiDirectionGet(decompositionLines,boundaryNumber,lineXiDirection,err,error,*999)
              xiDirection(2)=ABS(lineXiDirection)
              xiDirection(1)=OTHER_XI_DIRECTIONS2(xiDirection(2))
              orientation=SIGN(1,OTHER_XI_ORIENTATIONS2(xiDirection(1))*lineXiDirection)
            CASE(BASIS_SIMPLEX_TYPE)
              orientation=1
            CASE DEFAULT
              localError="Line integration for basis type "//TRIM(NumberToVString(dependentBasisType1,"*",err,error))// &
                & " is not yet implemented for Navier-Stokes boundary integration."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            NULLIFY(basis1)
            CALL DomainLines_LineBasisGet(domainLines1,boundaryNumber,basis1,err,error,*999)
            NULLIFY(basis2)
            CALL DomainLines_LineBasisGet(domainLines2,boundaryNumber,basis2,err,error,*999)
            CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,boundaryNumber,geometricInterpParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,boundaryNumber,dependentInterpParameters, &
              & err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,boundaryNumber,independentInterpParameters, &
                & err,error,*999)
            ENDIF
          ENDIF

          CALL Basis_NumberOfLocalNodesGet(basis1,numberOfNodes1,err,error,*999)
          CALL Basis_NumberOfLocalNodesGet(basis2,numberOfNodes2,err,error,*999)
          CALL Basis_NumberOfXiGet(basis2,numberOfXi2,err,error,*999)
          NULLIFY(quadratureScheme1)
          CALL Basis_QuadratureSchemeGet(basis1,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme1,err,error,*999)
          NULLIFY(quadratureScheme2)
          CALL Basis_QuadratureSchemeGet(basis2,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme2,err,error,*999)
          
          !Loop over gauss points
          CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme1,numberOfGauss,err,error,*999)
          DO gaussIdx=1,numberOfGauss
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,geometricInterpPoint, &
              & err,error,*999)
            normalTolerance=0.1_DP
            unitNormal = 0.0_DP
            velocity = 0.0_DP
            meshVelocity = 0.0_DP
            pressure = 0.0_DP
            IF(numberOfDimensions == 3) THEN
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_AREA_TYPE,geometricInterpPointMetrics, &
                & err,error,*999)
              !Make sure this is the boundary that corresponds with the provided normal (could be a wall rather than inlet/outlet)
              CALL CrossProduct(geometricInterpPointMetrics%dXdXi(:,1),geometricInterpPointMetrics%dXdXi(:,2),normalProjection, &
                & err,error,*999)
              normalProjection = normalProjection*orientation
              CALL Normalise(normalProjection,unitNormal,err,error,*999)
              CALL L2Norm(boundaryNormal-unitNormal,normalDifference,err,error,*999)
              IF(normalDifference>normalTolerance) EXIT
              CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,boundaryNumber,dependentInterpParameters, &
                & err,error,*999)
            ELSE
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_LINE_TYPE,geometricInterpPointMetrics, &
                & err,error,*999)
              !Make sure this is the boundary that corresponds with the provided normal (could be a wall rather than inlet/outlet)
              normalProjection = [geometricInterpPointMetrics%dXdXi(2,1),geometricInterpPointMetrics%dXdXi(1,1),0.0_DP]*orientation
              CALL Normalise(normalProjection,unitNormal,err,error,*999)
              CALL L2Norm(boundaryNormal-unitNormal,normalDifference,err,error,*999)
              IF(normalDifference>normalTolerance) EXIT
              CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,boundaryNumber,dependentInterpParameters, &
                & err,error,*999)
            ENDIF
            !Jacobian and Gauss weighting term
            CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
            CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme1,gaussIdx,gaussWeight,err,error,*999)
            jacobianGaussWeights=jacobian*gaussWeight

            !Get interpolated velocity and pressure
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,dependentInterpPoint,err,error,*999)
            velocity(1:numberOfDimensions)=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            IF(esSpecification(3)==EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE) THEN
              !Get interpolated mesh velocity
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,independentInterpPoint, &
                & err,error,*999)
              meshVelocity(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
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
            ENDIF

            ! Check for Neumann integrated boundary types rather than fixed pressure types
            IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
              & boundaryType==BOUNDARY_CONDITION_FIXED_CELLML .OR. &
              & boundaryType==BOUNDARY_CONDITION_PRESSURE) THEN
              !Get the pressure value interpolation parameters
             IF(numberOfDimensions==3) THEN
               CALL Field_InterpolationParametersFaceGet(FIELD_PRESSURE_VALUES_SET_TYPE,boundaryNumber, &
                 & dependentInterpParameters,err,error,*999)
              ELSE
                CALL Field_InterpolationParametersLineGet(FIELD_PRESSURE_VALUES_SET_TYPE,boundaryNumber, &
                  & dependentInterpParameters,err,error,*999)
              ENDIF
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,dependentInterpPoint, &
                & err,error,*999)
              pressure=dependentInterpPoint%values(numberOfDependentComponents,NO_PART_DERIV)
            ENDIF

            !Loop over field components
            DO rowComponentIdx=1,numberOfDimensions
              !Work out the first index of the vector for this element - (i.e. the number of previous)
              rowElementBaseDOFIdx=numberOfElementParameters1*(rowComponentIdx-1)
              DO rowNodeIdx=1,numberOfNodes1
                IF(numberOfDimensions==3) THEN
                  CALL Basis_FaceNodeNumberGet(dependentBasis1,rowNodeIdx,boundaryIdx,elementNodeIdx,err,error,*999)
                ELSE
                  CALL Basis_LineNodeNumberGet(dependentBasis1,rowNodeIdx,boundaryIdx,elementNodeIdx,err,error,*999)
                ENDIF
                CALL Basis_NodeNumberOfDerivativesGet(basis1,rowNodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO rowNodeDerivativeIdx=1,numberOfNodeDerivatives
                  IF(numberOfDimensions==3) THEN
                    CALL Basis_FaceNodeDerivativeNumberGet(dependentBasis1,rowNodeDerivativeIdx,rowNodeIdx,boundaryIdx, &
                      & globalNodeDerivativeIdx,err,error,*999)
                  ELSE
                    CALL Basis_LineNodeDerivativeNumberGet(dependentBasis1,rowNodeIdx,boundaryIdx,globalNodeDerivativeIdx, &
                      & err,error,*999)
                  ENDIF
                  CALL Basis_ElementParameterGet(dependentBasis1,globalNodeDerivativeIdx,elementNodeIdx,rowElementParameterIdx, &
                    & err,error,*999)
                  CALL Basis_ElementParameterGet(basis1,rowNodeDerivativeIdx,rowNodeIdx,rowParameterIdx,err,error,*999)
                  rowElementDOFIndex=rowElementBaseDOFIdx+rowElementParameterIdx
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme1,rowParameterIdx,NO_PART_DERIV,gaussIdx, &
                    & rowPhi,err,error,*999)

                  IF(.NOT.jacobianFlag) THEN
                    IF(boundaryType==BOUNDARY_CONDITION_PRESSURE .OR. &
                      & boundaryType==BOUNDARY_CONDITION_FIXED_CELLML .OR. &
                      & boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
                      !Integrated boundary pressure term
                      residualVector%elementResidual%vector(rowElementDOFIndex)= &
                        & residualVector%elementResidual%vector(rowElementDOFIndex)-&
                        &  pressure*unitNormal(rowComponentIdx)*rowPhi*jacobianGaussWeights
                    ENDIF
                    !Boundary stabilisation term (if necessary )
                    IF(ABS(beta) > ZERO_TOLERANCE) THEN
                      residualVector%elementResidual%vector(rowElementDOFIndex)= &
                        & residualVector%elementResidual%vector(rowElementDOFIndex)-&
                        & 0.5_DP*beta*density*rowPhi*(velocity(rowComponentIdx)-meshVelocity(rowComponentIdx))* &
                        & stabilisationTerm*jacobianGaussWeights
                    ENDIF
                  !Jacobian matrix term is the derivative of the nonlinear stabilisation term
                  !Loop over field components
                  ELSE
                    IF(ABS(beta) > ZERO_TOLERANCE) THEN
                      DO columnComponentIdx=1,numberOfDimensions
                        !Work out the first index of the rhs vector for this element - (i.e. the number of previous)
                        columnElementBaseDOFIdx=numberOfElementParameters2*(columnComponentIdx-1)
                        DO columnNodeIdx=1,numberOfNodes2
                          IF(numberOfDimensions == 3) THEN
                            CALL Basis_FaceNodeNumberGet(dependentBasis2,columnNodeIdx,boundaryIdx,elementNodeIdx2,err,error,*999)
                          ELSE
                            CALL Basis_LineNodeNumberGet(dependentBasis2,columnNodeIdx,boundaryIdx,elementNodeIdx2,err,error,*999)
                          ENDIF
                          CALL Basis_NodeNumberOfDerivativesGet(basis2,columnNodeIdx,numberOfNodeDerivatives,err,error,*999)
                          DO columnNodeDerivativeIdx=1,numberOfNodeDerivatives
                            IF(numberOfDimensions==3) THEN
                              CALL Basis_FaceNodeDerivativeNumberGet(dependentBasis2,columnNodeDerivativeIdx,columnNodeIdx, &
                                & boundaryIdx,globalNodeDerivativeIdx2,err,error,*999)
                            ELSE
                              CALL Basis_LineNodeDerivativeNumberGet(dependentBasis2,columnNodeIdx,boundaryIdx, &
                                & globalNodeDerivativeIdx,err,error,*999)
                              globalNodeDerivativeIdx2=dependentBasis2%derivativeNumbersInLocalLine(columnNodeIdx,boundaryIdx)
                            ENDIF
                            CALL Basis_ElementParameterGet(dependentBasis2,globalNodeDerivativeIdx2,elementNodeIdx2, &
                              & columnElementParameterIdx,err,error,*999)
                            CALL Basis_ElementParameterGet(basis2,columnNodeDerivativeIdx,columnNodeIdx,columnParameterIdx, &
                              & err,error,*999)
                            columnElementDOFIndex=columnElementBaseDOFIdx+columnElementParameterIdx
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme2,columnParameterIdx,NO_PART_DERIV, &
                              & gaussIdx,columnPhi,err,error,*999)
                            columnPhi=quadratureScheme2%gaussBasisFunctions(columnParameterIdx,NO_PART_DERIV,gaussIdx)
                            DO xiIdx=1,numberOfXi2
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme2,columnParameterIdx, &
                                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussIdx,columndPhidXi(xiIdx),err,error,*999)
                            ENDDO !xiIdx
                            !Jacobian term
                            IF(rowComponentIdx == columnComponentIdx) THEN
                              ! note that (u_j.n_j - |u_j.n_j|) term derivative will be zero
                              jacobianMatrix%elementJacobian%matrix(rowElementDOFIndex,columnElementDOFIndex)= &
                                & jacobianMatrix%elementJacobian%matrix(rowElementDOFIndex,columnElementDOFIndex) - &
                                & 0.5_DP*beta*density*rowPhi*columnPhi*stabilisationTerm*jacobianGaussWeights
                            ENDIF
                          ENDDO !columnNodeDerivativeIdx
                        ENDDO !columnNodeIdx
                      ENDDO !columnComponentIdx
                    ENDIF
                  ENDIF

                ENDDO !rowNodeDerivativeIdx
              ENDDO !rowNodeIdx
            ENDDO !rowComponentIdx
          ENDDO !gaussIdx
        ENDDO !boundaryIdx
      ENDIF

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

  !>Calculate the fluid flux through 3D boundaries for use in problems with coupled solutions (e.g. multidomain)
  SUBROUTINE NavierStokes_CalculateBoundaryFlux(equationsSet,coupledEquationsSet,iteration3D1D, &
    & convergedFlag,absolute3D0DTolerance,relative3D0DTolerance,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsSetType), POINTER :: coupledEquationsSet !<A pointer to the coupled equations set (for 3D-1D coupling)
    INTEGER(INTG), INTENT(IN) :: iteration3D1D !<iteration number for the 3D-1D loop if this is a coupled problem
    REAL(DP), INTENT(IN) :: absolute3D0DTolerance !<absolute convergence criteria for 3D-0D coupling
    REAL(DP), INTENT(IN) :: relative3D0DTolerance !<relative convergence criteria for 3D-0D coupling
    LOGICAL, INTENT(INOUT) :: convergedFlag !<convergence flag for 3D-0D coupling
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryID,boundaryType,componentIdx,computationNode,coupledNodeNumber,dependentBasisType, &
      & dependentBasisType2,dependentVariableType,elementIdx,elementNodeIdx,elementUserNumber,esSpecification(3),faceIdx, &
      & faceNodeIdx,faceNodeDerivativeIdx,faceNumber,gaussIdx,groupCommunicator,i,j,meshComponentNumber,mpiIError,nodeNumber, &
      & normalComponentIdx,numberOfBoundaries,numberOfDimensions,numberOfGauss,numberOfGlobalBoundaries, &
      & numberOfGroupComputationNodes,numberOfLocalElements,numberOfLocalFaces,numberOfLocalFaces2,numberOfNodes, &
      & numberOfNodeDerivatives,totalNumberOfLocal,versionNumber,xiNormalDirection
    REAL(DP) :: a1D,boundaryValue,cauchy(3,3),couplingFlow,couplingStress,courant,dUdXi(3,3),dXidX(3,3),elementNormal(3), &
      & faceArea,faceNormal(3),facePressure,faceTraction,faceVelocity,flowError,gaussWeight,globalBoundaryArea(10), &
      & globalBoundaryFlux(10),globalBoundaryMeanPressure(10),globalBoundaryMeanNormalStress(10),globalBoundaryNormalStress(10), &
      & globalBoundaryPressure(10),gradU(3,3),jacobian,localBoundaryArea(10),localBoundaryFlux(10),localBoundaryNormalStress(10), &
      & localBoundaryPressure(10),maxCourant,mu,muScale,normal,normalDifference,normalProjection,normalTolerance,normalWave(2), &
      & p0D,p1D,p3D,pressureError,pressureGauss,q0D,q1D,rho,stress1DPrevious,toleranceCourant,traction(3),unitNormal(3), &
      & velocityGauss(3),viscousTerm(3)
    LOGICAL :: boundary3D0DFound(10),boundaryFace,couple1DTo3D,couple3DTo1D
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(BasisType), POINTER :: dependentBasis,dependentBasis2,faceBasis
    TYPE(DecompositionType), POINTER :: decomposition3D
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(DecompositionElementsType), POINTER :: decompositionElements3D
    TYPE(DecompositionFaceType), POINTER :: face
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces3D
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology3D
    TYPE(DomainType), POINTER :: domain,domain3D
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainMappingType), POINTER :: elementsMapping3D
    TYPE(DomainMappingsType), POINTER :: domainMappings3D
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations3D
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField1D,dependentField3D,equationsSetField3D,geometricField3D, &
      & independentField1D,materialsField1D,materialsField3D
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable3D,geometricVariable3D,uDependentVariable1D,uDependentVariable3D, &
      & u1DependentVariable1D,u2DependentVariable1D,independentVariable1D,materialsVariable1D,uMaterialsVariable3D, &
      & uEquationsSetVariable3D,u1EquationsSetVariable3D,vEquationsSetVariable3D,vMaterialsVariable3D
    TYPE(QuadratureSchemeType), POINTER :: faceQuadratureScheme
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup
 
    ENTERS("NavierStokes_CalculateBoundaryFlux",err,error,*999)

    couple1DTo3D = .FALSE.
    couple3DTo1D = .FALSE.
    boundary3D0DFound = .FALSE.

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
 
    SELECT CASE(esSpecification(3))
    ! 3 D   t y p e s :   I n t e g r a t e   b o u n d a r y   v a l u e s
    ! ------------------------------------------------------------------------
    CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)

      !Get 3D field pointers
      NULLIFY(equations3D)
      CALL EquationsSet_EquationsGet(equationsSet,equations3D,err,error,*999)
      NULLIFY(geometricField3D)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField3D,err,error,*999)
      NULLIFY(geometricVariable3D)
      CALL Field_VariableGet(geometricField3D,FIELD_U_VARIABLE_TYPE,geometricVariable3D,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable3D,numberOfDimensions,err,error,*999)
      NULLIFY(dependentField3D)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField3D,err,error,*999)
      NULLIFY(materialsField3D)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField3D,err,error,*999)
      NULLIFY(uMaterialsVariable3D)
      CALL Field_VariableGet(materialsField3D,FIELD_U_VARIABLE_TYPE,uMaterialsVariable3D,err,error,*999)
      NULLIFY(vMaterialsVariable3D)
      CALL Field_VariableGet(materialsField3D,FIELD_V_VARIABLE_TYPE,vMaterialsVariable3D,err,error,*999)
      NULLIFY(equationsSetField3D)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField3D,err,error,*999)
      NULLIFY(u1EquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_U1_VARIABLE_TYPE,u1EquationsSetVariable3D,err,error,*999)
      NULLIFY(vEquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_V_VARIABLE_TYPE,vEquationsSetVariable3D,err,error,*999)
      ! Check for a coupled equations set
      NULLIFY(dependentField1D)
      
      NULLIFY(independentField1D)
      IF(ASSOCIATED(coupledEquationsSet)) THEN
        CALL EquationsSet_DependentFieldGet(coupledEquationsSet,dependentField1D,err,error,*999)
        CALL EquationsSet_IndependentFieldGet(coupledEquationsSet,independentField1D,err,error,*999)
      ENDIF
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations3D,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(dependentVariable3D)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable3D,err,error,*999)
      CALL FieldVariable_VariableTypeGet(dependentVariable3D,dependentVariableType,err,error,*999)
      !Get the mesh decomposition and mapping
      NULLIFY(decomposition3D)
      CALL Field_DecompositionGet(dependentField3D,decomposition3D,err,error,*999)
      NULLIFY(decompositionTopology3D)
      CALL Decomposition_DecompositionTopologyGet(decomposition3D,decompositionTopology3D,err,error,*999)
      NULLIFY(decompositionElements3D)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology3D,decompositionElements3D,err,error,*999)
      NULLIFY(decompositionFaces3D)
      CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology3D,decompositionFaces3D,err,error,*999)
      NULLIFY(domain3D)
      CALL Decomposition_DomainGet(decomposition3D,0,domain3D,err,error,*999)
      NULLIFY(domainMappings3D)
      CALL Domain_DomainMappingsGet(domain3D,domainMappings3D,err,error,*999)
      NULLIFY(elementsMapping3D)
      CALL DomainMappings_ElementsMappingGet(domainMappings3D,elementsMapping3D,err,error,*999)
      CALL DomainMapping_NumberOfLocalGet(elementsMapping3D,numberOfLocalElements,err,error,*999)

      !Get interpolations
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations3d,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,dependentVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,dependentVariableType,dependentInterpPoint, &
        & err,error,*999)
      
      !Get constant max Courant (CFL) number (default 1.0)
      CALL FieldVariable_ParameterSetGetConstant(u1equationsSetVariable3D,FIELD_VALUES_SET_TYPE,2,toleranceCourant,err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable3D,FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable3D,FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      ENDIF

      !Loop over elements to locate boundary elements
      maxCourant = 0.0_DP
      numberOfBoundaries = 0
      localBoundaryFlux = 0.0_DP
      localBoundaryArea = 0.0_DP
      localBoundaryPressure = 0.0_DP
      localBoundaryNormalStress = 0.0_DP
      DO elementIdx=1,numberOfLocalElements
        CALL FieldVariable_ComponentMeshComponentGet(dependentVariable3D,1,meshComponentNumber,err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable3D,1,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(domainFaces)
        CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        NULLIFY(dependentBasis)
        CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis,err,error,*999)
        CALL Basis_TypeGet(dependentBasis,dependentBasisType,err,error,*999)
        CALL Basis_NumberOfLocalFacesGet(dependentBasis,numberOfLocalFaces,err,error,*999)
        CALL DecompositionElements_ElementUserNumberGet(decompositionElements3D,elementIdx,elementUserNumber,err,error,*999)

        !Note: if CFL tolerance = 0, we'll skip this step, which speeds things up a bit
        IF(toleranceCourant > ZERO_TOLERANCE) THEN
          ! C F L  c o n d i t i o n   c h e c k
          ! ------------------------------------
          !Calculate element metrics (courant #, cell Reynolds number)
          CALL NavierStokes_CalculateElementMetrics(equationsSet,elementIdx,err,error,*999)
          !Get element metrics
          CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,3,courant, &
            & err,error,*999)
          IF(courant < -ZERO_TOLERANCE) CALL FlagWarning("Negative Courant (CFL) number.",err,error,*999)
          IF(courant > maxCourant) maxCourant = courant
          !Check if element CFL number below specified tolerance
          
          IF(courant > toleranceCourant) THEN
            localError="Element "//TRIM(NumberToVString(elementUserNumber,"*",err,error))//" has violated the CFL condition "// &
              & TRIM(NumberToVString(courant,"*",err,error))//" <= "//TRIM(NumberToVString(toleranceCourant,"*",err,error))// &
              & ". Decrease timestep or increase CFL tolerance for the 3D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF

        ! B o u n d a r y   n o r m a l   a n d   I D
        ! ----------------------------------------------
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,5, &
          & elementNormal(1),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,6, &
          & elementNormal(2),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,7, &
          & elementNormal(3),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,8, &
          & boundaryValue,err,error,*999)
        boundaryID=NINT(boundaryValue)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,9, &
          & boundaryValue,err,error,*999)
        boundaryType=NINT(boundaryValue)
        !Check if is a non-wall boundary element
        IF(boundaryID > numberOfBoundaries) numberOfBoundaries=boundaryID
        IF(boundaryID>1) THEN
          faceArea=0.0_DP
          faceVelocity=0.0_DP
          facePressure=0.0_DP
          faceTraction=0.0_DP
          ! Loop over faces to determine the boundary face contribution
          DO faceIdx=1,numberOfLocalFaces

            !Get the face normal and quadrature information
            CALL DecompositionElements_ElementFaceNumberGet(decompositionElements3D,faceIdx,elementIdx,faceNumber,err,error,*999)
            CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces3D,faceNumber,boundaryFace,err,error,*999)
           
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(boundaryFace)) CYCLE

            SELECT CASE(dependentBasisType)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces3D,faceNumber,xiNormalDirection,err,error,*999)
              normalComponentIdx=ABS(xiNormalDirection)
            CASE(BASIS_SIMPLEX_TYPE)
              CALL FlagWarning("Boundary flux calculation not yet set up for simplex element types.",err,error,*999)
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType,"*",err,error))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(localError,err,error,*999)
            END SELECT

            NULLIFY(faceBasis)
            CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)
            NULLIFY(faceQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_NumberOfGaussGet(faceQuadratureScheme,numberOfGauss,err,error,*999)

            ! Loop over face gauss points
            DO gaussIdx=1,numberOfGauss
              !Use the geometric field to find the face normal and Jacobian for the face integral
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters, &
                & err,error,*999)
              CALL Field_InterpolateLocalFaceGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                & geometricInterpPoint,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE,geometricInterpPointMetrics, &
                & err,error,*999)
              CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)              

              !TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
              SELECT CASE(dependentBasisType)
              CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
                ! Make sure this is the boundary face that corresponds with boundaryID (could be a wall rather than inlet/outlet)
                DO componentIdx=1,numberOfDimensions
                  normalProjection=DOT_PRODUCT(geometricInterpPointMetrics%gu(normalComponentIdx,:), &
                    & geometricInterpPointMetrics%dXdXi(componentIdx,:))
                  IF(xiNormalDirection<0) THEN
                    normalProjection=-normalProjection
                  ENDIF
                  faceNormal(componentIdx)=normalProjection
                END DO !componentIdx
                CALL Normalise(faceNormal,unitNormal,err,error,*999)
                CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
                normalTolerance=0.1_DP
                IF(normalDifference>normalTolerance) EXIT
              CASE(BASIS_SIMPLEX_TYPE)
                faceNormal=unitNormal
              CASE DEFAULT
                localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType,"*",err,error))// &
                  & " is not yet implemented for Navier-Stokes."
                CALL FlagError(localError,err,error,*999)
              END SELECT

              ! C a l c u l a t e   C a u c h y   s t r e s s   o n   c o u p l e d    n o d e s
              ! -----------------------------------------------------------------------------------
              CALL BasisQuadratureScheme_GaussWeightGet(faceQuadratureScheme,gaussIdx,gaussWeight,err,error,*999)
              ! Get the first partial velocity derivatives and pressure; calculate the Cauchy stress tensor if necessary
              IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
               & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                couple3DTo1D = .TRUE.
                ! Get the constitutive law (non-Newtonian) viscosity based on shear rate if needed
                IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
                  ! Get the gauss point based value returned from the CellML solver
                  CALL FieldVariable_ParameterSetGetLocalGaussPoint(vMaterialsVariable3D,FIELD_VALUES_SET_TYPE,gaussIdx, &
                    & elementIdx,1,mu,err,error,*999)
                  mu=mu*muScale
                ENDIF
                !Get the pressure and velocity interpolation parameters
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,dependentInterpParameters, &
                  & err,error,*999)
                CALL Field_InterpolateLocalFaceGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                 & dependentInterpPoint,err,error,*999)

                velocityGauss=0.0_DP
                pressureGauss=0.0_DP
                !Interpolated values at gauss point
                velocityGauss=dependentInterpPoint%values(1:3,NO_PART_DERIV)
                pressureGauss=dependentInterpPoint%values(4,NO_PART_DERIV)
                dUdXi = 0.0_DP
                dUdXi(1:3,1)=dependentInterpPoint%values(1:3,PART_DERIV_S1)
                dUdXi(1:3,2)=dependentInterpPoint%values(1:3,PART_DERIV_S2)
                dUdXi(1:3,3)=dependentInterpPoint%values(1:3,PART_DERIV_S3)
                ! Assemble viscous term
                dXidX=0.0_DP
                dXidX=geometricInterpPointMetrics%dXidX(:,:)
                CALL MatrixProduct(dUdXi,dXidX,gradU,err,error,*999)
                cauchy = 0.0_DP
                traction = 0.0_DP
                ! Calculate Cauchy stress tensor
                DO i = 1,3
                  DO j = 1,3
                    ! DEBUG: pressure only?
                    cauchy(i,j) = mu*(gradU(i,j)+gradU(j,i))
                  ENDDO !j
                ENDDO !i
                !DEBUG
                viscousTerm = MATMUL(cauchy,faceNormal)
              ELSE
                !Get interpolated velocity
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,dependentInterpParameters, &
                  & err,error,*999)
                CALL Field_InterpolateLocalFaceGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,gaussIdx, &
                 & dependentInterpPoint,err,error,*999)
                velocityGauss=dependentInterpPoint%values(1:3,NO_PART_DERIV)
                pressureGauss=dependentInterpPoint%values(4,NO_PART_DERIV)
              ENDIF

              ! I n t e g r a t e    f a c e   a r e a ,   v e l o c i t y   a n d   t r a c t i o n
              ! ----------------------------------------------------------------------------------------
              DO componentIdx=1,numberOfDimensions
                faceArea=faceArea + ABS(faceNormal(componentIdx)*gaussWeight*geometricInterpPointMetrics%jacobian)
                faceVelocity=faceVelocity+velocityGauss(componentIdx)*faceNormal(componentIdx)*gaussWeight*jacobian
                facePressure=facePressure+pressureGauss*faceNormal(componentIdx)*gaussWeight*jacobian
                IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
                 & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                  faceTraction = faceTraction + (viscousTerm(componentIdx) - pressureGauss)* &
                    & faceNormal(componentIdx)*gaussWeight*jacobian
                ENDIF
              ENDDO !componentIdx
            ENDDO !gaussIdx
          ENDDO !faceIdx
          localBoundaryFlux(boundaryID) = localBoundaryFlux(boundaryID) + faceVelocity
          localBoundaryArea(boundaryID) = localBoundaryArea(boundaryID) + faceArea
          localBoundaryPressure(boundaryID) = localBoundaryPressure(boundaryID) + facePressure
          localBoundaryNormalStress(boundaryID) = localBoundaryNormalStress(boundaryID)+ faceTraction
        ENDIF !boundaryIdentifier
      END DO !elementIdx
      ! Distribute any updated element fields
      CALL FieldVariable_ParameterSetUpdateStart(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)

      ! G a t h e r   v a l u e s   o v e r   t h r e a d s
      ! ------------------------------------------------------
      ! Need to add boundary flux for any boundaries split accross computation nodes
      numberOfGlobalBoundaries = 0
      globalBoundaryFlux = 0.0_DP
      globalBoundaryArea = 0.0_DP
      globalBoundaryPressure = 0.0_DP
      globalBoundaryNormalStress = 0.0_DP
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition3D,workGroup,err,error,*999)
      CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
      CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
      CALL WorkGroup_GroupNodeNumberGet(workGroup,computationNode,err,error,*999)
      IF(numberOfGroupComputationNodes>1) THEN !use mpi
        CALL MPI_ALLREDUCE(localBoundaryFlux,globalBoundaryFlux,10,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryArea,globalBoundaryArea,10,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryNormalStress,globalBoundaryNormalStress,10,MPI_DOUBLE_PRECISION,MPI_SUM,  &
	 & groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryPressure,globalBoundaryPressure,10,MPI_DOUBLE_PRECISION,MPI_SUM,  &
	 & groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(numberOfBoundaries,numberOfGlobalBoundaries,1,MPI_INTEGER,MPI_MAX,  &
	 & groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
      ELSE
        numberOfGlobalBoundaries = numberOfBoundaries
        globalBoundaryFlux = localBoundaryFlux
        globalBoundaryArea = localBoundaryArea
        globalBoundaryPressure = localBoundaryPressure
        globalBoundaryNormalStress = localBoundaryNormalStress
      ENDIF
      globalBoundaryArea=ABS(globalBoundaryArea)
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          globalBoundaryMeanNormalStress(boundaryID)=globalBoundaryNormalStress(boundaryID)/globalBoundaryArea(boundaryID)
          globalBoundaryMeanPressure(boundaryID)=globalBoundaryPressure(boundaryID)/globalBoundaryArea(boundaryID)
        ENDIF
      ENDDO !boundaryID
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          IF(diagnostics1) THEN
            IF(computationNode==0) THEN
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
              ENDIF
            ENDIF
          ENDIF
        ELSE
          localError="Zero or negative area boundary detected on boundary "// &
            & TRIM(NumberToVString(boundaryID,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !boundaryID

    ! 1 D   t y p e s :   p r e p a r e   t o   c o p y    a n y   c o u p l e d    v a l u e s
    ! --------------------------------------------------------------------------------------------
    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
      couple1DTo3D = .TRUE.
      ! Get coupled 3D field pointers
      IF(.NOT.ASSOCIATED(coupledEquationsSet)) CALL FlagError("Coupled 3D equations set is not associated.",err,error,*999)

      NULLIFY(equations3D)
      CALL EquationsSet_EquationsGet(coupledEquationsSet,equations3D,err,error,*999)
      NULLIFY(dependentField3D)
      CALL EquationsSet_DependentFieldGet(coupledEquationsSet,dependentField3D,err,error,*999)
      NULLIFY(materialsField3D)
      CALL EquationsSet_MaterialsFieldGet(coupledEquationsSet,materialsField3D,err,error,*999)
      NULLIFY(equationsSetField3D)
      CALL EquationsSet_EquationsSetFieldGet(coupledEquationsSet,equationsSetField3D,err,error,*999)
      NULLIFY(vEquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_V_VARIABLE_TYPE,vEquationsSetVariable3D,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations3D,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(dependentVariable3D)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable3D,err,error,*999)
      CALL FieldVariable_VariableTypeGet(dependentVariable3D,dependentVariableType,err,error,*999)
      NULLIFY(uDependentVariable3D)
      CALL Field_VariableGet(dependentField3D,FIELD_U_VARIABLE_TYPE,uDependentVariable3D,err,error,*999)
      NULLIFY(decomposition3D)
      CALL Field_DecompositionGet(dependentField3D,decomposition3D,err,error,*999)
      NULLIFY(decompositionTopology3D)
      CALL Decomposition_DecompositionTopologyGet(decomposition3D,decompositionTopology3D,err,error,*999)
      NULLIFY(decompositionElements3D)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology3D,decompositionElements3D,err,error,*999)
      NULLIFY(decompositionFaces3D)
      CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology3D,decompositionFaces3D,err,error,*999)
      NULLIFY(domain3D)
      CALL Decomposition_DomainGet(decomposition3D,0,domain3D,err,error,*999)
      NULLIFY(domainMappings3D)
      CALL Domain_DomainMappingsGet(domain3D,domainMappings3D,err,error,*999)
      NULLIFY(elementsMapping3D)
      CALL DomainMappings_ElementsMappingGet(domainMappings3D,elementsMapping3D,err,error,*999)
      CALL DomainMapping_TotalNumberOfLocalGet(elementsMapping3D,totalNumberOfLocal,err,error,*999)

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations3D,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      
      ! Check for a 1D equations set
      IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("1D equations set is not associated.",err,error,*999)
      NULLIFY(dependentField1D)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField1D,err,error,*999)
      NULLIFY(uDependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U_VARIABLE_TYPE,uDependentVariable1D,err,error,*999)
      NULLIFY(u1DependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U1_VARIABLE_TYPE,u1DependentVariable1D,err,error,*999)
      NULLIFY(u2DependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U2_VARIABLE_TYPE,u2DependentVariable1D,err,error,*999)
      NULLIFY(materialsField1D)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField1D,err,error,*999)
      NULLIFY(materialsVariable1D)
      CALL Field_VariableGet(materialsField1D,FIELD_U_VARIABLE_TYPE,materialsVariable1D,err,error,*999)
      NULLIFY(independentField1D)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField1D,err,error,*999)
      NULLIFY(independentVariable1D)
      CALL Field_VariableGet(independentField1D,FIELD_U_VARIABLE_TYPE,independentVariable1D,err,error,*999)
    CASE DEFAULT
      localError="Boundary flux calcluation for equations type "//TRIM(NumberToVString(esSpecification(3),"*", &
        & err,error))//" is not yet implemented for Navier-Stokes."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    ! ------------------------------------------------------------------------------------
    ! C o p y    i n t e g r a t e d   v a l u e s    t o    t a r g e t    f i e l d s
    ! ------------------------------------------------------------------------------------
    convergedFlag = .TRUE.
    ! Loop over elements again to allocate flux terms to boundary nodes
    DO elementIdx=1,totalNumberOfLocal
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,5, &
        & elementNormal(1),err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,6, &
        & elementNormal(2),err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,7, &
        & elementNormal(3),err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,8, &
        & boundaryValue,err,error,*999)
      boundaryID=NINT(boundaryValue)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,9, &
        & boundaryValue,err,error,*999)
      boundaryType=NINT(boundaryValue)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,11, &
        & boundaryValue,err,error,*999)
      coupledNodeNumber=NINT(boundaryValue)
      IF(boundaryID>1) THEN
        meshComponentNumber=2
        
        decompositionElement=>decompositionElements3D%elements(elementIdx)
        
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition3D,meshComponentNumber,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(dependentBasis2)
        CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis2,err,error,*999)
        CALL Basis_TypeGet(dependentBasis2,dependentBasisType2,err,error,*999)
        CALL Basis_NumberOfLocalFacesGet(dependentBasis2,numberOfLocalFaces2,err,error,*999)
        NULLIFY(domainFaces)
        CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)

        ! M a p   3 D - 1 D    b o u n d a r y    v a l u e s
        ! --------------------------------------------------------
        IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS .OR. &
          & boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
          ! Note: the coupled node number here is the user number- not the local number
          normal = 0.0_DP
          DO componentIdx =1,2
            CALL FieldVariable_ParameterSetGetNode(independentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber, &
              & componentIdx,normalWave(componentIdx),err,error,*999)
             IF(ABS(normalWave(componentIdx)) > ZERO_TOLERANCE) normal=normalWave(componentIdx)
          ENDDO !componentIdx
          IF(ABS(normal)<ZERO_TOLERANCE) THEN
            localError="Characteristic normal wave not specified for couping node " &
              & //TRIM(NumberToVString(coupledNodeNumber,"*",err,error))//"."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          IF(couple3DTo1D) THEN
            ! Convert integrated 3D values to 1D flow and pressure.
            couplingStress = globalBoundaryMeanNormalStress(boundaryID)
            couplingFlow = globalBoundaryFlux(boundaryID)
            ! Map the coupling flow and stress values from the 3D equations set to the coupled 1D field
            CALL FieldVariable_ParameterSetUpdateNode(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,1, &
              & couplingFlow,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateNode(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,2, &
              & couplingStress,err,error,*999)
          ELSE IF(couple1DTo3D) THEN
            ! Get the flow and pressure values from the coupled 1D equations set
            CALL FieldVariable_ParameterSetGetNode(uDependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,1,q1D, &
              & err,error,*999)
            CALL FieldVariable_ParameterSetGetNode(uDependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,2,a1D, &
              & err,error,*999)
            CALL FieldVariable_ParameterSetGetNode(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,1,p1D, &
              & err,error,*999)
            !DEBUG
            CALL FieldVariable_ParameterSetEnsureCreated(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
              & err,error,*999)
            CALL FieldVariable_ParameterSetGetNode(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,1,1, &
              & coupledNodeNumber,2,stress1DPrevious,err,error,*999)
            CALL FieldVariable_ParameterSetGetConstant(materialsVariable1D,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
            !Convert 1D flow and pressure to coupling stress and flow
            couplingStress = -p1D
            couplingFlow = q1D
            p3D = couplingStress
            ! Update the coupling flow and stress values from the 3D equations set to the coupled 1D field
            CALL FieldVariable_ParameterSetUpdateNode(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,2, &
              & couplingStress,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateNode(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,1,1,coupledNodeNumber,3, &
              & couplingFlow,err,error,*999)
          ENDIF
        ENDIF

        ! B o u n d a r y   F a c e    N o r m a l s
        ! --------------------------------------------------
        DO faceIdx=1,numberOfLocalFaces
          !Get the face normal and quadrature information
          CALL DecompositionElements_ElementFaceNumberGet(decompositionElements3D,faceIdx,elementIdx,faceNumber,err,error,*999)
          CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces3D,faceNumber,boundaryFace,err,error,*999)
          
          IF(.NOT.(boundaryFace)) CYCLE
          
          NULLIFY(faceBasis)
          CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)
          NULLIFY(faceQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)

          !TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
          SELECT CASE(dependentBasisType2)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces3D,faceNumber,xiNormalDirection,err,error,*999)
            normalComponentIdx=ABS(xiNormalDirection)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
            CALL Field_InterpolateLocalFaceGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,faceIdx,1,geometricInterpPoint, &
              & err,error,*999)
            !Should this be an area jacobian? Maybe volume needed for the projection below?
            CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE,geometricInterpPointMetrics,err,error,*999)
            DO componentIdx=1,numberOfDimensions
              normalProjection=DOT_PRODUCT(geometricInterpPointMetrics%gu(normalComponentIdx,:), &
                & geometricInterpPointMetrics%dXdXi(componentIdx,:))
              IF(xiNormalDirection<0) THEN
                normalProjection=-normalProjection
              ENDIF
              faceNormal(componentIdx)=normalProjection
            ENDDO !componentIdx
            CALL Normalise(faceNormal,unitNormal,err,error,*999)
          CASE(BASIS_SIMPLEX_TYPE)
            !still have faceNormal/unitNormal
          CASE DEFAULT
            localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType2,"*",err,error))// &
              & " is not yet implemented for Navier-Stokes."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
          normalTolerance=0.1_DP
          IF(normalDifference>normalTolerance) CYCLE

          ! U p d a t e    N o d a l   V a l u e s
          ! --------------------------------------------------
          ! Update local nodes with integrated boundary flow values
          CALL Basis_NumberOfLocalNodesGet(faceBasis,numberOfNodes,err,error,*999)
          DO faceNodeIdx=1,numberOfNodes
            CALL Basis_FaceNodeNumberGet(dependentBasis2,faceNodeIdx,faceIdx,elementNodeIdx,err,error,*999)
            CALL DomainElements_ElementNodeGet(domainElements,elementIdx,elementNodeIdx,nodeNumber,err,error,*999)
            CALL Basis_NodeNumberOfDerivativesGet(faceBasis,faceNodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO faceNodeDerivativeIdx=1,numberOfNodeDerivatives
              versionNumber=1
              IF(couple1DTo3D) THEN
                IF(boundaryType==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
                  ! Copy coupling stress from 1D to pressure values type on 3D
                  !CALL Field_ParameterSetUpdateLocalNode(dependentField3D,FIELD_U_VARIABLE_TYPE, &
                  !  & FIELD_VALUES_SET_TYPE,1,1,nodeNumber,4,p3D,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,1,1, &
                    & nodeNumber,4,p3D,err,error,*999)
                ELSE IF(boundaryType==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
                  localError="Coupling flow boundary type not yet implemented for 3D side of coupled 3D-1D problems."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ELSE
                IF(boundaryType==BOUNDARY_CONDITION_FIXED_CELLML) THEN
                  ! Check current values against those passed to the CellML solver
                  CALL FieldVariable_ParameterSetGetLocalNode(uEquationsSetVariable3D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
                    & 1,1,nodeNumber,1,q0D,err,error,*999)
                  CALL FieldVariable_ParameterSetGetLocalNode(uEquationsSetVariable3D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
                    & 1,1,nodeNumber,2,p0D,err,error,*999)
                  flowError = globalBoundaryFlux(boundaryID)-q0D
                  pressureError = globalBoundaryMeanPressure(boundaryID)+p0D
                  IF((ABS(flowError) < absolute3D0DTolerance .AND. ABS(pressureError) < absolute3D0DTolerance)) THEN
                    ! CONVERGED ABSOLUTE TOLERANCE
                  ELSE IF(ABS(flowError/globalBoundaryFlux(boundaryID)) < relative3D0DTolerance .AND. &
                    &  ABS(pressureError/globalBoundaryMeanPressure(boundaryID)) < relative3D0DTolerance) THEN
                    ! CONVERGED RELATIVE TOLERANCE
                  ELSE
                    IF(.NOT.boundary3D0DFound(boundaryID)) THEN
                      IF(diagnostics1) THEN
                        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  flow: ", &
                          & q0D,err,error,*999)
                        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  pressure: ", &
                          & p0D,err,error,*999)
                        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  flow error: ", &
                          & flowError,err,error,*999)
                        CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"  0D boundary ",boundaryID,"  pressure error: ", &
                          & pressureError,err,error,*999)
                      ENDIF
                      boundary3D0DFound(boundaryID) = .TRUE.
                    ENDIF
                    convergedFlag = .FALSE.
                  ENDIF
                ENDIF
                ! Update flow rates on 3D boundary equations set field
                CALL FieldVariable_ParameterSetUpdateLocalNode(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,1,globalBoundaryFlux(boundaryID),err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalNode(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,2,globalBoundaryMeanPressure(boundaryID),err,error,*999)
              ENDIF
            END DO !rowNodeDerivativeIdx
          END DO !faceNodeIdx
        END DO !faceIdx
      ENDIF !boundaryIdentifier
    END DO !elementIdx

    !allocate array for mpi communication
    IF(numberOfGroupComputationNodes>1) THEN !use mpi
      ALLOCATE(globalConverged(numberOfGroupComputationNodes),STAT=ERR) 
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",err,error,*999)
      CALL MPI_ALLGATHER(convergedFlag,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL,groupCommunicator,mpiIError)
      CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
      IF(ALL(globalConverged)) THEN
        convergedFlag = .TRUE.
      ELSE
        convergedFlag = .FALSE.
      ENDIF
    ENDIF

    !Distribute any updated ields
    IF(ASSOCIATED(equationsSetField3D)) THEN
      CALL FieldVariable_ParameterSetUpdateStart(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
      IF(convergedFlag) THEN
        ! If converged, update 0D initial conditions for the next step.
        CALL FieldVariable_ParameterSetsCopy(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE, &
          & 1.0_DP,err,error,*999)
      ENDIF
    ENDIF
    IF(ASSOCIATED(dependentField3D)) THEN
      CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
    ENDIF
    
    IF(ASSOCIATED(dependentField1D)) THEN
      CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDIF

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
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: boundaryConditionType,boundaryNumber,componentIdx,dependentDof,derivativeIdx,groupCommunicator, &
      & inputIteration,iteration,mpiIError,nodeIdx,nodeNumber,numberOfBoundaries,numberOfGroupComputationNodes, &
      & numberOfLocalNodes1D,outputIteration,pSpecification(3),solverNumber,solver1dNavierStokesNumber,timestep,versionIdx
    REAL(DP) :: aError,a0Param,aPrevious,a1d,beta,couplingTolerance,currentTime,eParam,hParam,normalWave(2),pCellML,pPrevious, &
      & q1d,qError,qPrevious,startTime,stopTime,timeIncrement
    LOGICAL :: boundaryNode,boundaryConverged(30),localConverged,MPI_LOGICAL,coupled1D0DBoundary,continueLoop
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: parentLoop,parentLoop2,subLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,materialsField,independentField
    TYPE(FieldParameterSetType), POINTER :: previousParameterSet
    TYPE(FieldVariableType), POINTER :: dependentVariable,uDependentVariable,u1DependentVariable,dynamicVariable, &
      & independentVariable,vMaterialsVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverType), POINTER :: solver1D
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError    
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("NavierStokes_Couple1D0D",err,error,*999)

    IF(.NOT.ASSOCIATED(controlLoop)) CALL FlagError("Control loop is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    continueLoop = .TRUE.
    
    !Get solvers based on the problem type
    CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      ! In the Navier-Stokes/Characteristic subloop, the Navier-Stokes solver should be the second solver
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      IF(solverNumber == solver1dNavierStokesNumber) THEN
        NULLIFY(subLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,subLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
        NULLIFY(solver1D)
        CALL Solvers_SolverGet(solvers,solver1DNavierStokesNumber,solver1D,err,error,*999)
        CALL ControlLoop_IterationNumberGet(controlLoop,iteration,err,error,*999)
        NULLIFY(parentLoop)
        CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
        CALL ControlLoop_CurrentTimeInformationGet(parentLoop,startTime,stopTime,currentTime,timeIncrement, &
          & timestep,outputIteration,inputIteration,err,error,*999)
      ELSE
        localError="The solver number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
         & " does not correspond with the Navier-Stokes solver number for 1D-0D fluid coupling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      ! In the Navier-Stokes/Characteristic subloop, the Navier-Stokes solver should be the second solver
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      IF(solverNumber == solver1dNavierStokesNumber) THEN
        NULLIFY(subLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,2,subLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
        NULLIFY(solver1D)
        CALL Solvers_SolverGet(solvers,solver1DNavierStokesNumber,solver1D,err,error,*999)
        CALL ControlLoop_IterationNumberGet(controlLoop,iteration,err,error,*999)
        NULLIFY(parentLoop)
        CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
        NULLIFY(parentLoop2)
        CALL ControlLoop_ParentLoopGet(parentLoop,parentLoop2,err,error,*999)
        CALL ControlLoop_CurrentTimeInformationGet(parentLoop,startTime,stopTime,currentTime,timeIncrement, &
          & timestep,outputIteration,inputIteration,err,error,*999)
      ELSE
        localError="The solver number of "//TRIM(NumberToVString(solverNumber,"*",err,error))// &
         & " does not correspond with the Navier-Stokes solver number for 1D-0D fluid coupling."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    CALL ControlLoop_AbsoluteToleranceGet(controlLoop,couplingTolerance,err,error,*999)

    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver1D,solverEquations,err,error,*999)
    NULLIFY(boundaryConditions)
    CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(uDependentVariable)
    CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,uDependentVariable,err,error,*999)
    NULLIFY(u1DependentVariable)
    CALL Field_VariableGet(dependentField,FIELD_U1_VARIABLE_TYPE,u1DependentVariable,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(vMaterialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    NULLIFY(independentVariable)
    CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(dynamicVariable)
    CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)

    !Get the number of local nodes
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfLocalNodes1D,err,error,*999)
    numberOfLocalNodes1D=domainNodes%numberOfNodes
 
    boundaryNumber = 0
    boundaryConverged = .TRUE.
    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes1D
      !Check for the boundary node
      CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)

      !Get node characteristic wave direction (specifies inlet/outlet)
      coupled1D0DBoundary = .FALSE.
      DO componentIdx=1,2
        CALL FieldVariable_ParameterSetGetLocalNode(independentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,componentIdx,normalWave(componentIdx),err,error,*999)

        ! Check that this is a 1D-0D boundary
        CALL FieldVariable_LocalNodeDOFGet(uDependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
          & dependentDOF,err,error,*999)
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,uDependentVariable,boundaryConditionsVariable,err,error,*999)
        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dependentDOF,boundaryConditionType, &
          & err,error,*999)
        IF(boundaryConditionType == BOUNDARY_CONDITION_FIXED_CELLML) coupled1D0DBoundary = .TRUE.
      ENDDO !nodeIdx

      !!!-- F i n d   B o u n d a r y   N o d e s --!!!
      IF(ABS(normalWave(1))>ZERO_TOLERANCE.AND.boundaryNode.AND.coupled1D0DBoundary) THEN

        boundaryNumber = boundaryNumber + 1
        boundaryConverged(boundaryNumber) = .FALSE.
        !Get material parameters
        CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,1,a0Param,err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,2,eParam,err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,3,hParam,err,error,*999)
        beta=(4.0_DP*SQRT(PI)*eParam*hParam)/(3.0_DP*a0Param)

        ! C u r r e n t   Q 1 D , A 1 D , p C e l l M L   V a l u e s
        ! ------------------------------------------------------------
        !Get Q1D
        CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,1,q1d,err,error,*999)
        !Get A1D
        CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,2,a1d,err,error,*999)
        !Get pCellML
        CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,1,pCellML,err,error,*999)

        ! C h e c k  1 D / 0 D   C o n v e r g e n c e   f o r   t h i s   n o d e
        ! -------------------------------------------------------------------------
        IF(iteration == 1 .AND. timestep == 0) THEN
          ! Create the previous iteration field values type on the dependent field if it does not exist
          CALL FieldVariable_ParameterSetEnsureCreated(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetEnsureCreated(u1DependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
        ELSE
          !Get previous Q1D
          CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeIdx,1,qPrevious,err,error,*999)
          !Get previous A1D
          CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeIdx,2,aPrevious,err,error,*999)
          !Get previous pCellML value
          CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
            & versionIdx,derivativeIdx,nodeIdx,1,pPrevious,err,error,*999)
          ! Check if the boundary interface values have converged
          qError = ABS(qPrevious - q1d)
          aError = ABS(aPrevious - a1d)
          IF( qError < couplingTolerance .AND. aError < couplingTolerance) THEN
            boundaryConverged(boundaryNumber) = .TRUE.
          ENDIF
        ENDIF

        ! store current Q and p Boundary values as previous iteration value
        CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,1,q1d,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,2,a1d,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalNode(u1DependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,1,pCellML,err,error,*999)

      ENDIF !Find boundary nodes
    ENDDO !Loop over nodes
    CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateStart(u1DependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(u1DependentVariable,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, err,error,*999)
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
      ENDIF
      ! Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
      CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
      CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
      IF(numberOfGroupComputationNodes>1) THEN !use mpi
        !allocate array for mpi communication
        ALLOCATE(globalConverged(numberOfGroupComputationNodes),STAT=ERR) 
        IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",err,error,*999)
        CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
        IF(ALL(globalConverged)) THEN
          !CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
        ENDIF
        DEALLOCATE(globalConverged)
      ELSE
        IF(localConverged) THEN
          !CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"1D/0D coupling converged; # iterations: ", &
          !  & iteration,err,error,*999)
          CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
        ENDIF
      ENDIF

      ! If the solution hasn't converged, need to revert field values to pre-solve state
      ! before continued iteration. This will counteract the field updates that occur
      ! in Solver_DynamicMeanPredictedCalculate. Ignore for initialisation
      IF(timestep == 0) CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
      CALL ControlLoop_ContinueLoopGet(controlLoop,continueLoop,err,error,*999)
      IF(continueLoop) THEN
        CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP, &
          & err,error,*999)
        CALL FieldVariable_ParameterSetsCopy(dynamicVariable,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP, &
          & err,error,*999)
      ENDIF
    ENDIF

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
    TYPE(ControlLoopType), POINTER :: controlLoop
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: boundaryConditionType,boundaryNumber,boundaryType1D,componentIdx,computationNode,derivativeIdx,globalDof, &
      & groupCommunicator,inputIteration,iteration,localDof,maxIterations,mpiIError,nodeIdx,nodeNumber,numberOfLocalNodes1D, &
      & numberOfBoundaries,numberOfGroupComputationNodes,outputIteration,pSpecification(3),solver1dNavierStokesNumber, &
      & solver3dNavierStokesNumber,timestep,userNodeNumber,versionIdx
    REAL(DP) :: absoluteCouplingTolerance,absoluteCouplingTolerance2,currentTime,flow1D,flow3D,flowError,flowTolerance, &
      & maxFlowError,maxStressError,normalWave(2),relativeCouplingTolerance,startTime,stopTime,stress1D,stress3D,stressError, &
      & stressTolerance,timeIncrement
    LOGICAL :: boundaryNode,boundaryConverged(30),localConverged,continueLoop
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: subLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: equationsSet1D,equationsSet3D
    TYPE(EquationsType), POINTER :: equations1D,equations3D
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping1D,vectorMapping3D
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping1D,dynamicMapping3D
    TYPE(EquationsVectorType), POINTER :: vectorEquations1D,vectorEquations3D
    TYPE(FieldType), POINTER :: dependentField1D,dependentField3D,independentField,independentField1D
    TYPE(FieldVariableType), POINTER :: dependentVariable,uDependentVariable1D,u1DependentVariable1D,u2DependentVariable1D, &
      & dynamicVariable1D,dynamicVariable3D,independentVariable1D
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver1D,solver3D
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("NavierStokes_Couple3D1D",err,error,*999)

    continueLoop=.TRUE.
    
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    !Get solvers based on the problem type
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      CALL ControlLoop_CurrentWhileInformationGet(controlLoop,iteration,maxIterations,absoluteCouplingTolerance, &
        & relativeCouplingTolerance,continueLoop,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,timestep, &
        & outputIteration,inputIteration,err,error,*999)
      absoluteCouplingTolerance2=absoluteCouplingTolerance
      ! 1D solver & equations pointers
      solver1dNavierStokesNumber=2
      versionIdx=1
      derivativeIdx=1
      !TODO: make this more general!
      NULLIFY(subLoop)
      CALL ControlLoop_SubLoopGet(controlLoop,1,subLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
      NULLIFY(solver1D)
      CALL Solvers_SolverGet(solvers,solver1DNavierStokesNumber,solver1D,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver1D,solverEquations,err,error,*999)
      NULLIFY(boundaryConditions)
      CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet1D)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet1D,err,error,*999)
      NULLIFY(dependentField1D)
      CALL EquationsSet_DependentFieldGet(equationsSet1D,dependentField1D,err,error,*999)
      NULLIFY(uDependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U_VARIABLE_TYPE,uDependentVariable1D,err,error,*999)
      NULLIFY(u1DependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U1_VARIABLE_TYPE,u1DependentVariable1D,err,error,*999)
      NULLIFY(u2DependentVariable1D)
      CALL Field_VariableGet(dependentField1D,FIELD_U2_VARIABLE_TYPE,u2DependentVariable1D,err,error,*999)
      NULLIFY(independentField1D)
      CALL EquationsSet_IndependentFieldGet(equationsSet1D,independentField1D,err,error,*999)
      NULLIFY(independentVariable1D)
      CALL Field_VariableGet(independentField1D,FIELD_U_VARIABLE_TYPE,independentVariable1D,err,error,*999)
      ! 3D equations pointers
      solver3dNavierStokesNumber = 1
      ! TODO: make this more general!
      NULLIFY(subLoop)
      CALL ControlLoop_SubLoopGet(controlLoop,2,subLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
      CALL Solvers_SolverGet(solvers,solver3DNavierStokesNumber,solver3D,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver3D,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet3D)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet3D,err,error,*999)
      NULLIFY(dependentField3D)
      CALL EquationsSet_DependentFieldGet(equationsSet3D,dependentField3D,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for 3D-1D Navier-Stokes fluid coupling."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations1D)
    CALL EquationsSet_EquationsGet(equationsSet1D,equations1D,err,error,*999)
    NULLIFY(vectorEquations1D)
    CALL Equations_VectorEquationsGet(equations1D,vectorEquations1D,err,error,*999)
    NULLIFY(vectorMapping1D)
    CALL EquationsVector_VectorMappingGet(vectorEquations1D,vectorMapping1D,err,error,*999)
    NULLIFY(dynamicMapping1D)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping1D,dynamicMapping1D,err,error,*999)
    NULLIFY(dynamicVariable1D)
    CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping1D,dynamicVariable1D,err,error,*999)
    NULLIFY(equations3D)
    CALL EquationsSet_EquationsGet(equationsSet3D,equations3D,err,error,*999)
    NULLIFY(vectorEquations3D)
    CALL Equations_VectorEquationsGet(equations3D,vectorEquations3D,err,error,*999)
    NULLIFY(vectorMapping3D)
    CALL EquationsVector_VectorMappingGet(vectorEquations3D,vectorMapping3D,err,error,*999)
    NULLIFY(dynamicMapping3D)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping3D,dynamicMapping3D,err,error,*999)
    NULLIFY(dynamicVariable3D)
    CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping3D,dynamicVariable3D,err,error,*999)

    !Get the number of local nodes
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField1D,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfLocalNodes1D,err,error,*999)

    boundaryNumber = 0
    maxStressError = 0.0_DP
    maxFlowError = 0.0_DP
    boundaryConverged = .TRUE.
    !!!--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    DO nodeIdx=1,numberOfLocalNodes1D
      CALL DomainNodes_NodeUserNumberGet(domainNodes,nodeIdx,userNodeNumber,err,error,*999)
      
      !Check for the boundary node- go to next if not boundary
      CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
      IF(.NOT. boundaryNode) CYCLE

      !Get node characteristic wave direction (specifies inlet/outlet)
      DO componentIdx=1,2
        CALL FieldVariable_ParameterSetGetLocalNode(independentVariable1D,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,componentIdx,normalWave(componentIdx),err,error,*999)
      ENDDO !componentIdx

      !!!-- F i n d   B o u n d a r y   N o d e s --!!!
      IF(ABS(normalWave(1))>ZERO_TOLERANCE .OR. ABS(normalWave(2))>ZERO_TOLERANCE) THEN
        ! Check that this is a coupled 3D-1D boundary
        boundaryType1D = 0
        NULLIFY(boundaryConditionsVariable)
        CALL BoundaryConditions_VariableGet(boundaryConditions,uDependentVariable1D,boundaryConditionsVariable,err,error,*999)
        CALL FieldVariable_ComponentDOFGetUserNode(uDependentVariable1D,versionIdx,derivativeIdx,userNodeNumber,1,localDOF, &
          & globalDOF,err,error,*999)
        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDOF,boundaryConditionType, &
          & err,error,*999)
        IF(boundaryConditionType==BOUNDARY_CONDITION_COUPLING_FLOW) boundaryType1D = BOUNDARY_CONDITION_COUPLING_FLOW
        CALL FieldVariable_ComponentDOFGetUserNode(uDependentVariable1D,versionIdx,derivativeIdx,userNodeNumber,2,localDOF, &
          & globalDOF,err,error,*999)
        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDOF,boundaryConditionType, &
          & err,error,*999)
        IF(boundaryConditionType==BOUNDARY_CONDITION_COUPLING_STRESS) THEN
          IF(boundaryType1D==BOUNDARY_CONDITION_COUPLING_FLOW) THEN
            localError="Boundary type for node number "//TRIM(NumberToVString(userNodeNumber,"*",err,error))// &
              & " is set as both FLOW and STRESS for 3D-1D Navier-Stokes fluid coupling."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          boundaryType1D = BOUNDARY_CONDITION_COUPLING_STRESS
        ENDIF
        ! If this is not a coupled 3D-1D boundary, go to the next node
        IF(boundaryType1D==0) CYCLE

        boundaryNumber = boundaryNumber + 1
        boundaryConverged(boundaryNumber) = .FALSE.
        CALL FieldVariable_ParameterSetEnsureCreated(u1DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
        !Get 1D flow and stress from current iteration
        CALL FieldVariable_ParameterSetGetLocalNode(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,2,stress1D,err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalNode(u2DependentVariable1D,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,3,flow1D,err,error,*999)
        !Get 3D flow and stress
        CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeIdx,1,flow3D,err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable1D,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
          & nodeNumber,2,stress3D,err,error,*999)
        IF(diagnostics1) THEN
          IF(boundaryNumber ==1) CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"------3D-1D----- iteration:  ",iteration, &
            & err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  node:  ",userNodeNumber,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    1D flow:  ",flow1D,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    3D flow:  ",flow3D,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    1D stress: ",stress1D,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    3D stress: ",stress3D,err,error,*999)
        ENDIF

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
        ENDIF
        maxFlowError = MAX(ABS(maxFlowError),ABS(flowError))
        maxStressError = MAX(ABS(maxStressError),ABS(stressError))

        ! Set current iteration values to previous
        CALL FieldVariable_ParameterSetUpdateLocalNode(u1DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,1,flow3D,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalNode(u1DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,2,stress3D,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalNode(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,2,stress1D,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateLocalNode(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
          & versionIdx,derivativeIdx,nodeIdx,3,flow1D,err,error,*999)

      ENDIF !Find boundary nodes
    END DO !Loop over nodes
    CALL FieldVariable_ParameterSetUpdateStart(u1DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateStart(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(u1DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(u2DependentVariable1D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE,err,error,*999)
    numberOfBoundaries = boundaryNumber

    ! ------------------------------------------------------------------
    ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
    ! ------------------------------------------------------------------
    ! Check whether all boundaries on the local process have converged
    localConverged=.FALSE.
    globalConverged=.FALSE.
    IF(numberOfBoundaries==0.OR.ALL(boundaryConverged(1:numberOfBoundaries))) THEN
      localConverged = .TRUE.
    ENDIF
    !Need to check that boundaries have converged globally (on all domains) if this is a MPI problem
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,computationNode,err,error,*999)
    IF(numberOfGroupComputationNodes>1) THEN !use mpi
      !allocate array for mpi communication
      ALLOCATE(globalConverged(numberOfGroupComputationNodes),STAT=err)
      IF(err/=0) CALL FlagError("Could not allocate global convergence check array.",err,error,*999)
      CALL MPI_ALLREDUCE(localConverged,globalConverged,1,MPI_LOGICAL,MPI_LAND,groupCommunicator,mpiIError)
      CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
      IF(ALL(globalConverged)) CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
      IF(diagnostics1) THEN
        IF(ALL(globalConverged)) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D coupling converged; # iterations: ",iteration,err,error,*999)
        ELSE
          CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"Rank ",computationNode," 3D/1D max flow error:  ",maxFlowError, &
            & err,error,*999)
          CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"Rank ",computationNode," 3D/1D max stress error:  ",maxStressError, &
            & err,error,*999)
        ENDIF
      ENDIF
      DEALLOCATE(globalConverged)
    ELSE
      IF(localConverged) CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
      IF(diagnostics1) THEN
        IF(localConverged) THEN
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D coupling converged; # iterations: ",iteration,err,error,*999)        
        ELSE
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D max flow error:  ",maxFlowError,err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"3D/1D max stress error: ",maxStressError,err,error,*999)
        ENDIF
      ENDIF
    ENDIF

    !If the solution hasn't converged, need to revert field values to pre-solve state
    !before continued iteration. This will counteract the field updates that occur
    !in Solver_DynamicMeanPredictedCalculate. Ignore for initialisation
    IF(timestep == 0) CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
    CALL ControlLoop_ContinueLoopGet(controlLoop,continueLoop,err,error,*999)
    IF(continueLoop) THEN
      !Reset 1D values
      CALL FieldVariable_ParameterSetsCopy(dynamicVariable1D,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP, &
        & err,error,*999)
      CALL FieldVariable_ParameterSetsCopy(dynamicVariable1D,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP, &
        & err,error,*999)
      !Reset 3D values
      CALL FieldVariable_ParameterSetsCopy(dynamicVariable3D,FIELD_PREVIOUS_VALUES_SET_TYPE,FIELD_VALUES_SET_TYPE,1.0_DP, &
        & err,error,*999)
      CALL FieldVariable_ParameterSetsCopy(dynamicVariable3D,FIELD_PREVIOUS_RESIDUAL_SET_TYPE,FIELD_RESIDUAL_SET_TYPE,1.0_DP, &
        & err,error,*999)
    ENDIF

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
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SolverType), POINTER :: solver
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: branchNumber,componentIdx,derivativeIdx,groupCommunicator,i,inputIteration,iteration,mpiIError,nodeIdx, &
      & nodeNumber,numberOfBranches,numberOfGroupComputationNodes,numberOfNodes,numberOfVersions,outputIteration, &
      & pSpecification(3),solver1dNavierStokesNumber,solverNumber,timestep,versionIdx
    REAL(DP) :: a0Param,aCharacteristic(7),alpha,aNavierStokes(7),aNew,beta,couplingTolerance,currentTime,eParam,hParam, &
      & l2ErrorA(100),l2ErrorQ(100),l2ErrorW(30),normalWave,penaltyCoeff,startTime,stopTime,qCharacteristic(7),qNavierStokes(7), &
      & rho,timeIncrement,totalErrorA,totalErrorMass,totalErrorMomentum,totalErrorQ,totalErrorW,totalErrorWPrevious, &
      & wCharacteristic(2,7),wError(2,7),wNavierStokes(2,7),wNext(2,7),wPrevious(2,7)
    LOGICAL :: branchConverged(100),localConverged,MPI_LOGICAL,boundaryNode,fluxDiverged
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(ControlLoopType), POINTER :: parentLoop,subLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,independentField,materialsField
    TYPE(FieldVariableType), POINTER :: uDependentVariable,vDependentVariable,uEquationsSetVariable,uIndependentVariable, &
      & uMaterialsVariable,vMaterialsVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolverType), POINTER :: solver1DNavierStokes
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("NavierStokes_CoupleCharacteristics",err,error,*999)

    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
      NULLIFY(solver1DNavierStokes)
      CALL Solvers_SolverGet(solvers,solver1dNavierStokesNumber,solver1DNavierStokes,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,startTime,stopTime,currentTime,timeIncrement, &
        & timestep,outputIteration,inputIteration,err,error,*999)
      CALL ControlLoop_IterationNumberGet(controlLoop,iteration,err,error,*999)
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
       & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      solver1dNavierStokesNumber=2
      NULLIFY(parentLoop)
      CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
      NULLIFY(subLoop)
      CALL ControlLoop_SubLoopGet(parentLoop,2,subLoop,err,error,*999)
      NULLIFY(solvers)
      CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
      NULLIFY(solver1DNavierStokes)
      CALL Solvers_SolverGet(solvers,solver1DNavierStokesNumber,solver1DNavierStokes,err,error,*999)
      CALL ControlLoop_IterationNumberGet(controlLoop,iteration,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(parentLoop,startTime,stopTime,currentTime,timeIncrement,timestep, &
        & outputIteration,inputIteration,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for 1D-0D Navier-Stokes fluid coupling."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
    CALL ControlLoop_AbsoluteToleranceGet(controlLoop,couplingTolerance,err,error,*999)
    NULLIFY(solverEquations)
    CALL Solver_SolverEquationsGet(solver1DNavierStokes,solverEquations,err,error,*999)
    NULLIFY(solverMapping)
    CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
    NULLIFY(equationsSet)
    CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(uDependentVariable)
    CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,uDependentVariable,err,error,*999)
    NULLIFY(vDependentVariable)
    CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,vDependentVariable,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(uMaterialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
    NULLIFY(vMaterialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    NULLIFY(uIndependentVariable)
    CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,uIndependentVariable,err,error,*999)
    NULLIFY(equationsSetField)
    CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
    NULLIFY(uEquationsSetVariable)
    CALL Field_VariableGet(equationsSetField,FIELD_U_VARIABLE_TYPE,uEquationsSetVariable,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainTopology)
    CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
    NULLIFY(domainNodes)
    CALL DomainTopology_DomainNodesget(domainTopology,domainNodes,err,error,*999)
    !Get the number of local nodes
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

    !Get material constants
    CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
    CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,3,alpha,err,error,*999)
    CALL FieldVariable_ParameterSetGetConstant(uEquationsSetVariable,FIELD_VALUES_SET_TYPE,1,penaltyCoeff,err,error,*999)

    !--  L o o p   O v e r   L o c a l  N o d e s  --!!!
    CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
    DO nodeIdx=1,numberOfNodes
      CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
      derivativeIdx = 1
      CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)

      !DEBUG
      CALL FieldVariable_ParameterSetGetLocalNode(uIndependentVariable,FIELD_VALUES_SET_TYPE,1,1,nodeIdx,1,normalWave, &
        & err,error,*999)
      !Find branch nodes
      IF(numberOfVersions>1) THEN
        branchNumber = branchNumber + 1
        branchConverged(branchNumber) = .FALSE.

        wError = 0.0_DP
        i = 0
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            i = i + 1
            CALL FieldVariable_ParameterSetGetLocalNode(uIndependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
              & nodeIdx,componentIdx,normalWave,err,error,*999)
            IF(ABS(normalWave)>ZERO_TOLERANCE) THEN

              ! Get the previously set characteristic (W) for this timestep-
              !  if this is the first iteration it will be based on extrapolated values
              !  otherwise it will come from the last iteration of this subroutine.
              CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,componentIdx,wPrevious(componentIdx,versionIdx),err,error,*999)

              !Get material parameters
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,a0Param,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,eParam,err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,3,hParam,err,error,*999)
              beta=(4.0_DP*SQRT(PI)*eParam*hParam)/(3.0_DP*a0Param)

              ! Get current Q,A values based on N-S solve
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,1,qNavierStokes(versionIdx),err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,2,aNavierStokes(versionIdx),err,error,*999)

              ! Calculate the characteristic based on the values converged upon by the
              !  N-S solver at this iteration.
              wNavierStokes(componentIdx,versionIdx)= ((qNavierStokes(versionIdx)/aNavierStokes(versionIdx))+ &
               & normalWave*4.0_DP*SQRT(beta/(2.0_DP*rho))*(aNavierStokes(versionIdx)**(0.25_DP) - (a0Param)**(0.25_DP)))

              IF(boundaryNode) THEN
                aNew = (1.0_DP/(beta/(2.0_DP*rho)))**2.0_DP*((wNavierStokes(componentIdx,versionIdx))/8.0_DP+ &
                 & SQRT(beta/(2.0_DP*rho))*((a0Param)**0.25_DP))**4.0_DP
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
                  & versionIdx,derivativeIdx,nodeIdx,2,aNew,err,error,*999)
              ENDIF

              ! Get characteristic (flux conserving) Q,A values
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,versionIdx, &
                & derivativeIdx,nodeIdx,1,qCharacteristic(versionIdx),err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_UPWIND_VALUES_SET_TYPE,versionIdx, &
                & derivativeIdx,nodeIdx,2,aCharacteristic(versionIdx),err,error,*999)

              ! Calculate the characteristic based on the upwind values
              wCharacteristic(componentIdx,versionIdx)= ((qCharacteristic(versionIdx)/aCharacteristic(versionIdx))+ &
               & normalWave*4.0_DP*SQRT((beta/(2.0_DP*rho)))*(aCharacteristic(versionIdx)**(0.25_DP) - (a0Param)**(0.25_DP)))
            ENDIF
          ENDDO !versionIdx
        ENDDO !componentIdx

        ! Evaluate error between current and previous Q,A values
        IF(numberOfVersions > 1 ) THEN
          CALL L2Norm(qNavierStokes-qCharacteristic,l2ErrorQ(branchNumber),err,error,*999)
          CALL L2Norm(aNavierStokes-aCharacteristic,l2ErrorA(branchNumber),err,error,*999)
        ENDIF
        ! Check if the branch values have converged
        IF((ABS(l2ErrorQ(branchNumber)) < couplingTolerance).AND.(ABS(l2ErrorA(branchNumber)) < couplingTolerance)) THEN
          branchConverged(branchNumber) = .TRUE.
        ENDIF
        totalErrorQ = totalErrorQ + l2ErrorQ(branchNumber)
        totalErrorA = totalErrorA + l2ErrorA(branchNumber)

        wNext = ((wNavierStokes + wCharacteristic)/2.0_DP)
        ! If N-S/C w values did not converge re-solve with new w.
        IF(numberOfVersions > 1) THEN
          DO componentIdx=1,2
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_ParameterSetGetLocalNode(uIndependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                & nodeIdx,componentIdx,normalWave,err,error,*999)
              IF(ABS(normalWave)>ZERO_TOLERANCE) THEN
                !Update W value
                CALL FieldVariable_ParameterSetUpdateLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,componentIdx,wNext(componentIdx,versionIdx),err,error,*999)
              ENDIF
            ENDDO !versionIdx
          ENDDO !componentIdx
        ENDIF

      ENDIF !Find boundary nodes
    ENDDO !Loop over nodes
    CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateStart(vDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE,err,error,*999)
    CALL FieldVariable_ParameterSetUpdateFinish(vDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    numberOfBranches = branchNumber

    ! ------------------------------------------------------------------
    ! C h e c k   G l o b a l   C o u p l i n g   C o n v e r g e n c e
    ! ------------------------------------------------------------------
    ! Check whether all branches on the local process have converged
    IF(numberOfBranches == 0.OR.ALL(branchConverged(1:numberOfBranches))) THEN
      localConverged = .TRUE.
    ELSE
      localConverged = .FALSE.
    ENDIF
    ! Need to check that boundaries have converged globally (on all domains) if this is a parallel problem
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
    CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
    IF(numberOfGroupComputationNodes>1) THEN !use mpi
      !allocate array for mpi communication
      ALLOCATE(globalConverged(numberOfGroupComputationNodes),STAT=ERR) 
      IF(ERR/=0) CALL FlagError("Could not allocate global convergence check array.",err,error,*999)
      CALL MPI_ALLGATHER(localConverged,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL,groupCommunicator,mpiIError)
      CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
      IF(ALL(globalConverged)) THEN
        !CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
       ENDIF
      DEALLOCATE(globalConverged)
    ELSE
      IF(localConverged) THEN
        !CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Navier-Stokes/Characteristic converged; # iterations: ", &
        !  & iteration,err,error,*999)
        CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
      ENDIF
    ENDIF

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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,decompositionLocalElementNumber,esSpecification(3),gaussIdx,i,j,localElementNumber, &
      & numberOfDimensions,numberOfGauss,numberOfXi,outputType,startElement,stopElement,userElementNumber,variableType,xiIdx
    REAL(DP) :: dUdXi(3,3),dXidX(3,3),dUdX(3,3),dUdXTrans(3,3),D(3,3),gaussWeight,rateOfStrainMag,shearRateDefault, &
      & shearRate,strainRate,velocityGauss(3)
    LOGICAL :: calculateFaces,calculateLines,defaultUpdate,elementExists,ghostElement
    TYPE(BasisType), POINTER :: dependentBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionElementsType), POINTER :: decompositionElements
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters
    TYPE(FieldVariableType), POINTER :: dependentVariable,vMaterialsVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_ShearRateCalculate",err,error,*999)

    CALL EquationsSet_OutputTypeGet(equationsSet,outputType,err,error,*999)
    IF(outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) &
      & CALL WriteString(GENERAL_OUTPUT_TYPE,"...Calculating shear rate...",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(vMaterialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
      NULLIFY(dependentVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(dependentVariable,variableType,err,error,*999)
      !Get the mesh decomposition and mapping
      NULLIFY(dependentField)
      CALL FieldVariable_FieldGet(dependentVariable,dependentField,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
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
      NULLIFY(domainMappings)
      CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
      NULLIFY(elementsMapping)
      CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
      CALL DomainMapping_InternalStartGet(elementsMapping,startElement,err,error,*999)
      CALL DomainMapping_BoundaryFinishGet(elementsMapping,stopElement,err,error,*999)

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,variableType,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,variableType,geometricInterpPoint,err,error,*999)
      
      defaultUpdate=.FALSE.

      ! Loop over internal and boundary elements, skipping ghosts
      DO elementIdx=startElement,stopElement
        CALL DomainMapping_NumberGet(elementsMapping,elementIdx,localElementNumber,err,error,*999)
        CALL DecompositionElements_ElementUserNumberGet(decompositionElements,localElementNumber,userElementNumber,err,error,*999)
        !Check computation node for elementIdx
        elementExists=.FALSE.
        ghostElement=.TRUE.
        CALL DecompositionElements_ElementCheckExists(decompositionElements,userElementNumber,elementExists, &
          & decompositionLocalElementNumber,ghostElement,err,error,*999)
        IF(diagnostics1) THEN
          IF(ghostElement) CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Ghost: ",userElementNumber,err,error,*999)
        ENDIF

        IF(elementExists) THEN
          NULLIFY(dependentBasis)
          CALL DomainElements_ElementBasisGet(domainElements,localElementNumber,dependentBasis,err,error,*999)
          NULLIFY(quadratureScheme)
          CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)

          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber,dependentInterpParameters, &
            & err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,localElementNumber,geometricInterpParameters, &
            & err,error,*999)

          !Loop over gauss points
          DO gaussIdx=1,numberOfGauss
            !Get interpolated velocity
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & dependentInterpPoint,err,error,*999)
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx, &
              & geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_VOLUME_TYPE,geometricInterpPointMetrics, &
              & err,error,*999)

            CALL BasisQuadratureScheme_GaussWeightGet(quadratureScheme,gaussIdx,gaussWeight,err,error,*999)
            
            !Interpolated values at gauss point
            dXidX=0.0_DP
            dUdXi=0.0_DP
            velocityGauss=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
            DO xiIdx=1,numberOfXi              
              dXidX(xiIdx,1:numberOfDimensions)=geometricInterpPointMetrics%dXidX(xiIdx,1:numberOfDimensions)
              dUdXi(1:numberOfDimensions,xiIdx)=dependentInterpPoint%values(1:numberOfDimensions, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
            ENDDO !xiIdx
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
            CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(vMaterialsVariable,FIELD_VALUES_SET_TYPE,gaussIdx, &
              & localElementNumber,2,shearRate,err,error,*999)

          ENDDO !gaussIdx
        ENDIF ! check for ghost element
        IF(defaultUpdate) EXIT
      ENDDO !elementIdx
      CALL FieldVariable_ParameterSetUpdateStart(vMaterialsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(vMaterialsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)

      IF(defaultUpdate) THEN
        shearRateDefault=1.0E-10_DP
        CALL FieldVariable_ComponentValuesInitialise(vMaterialsVariable,FIELD_VALUES_SET_TYPE,1,shearRateDefault,err,error,*999)
        IF(diagnostics1) CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Setting default shear field values...", &
          & shearRateDefault,err,error,*999)
      ENDIF

    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    LOGICAL :: convergedFlag
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("NavierStokes_FiniteElementPreResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
     
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
      ! Shear rate should either be calculated here to update at each minor iteration
      ! or during post solve so it is updated once per timestep
      !CALL NavierStokes_ShearRateCalculate(equationsSet,err,error,*999)
    CASE(EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE)
!!TODO: what is setting the converged flag doing???
      convergedFlag = .FALSE.
      CALL NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*999)
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_OPTIMISED_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
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
        & TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

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
  SUBROUTINE NavierStokes_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,equationsSetIdx2,esSpecification(3),esSpecification2(3),iteration3D1D,loopType,loopLevel, &
      & numberOfEquationsSets,numberOfEquationsSets2,numberOfSolvers,numberOfSolvers2,numberOfSubLoops,numberOfSubLoops2, &
      & numberOfSubLoops3,numberOfSubLoops4,pSpecification(3),solverIdx,solverIdx2,solveType,subloopIdx,subloopIdx2,subloopIdx3
    REAL(DP) :: absolute3D0DTolerance,relative3D0DTolerance
    LOGICAL :: calculateFaces,calculateLines,continueLoop,convergedFlag
    CHARACTER(70) :: label
    TYPE(ControlLoopType), POINTER :: subloop,subloop2,subloop3,iterativeWhileLoop2,iterativeWhileLoop3,parentLoop
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(EquationsSetType), POINTER :: equationsSet,equationsSet2,coupledEquationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: navierStokesSolver,navierStokesSolver3D,navierStokesSolver1D,solver,solver2
    TYPE(SolversType), POINTER :: solvers,solvers2
    TYPE(SolverEquationsType), POINTER :: solverEquations,solverEquations2,solverEquations1D,solverEquations3D
    TYPE(SolverMappingType), POINTER :: solverMapping,solverMapping1D,solverMapping2,solverMapping3D
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_PostLoop",err,error,*999)

    convergedFlag = .FALSE.
    absolute3D0DTolerance = 0.0_DP
    relative3D0DTolerance = 0.0_DP
    continueLoop = .TRUE.
    
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_NAVIER_STOKES_SUBTYPE, &
      &  PROBLEM_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      &  PROBLEM_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
      &  PROBLEM_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      &  PROBLEM_ALE_NAVIER_STOKES_SUBTYPE)
      !Do nothing
    CASE(PROBLEM_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_CONSTITUTIVE_RBS_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(loopType)
      CASE(CONTROL_TIME_LOOP_TYPE)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(navierStokesSolver)
        CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(navierStokesSolver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        CALL Decomposition_CalculateFacesGet(decomposition,calculateFaces,err,error,*999)
        IF(calculateFaces) CALL NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*999)
        CALL NavierStokes_PostSolveOutputData(navierStokesSolver,err,error,*999)
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for the specified Navier-Stokes problem type."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
      &  PROBLEM_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        !Do nothing
      CASE(CONTROL_TIME_LOOP_TYPE)
        !Global time loop - export data
        NULLIFY(subLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,subLoop,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
        NULLIFY(navierStokesSolver)
        CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
        CALL NavierStokes_PostSolveOutputData(navierStokesSolver,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(navierStokesSolver)
        CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
        CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,err,error,*999)
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for a Coupled 1D0D Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_STREE1D0D_ADV_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        !CellML simple loop - do nothing
      CASE(CONTROL_TIME_LOOP_TYPE)
        !Global time loop - export data
        NULLIFY(subLoop)
        CALL ControlLoop_SubLoopGet(controlLoop,1,subLoop,err,error,*999)
        NULLIFY(subLoop2)
        CALL ControlLoop_SubLoopGet(subLoop,2,subLoop2,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(subLoop2,solvers,err,error,*999)
        NULLIFY(navierStokesSolver)
        CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
        CALL NavierStokes_PostSolveOutputData(navierStokesSolver,err,error,*999)
      CASE(CONTROL_WHILE_LOOP_TYPE)
        !Couple 1D/0D loop
        CALL ControlLoop_LoopLevelGet(controlLoop,loopLevel,err,error,*999)
        IF(loopLevel==2) THEN
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,subLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
          NULLIFY(navierStokesSolver)
          CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
          !update 1D/0D coupling parameters and check convergence
          CALL NavierStokes_Couple1D0D(controlLoop,navierStokesSolver,err,error,*999)
          !Couple Navier-Stokes/Characteristics loop
        ELSE IF(loopLevel==3) THEN
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
          NULLIFY(navierStokesSolver)
          CALL Solvers_SolverGet(solvers,2,navierStokesSolver,err,error,*999)
          CALL NavierStokes_CoupleCharacteristics(controlLoop,navierStokesSolver,err,error,*999)
        ELSE
          localError="The while loop level of "//TRIM(NumberToVString(loopLevel,"*",err,error))// &
            & " is invalid for a Coupled 1D0D Navier-Stokes problem."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for a Coupled 1D0D Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE, &
      & PROBLEM_COUPLED3D0D_NAVIER_STOKES_SUBTYPE)
      SELECT CASE(loopType)
        !Simple loops- could be 3D or 0D
      CASE(CONTROL_SIMPLE_TYPE)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
        DO solverIdx=1,numberOfSolvers
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
          CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
          SELECT CASE(solveType)
          CASE(SOLVER_DYNAMIC_TYPE)
            NULLIFY(solverEquations)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            NULLIFY(solverMapping)
            CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
            CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
            DO equationsSetIdx = 1,numberOfEquationsSets
              NULLIFY(equationsSet)
              CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
              CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
              SELECT CASE(esSpecification(3))
              CASE(EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
                & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE)
                !--- 3 D   T r a n s i e n t   N a v i e r - S t o k e s   E q u a t i o n s---
                NULLIFY(coupledEquationsSet)
                !TODO: This is a bit of a hack- make more robust
                !3D-1D iteration
                NULLIFY(parentLoop)
                CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
                CALL ControlLoop_IterationNumberGet(parentLoop,iteration3D1D,err,error,*999)
                !1D-0D
                IF(pSpecification(3)==PROBLEM_MULTISCALE_NAVIER_STOKES_SUBTYPE) THEN
                  NULLIFY(iterativeWhileLoop2)
                  CALL ControlLoop_SubLoopGet(parentLoop,1,iterativeWhileLoop2,err,error,*999)
                  !1D
                  NULLIFY(iterativeWhileLoop3)
                  CALL ControlLoop_SubLoopGet(iterativeWhileLoop2,2,iterativeWhileLoop3,err,error,*999)
                  CALL ControlLoop_NumberOfSubLoopsGet(iterativeWhileLoop3,numberOfSubLoops,err,error,*999)
                  IF(numberOfSubLoops/=0) CALL FlagError("Unrecognized subloop pattern for 3D-1D-0D!.",err,error,*999)
                  NULLIFY(solvers2)
                  CALL ControlLoop_SolversGet(iterativeWhileLoop3,solvers2,err,error,*999)
                  CALL Solvers_NumberOfSolversGet(solvers2,numberOfSolvers2,err,error,*999)
                  DO solverIdx2=1,numberOfSolvers2
                    NULLIFY(solver2)
                    CALL Solvers_SolverGet(solvers2,solverIdx2,solver2,err,error,*999)
                    NULLIFY(solverEquations2)
                    CALL Solver_SolverEquationsGet(solver2,solverEquations2,err,error,*999)
                    NULLIFY(solverMapping2)
                    CALL SolverEquations_SolverMappingGet(solverEquations2,solverMapping2,err,error,*999)
                    CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping2,numberOfEquationsSets2,err,error,*999)
                    DO equationsSetIdx2 = 1,numberOfEquationsSets2
                      NULLIFY(equationsSet2)
                      CALL SolverMapping_EquationsSetGet(solverMapping2,equationsSetIdx2,equationsSet2,err,error,*999)
                      CALL EquationsSet_SpecificationGet(equationsSet2,3,esSpecification2,err,error,*999)
                      SELECT CASE(esSpecification2(3))
                      CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
                        &  EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE)
                        IF(ASSOCIATED(coupledEquationsSet)) &
                          & CALL FlagError("Coupled 3D-1D equations set already found for multiscale Navier-Stokes problem.", &
                          & err,error,*999)
                        coupledEquationsSet=>equationsSet2
                      CASE DEFAULT
                        ! Do nothing
                      END SELECT
                    ENDDO !equationsSetIdx2
                  ENDDO !solverIdx2
                ENDIF
                CALL ControlLoop_ContinueLoopSet(parentLoop,.FALSE.,err,error,*999)
              CASE DEFAULT
                localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                  & " is not valid for a dynamic solver in a simple loop for a multiscale Navier-Stokes problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !equationsSetIdx
          CASE(SOLVER_DAE_TYPE)
            !CellML solver simple loop - do nothing
          CASE DEFAULT
            localError="The solve type of "//TRIM(NumberToVString(solveType,"*",err,error))// &
              & " is invalid for a simple loop in a Navier-Stokes multiscale problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO ! solverIdx
      CASE(CONTROL_TIME_LOOP_TYPE)
        !Global time loop - export data from all dynamic solvers and equations sets
        !Check for subloops two layers down
        CALL ControlLoop_NumberOfSubLoopsGet(controlLoop,numberOfSubLoops,err,error,*999)
        IF(numberOfSubLoops==0) CALL FlagError("Navier-Stokes Multiscale problem time loop has no subloops.",err,error,*999)
        DO subloopIdx=1,numberOfSubLoops
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,subLoopIdx,subLoop,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsGet(subLoop,numberOfSubLoops2,err,error,*999)
          IF(numberOfSubLoops2==0) CALL FlagError("Navier-Stokes Multiscale problem level 2 loop has no subloops.",err,error,*999)
          DO subloopIdx2=1,numberOfSubLoops2
            NULLIFY(subLoop2)
            CALL ControlLoop_SubLoopGet(subLoop,subLoopIdx2,subLoop2,err,error,*999)
            CALL ControlLoop_NumberOfSubLoopsGet(subLoop2,numberOfSubLoops3,err,error,*999)
            IF(numberOfSubLoops3>0) THEN
              !1D output
              DO subloopIdx3=1,numberOfSubLoops3
                NULLIFY(subLoop3)
                CALL ControlLoop_SubLoopGet(subLoop2,subLoopIdx3,subLoop3,err,error,*999)
                CALL ControlLoop_NumberOfSubLoopsGet(subLoop3,numberOfSubLoops4,err,error,*999)
                IF(numberOfSubLoops4>0) CALL FlagError("Unrecognized number of subloops in 3D-1D-0D.",err,error,*999)
                NULLIFY(solvers)
                CALL ControlLoop_SolversGet(subLoop3,solvers,err,error,*999)
                CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
                DO solverIdx=1,numberOfSolvers
                  NULLIFY(solver)
                  CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
                  CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
                  IF(solveType==SOLVER_DYNAMIC_TYPE) CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
                ENDDO !solverIdx
              ENDDO !subLoopIdx3
            ELSE
              !3D output
              NULLIFY(solvers)
              CALL ControlLoop_SolversGet(subLoop2,solvers,err,error,*999)
              CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
              DO solverIdx=1,numberOfSolvers
                NULLIFY(solver)
                CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
                CALL Solver_SolverTypeGet(solver,solveType,err,error,*999)
                IF(solveType==SOLVER_DYNAMIC_TYPE) CALL NavierStokes_PostSolveOutputData(solver,err,error,*999)
              ENDDO !solverIdx
            ENDIF
          ENDDO ! subloop2
        ENDDO ! subloop 1
      CASE(CONTROL_WHILE_LOOP_TYPE)
        !TODO: here we get loop type by label- may need to think of a more robust way to do this
        CALL ControlLoop_LabelGet(controlLoop,label,err,error,*999)
        SELECT CASE(label)
        CASE("3D-0D Iterative Loop")
          !Will handle 3D-0D coupling in calc boundary flux
          !update 3D/D coupling parameters and check convergence
          !CALL NavierStokes_Couple3D0D(controlLoop,err,error,*999)
        CASE("3D-1D Iterative Loop")
          !update 1D/3D coupling parameters and check convergence
          CALL NavierStokes_Couple3D1D(controlLoop,err,error,*999)
        CASE("1D-0D Iterative Coupling Convergence Loop")
          NULLIFY(navierStokesSolver1D,navierStokesSolver3D)
          !TODO: make this more general!e
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoop,2,subLoop,err,error,*999)
          NULLIFY(solvers)
          CALL ControlLoop_SolversGet(subLoop,solvers,err,error,*999)
          NULLIFY(navierStokesSolver1D)
          CALL Solvers_SolverGet(solvers,2,navierStokesSolver1D,err,error,*999)
          !Update 1D/0D coupling parameters and check convergence
          CALL NavierStokes_Couple1D0D(controlLoop,navierStokesSolver1D,err,error,*999)
          CALL ControlLoop_ContinueLoopGet(controlLoop,continueLoop,err,error,*999)
          IF(.NOT.continueLoop) THEN
            !TODO: make this more general!
            NULLIFY(parentLoop)
            CALL ControlLoop_ParentLoopGet(controlLoop,parentLoop,err,error,*999)
            NULLIFY(subLoop2)
            CALL ControlLoop_SubLoopGet(parentLoop,2,subLoop2,err,error,*999)
            NULLIFY(solvers2)
            CALL ControlLoop_SolversGet(subLoop2,solvers2,err,error,*999)
            NULLIFY(navierStokesSolver3D)
            CALL Solvers_SolverGet(solvers2,1,navierStokesSolver3D,err,error,*999)
            CALL ControlLoop_IterationNumberGet(parentLoop,iteration3D1D,err,error,*999)
            NULLIFY(solverEquations1D)
            CALL Solver_SolverEquationsGet(navierStokesSolver1D,solverEquations1D,err,error,*999)
            NULLIFY(solverEquations3D)
            CALL Solver_SolverEquationsGet(navierStokesSolver3D,solverEquations3D,err,error,*999)
            NULLIFY(solverMapping1D)
            CALL SolverEquations_SolverMappingGet(solverEquations1D,solverMapping1D,err,error,*999)
            NULLIFY(solverMapping3D)
            CALL SolverEquations_SolverMappingGet(solverEquations3D,solverMapping3D,err,error,*999)
            NULLIFY(equationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping1D,1,equationsSet,err,error,*999)
            NULLIFY(coupledEquationsSet)
            CALL SolverMapping_EquationsSetGet(solverMapping3D,1,coupledEquationsSet,err,error,*999)
            CALL NavierStokes_CalculateBoundaryFlux(equationsSet,coupledEquationsSet,iteration3D1D, &
              & convergedFlag,absolute3D0DTolerance,relative3D0DTolerance,err,error,*999)
          ENDIF
        CASE("1D Iterative Loop")
          !No longer using the NS/C coupling loop
          CALL ControlLoop_ContinueLoopSet(controlLoop,.FALSE.,err,error,*999)
        CASE DEFAULT
          localError="The iterative loop label of "//label// &
            & " does not correspond to a recognised loop type for a Navier-Stokes multiscale problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))// &
          & " is invalid for a Coupled 1D0D Navier-Stokes problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("NavierStokes_PostLoop")
    RETURN
999 ERRORSEXITS("NavierStokes_PostLoop",err,error)
    RETURN 1

  END SUBROUTINE NavierStokes_PostLoop

  !
  !================================================================================================================================
  !

  !>Updates boundary conditions for multiscale fluid problems
  SUBROUTINE NavierStokes_UpdateMultiscaleBoundary(equationsSet,boundaryConditions,timeIncrement,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    REAL(DP), INTENT(IN) :: timeIncrement
    INTEGER(INTG), INTENT(OUT) :: err
    TYPE(VARYING_STRING), INTENT(OUT) :: error
    !Local Variables
    INTEGER(INTG) :: boundaryConditionType,componentIdx,dependentDof,derivativeIdx,esSpecification(3),k,nodeIdx, &
      & numberOfVersions,numberOfLocalNodes,versionIdx
    REAL(DP) :: A0,A3D,ABoundary,ACellML,beta,E,H0,lengthScale,massScale,norm,normalWave(2,7),p3D,pCellml,pExternal,Q3D, &
      & qCellml,qBoundary,rho,timeScale,W1,W2
    REAL(DP), POINTER :: impedance(:),flow(:)
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainNodesType), POINTER :: dependentDomainNodes
    TYPE(DomainTopologyType), POINTER :: dependentDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(FieldType), POINTER :: dependentField,materialsField,independentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,uDependentVariable,u1DependentVariable,vDependentVariable, &
      & uIndependentVariable,uMaterialsVariable,vMaterialsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_UpdateMultiscaleBoundary",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    
    SELECT CASE(esSpecification(3))
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
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(uDependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,uDependentVariable,err,error,*999)
      NULLIFY(u1DependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U1_VARIABLE_TYPE,u1DependentVariable,err,error,*999)
      NULLIFY(vDependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,vDependentVariable,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(uIndependentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,uIndependentVariable,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      NULLIFY(vMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
      NULLIFY(decomposition)
      CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(decomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainNodes)
      CALL DomainTopology_DomainNodesGet(dependentDomainTopology,dependentDomainNodes,err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Navier-Stokes multiscale boundary update."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    SELECT CASE(esSpecification(3))
    !!!-- 1 D    E q u a t i o n s   S e t --!!!
    CASE(EQUATIONS_SET_TRANSIENT1D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_TRANSIENT1D_ADV_NAVIER_STOKES_SUBTYPE, &
       & EQUATIONS_SET_COUPLED1D0D_ADV_NAVIER_STOKES_SUBTYPE)

      CALL DomainNodes_NumberOfNodesGet(dependentDomainNodes,numberOfLocalNodes,err,error,*999)
      derivativeIdx=1
      !Get constant material parameters
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,2,rho,err,error,*999)
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,4,pExternal,err,error,*999)
      !Get materials scale factors
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,5,lengthScale,err,error,*999)
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,6,timeScale,err,error,*999)
      CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,7,massScale,err,error,*999)

      !!!--  L o o p   o v e r   l o c a l    n o d e s  --!!!
      DO nodeIdx=1,numberOfLocalNodes
        numberOfVersions=dependentDomainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions

        !Get normal wave direction
        normalWave=0.0_DP
        DO componentIdx=1,2
          DO versionIdx=1,numberOfVersions
            CALL FieldVariable_ParameterSetGetLocalNode(uIndependentVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,componentIdx,normalWave(componentIdx,versionIdx),err,error,*999)
          ENDDO !versionIdx
        ENDDO !componentIdx
        !!!-- F i n d   b o u n d a r y    n o d e s --!!!
        IF(ABS(normalWave(1,1)) > ZERO_TOLERANCE .OR. ABS(normalWave(2,1))> ZERO_TOLERANCE) THEN
          CALL L2Norm(normalWave(:,1),norm,err,error,*999)
          IF(numberOfVersions == 1 .AND. norm > ZERO_TOLERANCE) THEN
            versionIdx = 1
            !Get material parameters
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,1,A0,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,2,E,err,error,*999)
            CALL FieldVariable_ParameterSetGetLocalNode(vMaterialsVariable,FIELD_VALUES_SET_TYPE,versionIdx, &
             & derivativeIdx,nodeIdx,3,H0,err,error,*999)
            beta=(4.0_DP*(SQRT(PI))*E*H0)/(3.0_DP*A0)
            ! Get the boundary condition type for the dependent field primitive variables (Q,A)
            DO componentIdx=1,2
              CALL FieldVariable_LocalNodeDOFGet(uDependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx, &
                & dependentDOF,err,error,*999)
              NULLIFY(boundaryConditionsVariable)
              CALL BoundaryConditions_VariableGet(boundaryConditions,uDependentVariable,boundaryConditionsVariable,err,error,*999)
              CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,dependentDOF,boundaryConditionType, &
                & err,error,*999)
              SELECT CASE(boundaryConditionType)

              CASE(BOUNDARY_CONDITION_FIXED_NONREFLECTING)
                ! N o n - r e f l e c t i n g   B o u n d a r y
                ! ----------------------------------------------------
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  !Outlet - set W2 to 0, get W1 from the extrapolated value
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,1,W1,err,error,*999)
                  W2 = 0.0_DP
                ELSE
                  !Inlet - set W1 to 0, get W2 from the extrapolated value
                  W1 = 0.0_DP
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,2,W2,err,error,*999)
                ENDIF
                !Calculate new area value based on W1, W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                 & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area 1D non-reflecting boundary detected at node "// &
                    & TRIM(NumberToVString(nodeIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,ABoundary,err,error,*999)

              CASE(BOUNDARY_CONDITION_FIXED_CELLML)
                ! C o u p l e d   C e l l M L  ( 0 D )   B o u n d a r y
                ! ------------------------------------------------------------
                !Get qCellML used in pCellML calculation
!!TODO: should this be U1 variable???
                CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,1,QCellML,err,error,*999)
                !Get pCellML if this is a coupled problem
                CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,pCellml,err,error,*999)
                !Convert pCellML from SI base units specified in CellML file to scaled units (e.g., kg/(m.s^2) --> g/(mm.ms^2))
                !pCellml = pCellml*massScale/(lengthScale*(timeScale**2.0_DP))
                !Convert pCellML --> A0D
                ACellML=((pCellml-pExternal)/beta+SQRT(A0))**2.0_DP
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  !  O u t l e t
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,1,W1,err,error,*999)
                  !Calculate W2 from 0D domain
                  W2 = QCellml/ACellml - 4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP - A0**0.25_DP)
                ELSE
                  !  I n l e t
                  !Calculate W1 from 0D domain
                  W1 = QCellml/ACellml + 4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP - A0**0.25_DP)
                  !Calculate W2 from 1D domain
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,2,W2,err,error,*999)
                ENDIF
                !Calculate new area value based on W1,W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                  & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area coupled 1D0D boundary detected at node "// &
                    & TRIM(NumberToVString(nodeIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,ABoundary,err,error,*999)

              CASE(BOUNDARY_CONDITION_COUPLING_FLOW)
                ! C o u p l e d    3 D    B o u n d a r y
                ! ------------------------------------------------------------
                !Get q3D
                CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,1,Q3D,err,error,*999)
                QBoundary = Q3D
                !Set new Q value based on 3D value
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,1,QBoundary,err,error,*999)
              CASE(BOUNDARY_CONDITION_COUPLING_STRESS)
                !Get q3D
                CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,1,Q3D,err,error,*999)
                !Get p3D if this is a coupled problem
                CALL FieldVariable_ParameterSetGetLocalNode(u1DependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,p3D,err,error,*999)
                !Convert p3D --> A3D
                A3D=((p3D-pExternal)/beta+SQRT(A0))**2.0_DP
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  !  O u t l e t
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,1,W1,err,error,*999)
                  !Calculate W2 from 0D domain
                  W2 = Q3D/A3D - 4.0_DP*SQRT(beta/(2.0_DP*rho))*(A3D**0.25_DP - A0**0.25_DP)
                ELSE
                  !  I n l e t
                  !Calculate W1 from 0D domain
                  W1 = Q3D/A3D + 4.0_DP*SQRT(beta/(2.0_DP*rho))*(A3D**0.25_DP - A0**0.25_DP)
                  !Calculate W2 from 1D domain
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,2,W2,err,error,*999)
                ENDIF
                !Calculate new area value based on W1,W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                  & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                !DEBUG
                !ABoundary=A3D
                IF(ABoundary < ZERO_TOLERANCE) THEN
                  localError="Negative area coupled 3D1D boundary detected at node "// &
                    & TRIM(NumberToVString(nodeIdx,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Set new A value based on 3D boundary and update dof
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,ABoundary,err,error,*999)

              CASE(BOUNDARY_CONDITION_FIXED_STREE)
                ! S t r u c t u r e d   T r e e   B o u n d a r y
                ! ------------------------------------------------------------
               !Get qCellML used in pCellML calculation
                CALL FieldVariable_ParameterSetGetLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,1,QCellML,err,error,*999)
                !Get impedance function
                NULLIFY(impedance)
                CALL FieldVariable_ParameterSetDataGet(uMaterialsVariable,FIELD_VALUES_SET_TYPE,impedance,err,error,*999)
                !Get flow function
                NULLIFY(flow)
                CALL FieldVariable_ParameterSetDataGet(vMaterialsVariable,FIELD_VALUES_SET_TYPE,flow,err,error,*999)
                pCellml = 0.0_DP
                DO k=1,size(flow)
                  pCellml=pCellml+flow(k)*impedance(k)*timeIncrement
                ENDDO !k
                !Convert pCellML --> A0D
                ACellML=((pCellml-pExternal)/beta+SQRT(A0))**2.0_DP
                IF(normalWave(1,1) > ZERO_TOLERANCE) THEN
                  !  O u t l e t
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,1,W1,err,error,*999)
                  !Calculate W2 from 0D domain
                  W2 = QCellml/ACellml-4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP-A0**0.25_DP)
                ELSE
                  !  I n l e t
                  !Calculate W1 from 0D domain
                  W1 = QCellml/ACellml+4.0_DP*SQRT(beta/(2.0_DP*rho))*(ACellml**0.25_DP-A0**0.25_DP)
                  !Calculate W2 from 1D domain
                  CALL FieldVariable_ParameterSetGetLocalNode(vDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                    & nodeIdx,2,W2,err,error,*999)
                ENDIF
                !Calculate new area value based on W1, W2 and update dof
                ABoundary = (((2.0_DP*rho)/(beta))**2.0_DP)* &
                 & (((W1-W2)/8.0_DP+SQRT(beta/(2.0_DP*rho))*((A0)**0.25_DP))**4.0_DP)
                CALL FieldVariable_ParameterSetUpdateLocalNode(uDependentVariable,FIELD_VALUES_SET_TYPE,versionIdx,derivativeIdx, &
                  & nodeIdx,2,ABoundary,err,error,*999)

              CASE(BOUNDARY_CONDITION_NONE, &
                & BOUNDARY_CONDITION_FIXED, &
                & BOUNDARY_CONDITION_FIXED_INLET, &
                & BOUNDARY_CONDITION_FIXED_OUTLET, &
                & BOUNDARY_CONDITION_FIXED_FITTED)
                !Do nothing
              CASE DEFAULT
                localError="The boundary conditions type "//TRIM(NumberToVString(boundaryConditionType,"*",err,error))// &
                  & " is not valid for a coupled characteristic problem."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO ! componentIdx

          ENDIF ! boundary node
        ENDIF ! branch or boundary node
      ENDDO !Loop over nodes
      !Update distributed fields
      CALL FieldVariable_ParameterSetUpdateStart(uDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uDependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      
    CASE(EQUATIONS_SET_STATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_LAPLACE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_QUASISTATIC_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE)
      !-- 3 D    E q u a t i o n s   S e t --!!!
      ! Do nothing
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
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

  !>Calculate the fluid flux through 3D boundaries for use in problems with coupled solutions (e.g. multidomain)
  SUBROUTINE NavierStokes_CalculateBoundaryFlux3D0D(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryID,boundaryType,componentIdx,computationNode,coupledNodeNumber,dependentBasisType, &
      & dependentBasisType2,elementIdx,elementNodeIdx,esSpecification(3),faceIdx,faceNodeIdx,faceNodeDerivativeIdx, &
      & faceNumber,gaussIdx,groupCommunicator,meshComponentNumber,mpiIError,nodeNumber,numberOfBoundaries,numberOfDimensions, &
      & numberOfFaceGauss,numberOfGlobalBoundaries,numberOfGroupComputationNodes,numberOfLocalElements,numberOfLocalFaces, &
      & numberOfNodes,numberOfNodeDerivatives,orientation,totalNumberOfLocalElements,userElementNumber,versionNumber, &
      & xiDirection(3),xiNormalDirection
    REAL(DP) :: boundaryValue,boundaryValueTemp,courant,elementNormal(3),faceArea,faceNormal(3),facePressure,faceTraction, &
      & faceVelocity,gaussWeight,globalBoundaryArea(10),globalBoundaryFlux(10),globalBoundaryMeanPressure(10), &
      & globalBoundaryPressure(10),jacobian,jacobianGaussWeight,localBoundaryArea(10),localBoundaryFlux(10), &
      & localBoundaryNormalStress(10),localBoundaryPressure(10),maxCourant,mu,muScale,normalDifference,normalTolerance,p0D, &
      & pressureGauss,q0D,toleranceCourant,unitNormal(3),velocityGauss(3)
    LOGICAL :: boundaryFace,boundary3D0DFound(10),convergedFlag !<convergence flag for 3D-0D coupling
    LOGICAL, ALLOCATABLE :: globalConverged(:)
    TYPE(BasisType), POINTER :: dependentBasis,dependentBasis2,faceBasis
    TYPE(DecompositionType), POINTER :: decomposition3D
    TYPE(DecompositionElementType), POINTER :: decompositionElement
    TYPE(DecompositionElementsType), POINTER :: decompositionElements3D
    TYPE(DecompositionFaceType), POINTER :: face
    TYPE(DecompositionFacesType), POINTER :: decompositionFaces3D
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology3D
    TYPE(DomainType), POINTER :: domain,domain3D
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainMappingType), POINTER :: elementsMapping3D
    TYPE(DomainMappingsType), POINTER :: domainMappings3D
    TYPE(DomainTopologyType), POINTER :: domainTopology,domainTopology3D
    TYPE(EquationsType), POINTER :: equations3D
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField3D,equationsSetField3D,geometricField3D,materialsField
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable3D,uDependentVariable3D,uMaterialsVariable,uEquationsSetVariable3D, &
      & u1EquationsSetVariable3D,vEquationsSetVariable3D
    TYPE(QuadratureSchemeType), POINTER :: faceQuadratureScheme
    TYPE(VARYING_STRING) :: localError
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("NavierStokes_CalculateBoundaryFlux3D0D",err,error,*999)

    boundary3D0DFound = .FALSE.
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    ! 3 D   t y p e s :   I n t e g r a t e   b o u n d a r y   v a l u e s
    ! ------------------------------------------------------------------------
    CASE(EQUATIONS_SET_MULTISCALE3D_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_STATIC_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_TRANSIENT_RBS_NAVIER_STOKES_SUBTYPE, &
      & EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      ! Get 3D field pointers
      NULLIFY(equations3D)
      CALL EquationsSet_EquationsGet(equationsSet,equations3D,err,error,*999)
      NULLIFY(geometricField3D)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField3D,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField3D,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(equationsSetField3D)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField3D,err,error,*999)
      NULLIFY(uEquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_U_VARIABLE_TYPE,uEquationsSetVariable3D,err,error,*999)
      NULLIFY(u1EquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_U1_VARIABLE_TYPE,u1EquationsSetVariable3D,err,error,*999)
      NULLIFY(vEquationsSetVariable3D)
      CALL Field_VariableGet(equationsSetField3D,FIELD_V_VARIABLE_TYPE,vEquationsSetVariable3D,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations3D,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      NULLIFY(nonlinearMapping)
      CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
      NULLIFY(dependentVariable3D)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable3D,err,error,*999)
      NULLIFY(dependentField3D)
      CALL FieldVariable_FieldGet(dependentVariable3D,dependentField3D,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      NULLIFY(decomposition3D)
      CALL Field_DecompositionGet(dependentField3D,decomposition3D,err,error,*999)
      NULLIFY(decompositionTopology3D)
      CALL Decomposition_DecompositionTopologyGet(decomposition3D,decompositionTopology3D,err,error,*999)
      NULLIFY(decompositionElements3D)
      CALL DecompositionTopology_DecompositionElementsGet(decompositionTopology3D,decompositionElements3D,err,error,*999)
      NULLIFY(decompositionFaces3D)
      CALL DecompositionTopology_DecompositionFacesGet(decompositionTopology3D,decompositionFaces3D,err,error,*999)
      NULLIFY(domain3D)
      CALL Decomposition_DomainGet(decomposition3D,0,domain3D,err,error,*999)
      NULLIFY(domainTopology3D)
      CALL Domain_DomainTopologyGet(domain3D,domainTopology3D,err,error,*999)
      NULLIFY(domainMappings3D)
      CALL Domain_DomainMappingsGet(domain3D,domainMappings3D,err,error,*999)
      NULLIFY(elementsMapping3D)
      CALL DomainMappings_ElementsMappingGet(domainMappings3D,elementsMapping3D,err,error,*999)
      ! Get constant max Courant (CFL) number (default 1.0)
      CALL FieldVariable_ParameterSetGetConstant(u1EquationsSetVariable3D,FIELD_VALUES_SET_TYPE,2,toleranceCourant,err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_CONSTITUTIVE_MU_NAVIER_STOKES_SUBTYPE) THEN
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,1,muScale,err,error,*999)
      ELSE
        CALL FieldVariable_ParameterSetGetConstant(uMaterialsVariable,FIELD_VALUES_SET_TYPE,1,mu,err,error,*999)
      ENDIF

      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations3d,equationsInterpolation,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
        & err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
        & err,error,*999)

      ! Loop over elements to locate boundary elements
      maxCourant = 0.0_DP
      numberOfBoundaries = 0
      localBoundaryFlux = 0.0_DP
      localBoundaryArea = 0.0_DP
      localBoundaryPressure = 0.0_DP
      localBoundaryNormalStress = 0.0_DP
      CALL DomainMapping_NumberOfLocalGet(elementsMapping3D,numberOfLocalElements,err,error,*999)
      DO elementIdx=1,numberOfLocalElements
        CALL DecompositionElements_ElementUserNumberGet(decompositionElements3D,elementIdx,userElementNumber,err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable3D,1,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(domainFaces)
        CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        NULLIFY(dependentBasis)
        CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis,err,error,*999)
        CALL Basis_TypeGet(dependentBasis,dependentBasisType,err,error,*999)
 
        ! Note: if CFL tolerance = 0, we'll skip this step, which speeds things up a bit
        IF (toleranceCourant > ZERO_TOLERANCE) THEN
          ! C F L  c o n d i t i o n   c h e c k
          ! ------------------------------------
          !Calculate element metrics (courant #, cell Reynolds number)
          CALL NavierStokes_CalculateElementMetrics(equationsSet,elementIdx,err,error,*999)
          !Get element metrics
          CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,3,courant, &
            & err,error,*999)
          IF(courant < -ZERO_TOLERANCE) THEN
            localError="The Courant (CFL) number of "//TRIM(NumberToVString(courant,"*",err,error))// &
              & " is negative for element number "//TRIM(NumberToVString(userElementNumber,"*",err,error))//"."
            CALL FlagWarning(localError,err,error,*999)
          ENDIF
          IF(courant > maxCourant) maxCourant = courant
          !Check if element CFL number below specified tolerance
          IF(courant > toleranceCourant) THEN
            localError="Element "//TRIM(NumberToVString(userElementNumber,"*",err,error))// &
              & " has violated the CFL condition ("//TRIM(NumberToVString(courant,"*",err,error))//" <= "// &
              & TRIM(NumberToVString(toleranceCourant,"*",err,error))// &
              & "). Decrease timestep or increase CFL tolerance for the 3D Navier-Stokes problem."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        ENDIF

        ! B o u n d a r y   n o r m a l   a n d   I D
        ! ----------------------------------------------
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,5, &
          & elementNormal(1),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,6, &
          & elementNormal(2),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,7, &
          & elementNormal(3),err,error,*999)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,8, &
          & boundaryValueTemp,err,error,*999)
        boundaryID=NINT(boundaryValueTemp)
        CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,9, &
          & boundaryValue,err,error,*999)
        boundaryType=NINT(boundaryValue)
        !Check if is a non-wall boundary element
        IF(boundaryID > numberOfBoundaries) numberOfBoundaries=boundaryID
        IF(boundaryID>1) THEN
          faceArea=0.0_DP
          faceVelocity=0.0_DP
          facePressure=0.0_DP
          faceTraction=0.0_DP
          !Loop over faces to determine the boundary face contribution
          CALL Basis_NumberOfLocalFacesGet(dependentBasis,numberOfLocalFaces,err,error,*999)
          DO faceIdx=1,numberOfLocalFaces
            !Get the face normal and quadrature information
            CALL DecompositionElements_ElementFaceNumberGet(decompositionElements3D,faceIdx,elementIdx,faceNumber,err,error,*999)
            CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces3D,faceNumber,boundaryFace,err,error,*999)
            
            !This speeds things up but is also important, as non-boundary faces have an XI_DIRECTION that might
            !correspond to the other element.
            IF(.NOT.(boundaryFace)) CYCLE

            xiDirection = 0.0_DP
            SELECT CASE(dependentBasisType)
            CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
              CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces3D,faceNumber,xiNormalDirection,err,error,*999)
              xiDirection(3)=ABS(xiNormalDirection)
            CASE(BASIS_SIMPLEX_TYPE)
              CALL FlagWarning("Boundary flux calculation not yet set up for simplex element types.",err,error,*999)
            CASE DEFAULT
              localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType,"*",err,error))// &
                & " is not yet implemented for Navier-Stokes."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            NULLIFY(faceBasis)
            CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)
            NULLIFY(faceQuadratureScheme)
            CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_NumberOfGaussGet(faceQuadratureScheme,numberOfFaceGauss,err,error,*999)
 
            !Use the geometric field to find the face normal and Jacobian for the face integral
            CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,geometricInterpParameters,err,error,*999)
            CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,dependentInterpParameters,err,error,*999)
 
            xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
            xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
            orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*xiNormalDirection)
            !Loop over face gauss points
            DO gaussIdx=1,numberOfFaceGauss
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,geometricInterpPoint, &
                & err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_AREA_TYPE,geometricInterpPointMetrics, &
                & err,error,*999)

              CALL BasisQuadratureScheme_GaussWeightGet(faceQuadratureScheme,gaussIdx,gaussWeight,err,error,*999)
              CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
              jacobianGaussWeight=jacobian*gaussWeight
              
              !Make sure this is the boundary face that corresponds with boundaryID (could be a wall rather than inlet/outlet)
              CALL CrossProduct(geometricInterpPointMetrics%dXdXi(:,1),geometricInterpPointMetrics%dXdXi(:,2),faceNormal, &
                & err,error,*999)
              faceNormal = faceNormal*orientation
              CALL Normalise(faceNormal,unitNormal,err,error,*999)
              CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
              normalTolerance=0.1_DP
              IF(normalDifference>normalTolerance) EXIT

              !Get interpolated velocity and pressure
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussIdx,dependentInterpPoint, &
                & err,error,*999)
              velocityGauss(1:numberOfDimensions)=dependentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
              pressureGauss=dependentInterpPoint%values(numberOfDimensions+1,NO_PART_DERIV)

              ! I n t e g r a t e    f a c e   a r e a ,   v e l o c i t y   a n d   p r e s s u r e
              ! ----------------------------------------------------------------------------------------
              faceArea=faceArea + jacobianGaussWeight
              facePressure=facePressure + pressureGauss*jacobianGaussWeight
              DO componentIdx=1,numberOfDimensions
                faceVelocity=faceVelocity+velocityGauss(componentIdx)*unitNormal(componentIdx)*jacobianGaussWeight
              ENDDO !componentIdx
            ENDDO !gaussIdx
          ENDDO !faceIdx
          localBoundaryFlux(boundaryID) = localBoundaryFlux(boundaryID) + faceVelocity
          localBoundaryArea(boundaryID) = localBoundaryArea(boundaryID) + faceArea
          localBoundaryPressure(boundaryID) = localBoundaryPressure(boundaryID) + facePressure
        ENDIF !boundaryIdentifier
      END DO !elementIdx
      !Distribute any updated element fields
      CALL FieldVariable_ParameterSetUpdateStart(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)

      ! G a t h e r   v a l u e s   o v e r   t h r e a d s
      ! ------------------------------------------------------
      !Need to add boundary flux for any boundaries split accross computation nodes
      numberOfGlobalBoundaries = 0
      globalBoundaryFlux = 0.0_DP
      globalBoundaryArea = 0.0_DP
      globalBoundaryPressure = 0.0_DP
      NULLIFY(workGroup)
      CALL Decomposition_WorkGroupGet(decomposition3D,workGroup,err,error,*999)
      CALL WorkGroup_GroupCommunicatorGet(workGroup,groupCommunicator,err,error,*999)
      CALL WorkGroup_NumberOfGroupNodesGet(workGroup,numberOfGroupComputationNodes,err,error,*999)
      CALL WorkGroup_GroupNodeNumberGet(workGroup,computationNode,err,error,*999)
      IF(numberOfGroupComputationNodes>1) THEN !use mpi
        CALL MPI_ALLREDUCE(localBoundaryFlux,globalBoundaryFlux,10,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryArea,globalBoundaryArea,10,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(localBoundaryPressure,globalBoundaryPressure,10,MPI_DOUBLE_PRECISION,MPI_SUM,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
        CALL MPI_ALLREDUCE(numberOfBoundaries,numberOfGlobalBoundaries,1,MPI_INTEGER,MPI_MAX,groupCommunicator,mpiIError)
        CALL MPI_ErrorCheck("MPI_ALLREDUCE",mpiIError,err,error,*999)
      ELSE
        numberOfGlobalBoundaries = numberOfBoundaries
        globalBoundaryFlux = localBoundaryFlux
        globalBoundaryArea = localBoundaryArea
        globalBoundaryPressure = localBoundaryPressure
      ENDIF
      globalBoundaryArea=ABS(globalBoundaryArea)
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          globalBoundaryMeanPressure(boundaryID)=globalBoundaryPressure(boundaryID)/globalBoundaryArea(boundaryID)
        ENDIF
      ENDDO !boundaryID
      DO boundaryID=2,numberOfGlobalBoundaries
        IF(globalBoundaryArea(boundaryID) > ZERO_TOLERANCE) THEN
          IF(diagnostics1) THEN
            IF(computationNode==0) THEN
              CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  flow:  ", &
                & globalBoundaryFlux(boundaryID),err,error,*999)
              CALL WriteStringTwoValue(DIAGNOSTIC_OUTPUT_TYPE,"3D boundary ",boundaryID,"  mean pressure:  ", &
                & globalBoundaryMeanPressure(boundaryID),err,error,*999)
              IF(toleranceCourant > ZERO_TOLERANCE) THEN
                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Max Courant (CFL) number: ",maxCourant,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ELSE
          localError="Zero or negative area boundary detected on boundary "//TRIM(NumberToVString(boundaryID,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDDO !boundaryID

    CASE DEFAULT
      localError="Boundary flux calcluation for equations type "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not yet implemented for Navier-Stokes."
      CALL FlagError(localError,err,error,*999)
    END SELECT

!!TODO: WHY IS THE CODE BELOW NOT INSIDE THE CASE SELECTION ABOVE AS THE CASE DEFAULT IS A FLAG ERROR AND SO WE MUST HAVE THE SAME EQUATIONS SET TO GET HERE?    

    
    ! C o p y    i n t e g r a t e d   v a l u e s    t o    t a r g e t    f i e l d s
    ! ------------------------------------------------------------------------------------
    convergedFlag = .TRUE.
    !Loop over elements again to allocate flux terms to boundary nodes
    DO elementIdx=1,totalNumberOfLocalElements
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,5,elementNormal(1), &
        & err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,6,elementNormal(2), &
        & err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,7,elementNormal(3), &
        & err,error,*999)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,8,boundaryValue, &
        & err,error,*999)
      boundaryID=NINT(boundaryValue)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,9,boundaryValue, &
        & err,error,*999)
      boundaryType=NINT(boundaryValue)
      CALL FieldVariable_ParameterSetGetLocalElement(vEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,elementIdx,11,boundaryValue, &
        & err,error,*999)
      coupledNodeNumber=NINT(boundaryValue)
      IF(boundaryID>1) THEN
        meshComponentNumber=2
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition3D,meshComponentNumber,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(domainFaces)
        CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        NULLIFY(dependentBasis2)
        CALL DomainElements_ElementBasisGet(domainElements,elementIdx,dependentBasis2,err,error,*999)
        CALL Basis_TypeGet(dependentBasis2,dependentBasisType2,err,error,*999)
        
        decompositionElement=>decompositionElements3D%elements(elementIdx)
   
        ! B o u n d a r y   F a c e    N o r m a l s
        ! --------------------------------------------------
        CALL Basis_NumberOfLocalFacesGet(dependentBasis2,numberOfLocalFaces,err,error,*999)
        DO faceIdx=1,numberOfLocalFaces
          !Get the face normal and quadrature information
          CALL DecompositionElements_ElementFaceNumberGet(decompositionElements3D,faceIdx,elementIdx,faceNumber,err,error,*999)
          CALL DecompositionFaces_FaceBoundaryFaceGet(decompositionFaces3D,faceNumber,boundaryFace,err,error,*999)
          
          IF(.NOT.(boundaryFace)) CYCLE
          
          !TODO: this sort of thing should be moved to a more general Basis_FaceNormalGet (or similar) routine
          xiDirection = 0.0_DP
          SELECT CASE(dependentBasisType2)
          CASE(BASIS_LAGRANGE_HERMITE_TP_TYPE)
            CALL DecompositionFaces_FaceXiNormalDirectionGet(decompositionFaces3D,faceNumber,xiNormalDirection,err,error,*999)
            xiDirection(3)=ABS(xiNormalDirection)
          CASE(BASIS_SIMPLEX_TYPE)
            CALL FlagWarning("Boundary flux calculation not yet set up for simplex element types.",err,error,*999)
          CASE DEFAULT
            localError="Face integration for basis type "//TRIM(NumberToVString(dependentBasisType2,"*",err,error))// &
              & " is not yet implemented for Navier-Stokes."
            CALL FlagError(localError,err,error,*999)
          END SELECT

          NULLIFY(faceBasis)
          CALL DomainFaces_FaceBasisGet(domainFaces,faceNumber,faceBasis,err,error,*999)
          NULLIFY(faceQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(faceBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,faceQuadratureScheme,err,error,*999)
          !Use the geometric field to find the face normal and Jacobian for the face integral
          CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,faceNumber,geometricInterpParameters,err,error,*999)
 
          xiDirection(1)=OTHER_XI_DIRECTIONS3(xiDirection(3),2,1)
          xiDirection(2)=OTHER_XI_DIRECTIONS3(xiDirection(3),3,1)
          orientation=SIGN(1,OTHER_XI_ORIENTATIONS3(xiDirection(1),xiDirection(2))*xiNormalDirection)
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,1,geometricInterpPoint,err,error,*999)
          CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_AREA_TYPE,geometricInterpPointMetrics,err,error,*999)
          CALL CrossProduct(geometricInterpPointMetrics%dXdXi(:,1),geometricInterpPointMetrics%dXdXi(:,2),faceNormal, &
            & err,error,*999)
          faceNormal = faceNormal*orientation
          CALL Normalise(faceNormal,unitNormal,err,error,*999)
          CALL L2Norm(elementNormal-unitNormal,normalDifference,err,error,*999)
          normalTolerance=0.1_DP
          IF(normalDifference>normalTolerance) CYCLE

          ! U p d a t e    N o d a l   V a l u e s
          ! --------------------------------------------------
          ! Update local nodes with integrated boundary flow values
          CALL Basis_NumberOfLocalNodesGet(faceBasis,numberOfNodes,err,error,*999)
          DO faceNodeIdx=1,numberOfNodes
            CALL Basis_FaceNodeNumberGet(dependentBasis2,faceNodeIdx,faceIdx,elementNodeIdx,err,error,*999)
            CALL DomainElements_ElementNodeGet(domainElements,elementNodeIdx,elementIdx,nodeNumber,err,error,*999)
            CALL Basis_NodeNumberOfDerivativesGet(faceBasis,faceNodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO faceNodeDerivativeIdx=1,numberOfNodeDerivatives
              versionNumber=1
              IF(boundaryType==BOUNDARY_CONDITION_FIXED_CELLML) THEN
                !Check current values against those passed to the CellML solver
                CALL FieldVariable_ParameterSetGetLocalNode(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,1,1,nodeNumber,1,q0D, &
                  & err,error,*999)
                CALL FieldVariable_ParameterSetGetLocalNode(dependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,1,1,nodeNumber, &
                  & 4,p0D,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVariable3D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,4,p0D,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalNode(uEquationsSetVariable3D,FIELD_PREVIOUS_ITERATION_VALUES_SET_TYPE, &
                  & versionNumber,faceNodeDerivativeIdx,nodeNumber,1,q0D,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalNode(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,versionNumber, &
                  & faceNodeDerivativeIdx,nodeNumber,1,globalBoundaryFlux(boundaryID),err,error,*999)
              ENDIF
            ENDDO !faceNodeDerivativeIdx
          ENDDO !faceNodeIdx
        ENDDO !faceIdx
      ENDIF !boundaryIdentifier
    ENDDO !elementIdx

!!TODO: WHAT IS THE CONVERGED FLAG? IT IS NOT PASSED IN ANYWHERE?
    
    !allocate array for mpi communication
    IF(numberOfGroupComputationNodes>1) THEN !use mpi
      ALLOCATE(globalConverged(numberOfGroupComputationNodes),STAT=err) 
      IF(err/=0) CALL FlagError("Could not allocate global convergence check array.",err,error,*999)
      CALL MPI_ALLGATHER(convergedFlag,1,MPI_LOGICAL,globalConverged,1,MPI_LOGICAL,groupCommunicator,mpiIError)
      CALL MPI_ErrorCheck("MPI_ALLGATHER",mpiIError,err,error,*999)
      IF(ALL(globalConverged)) THEN
        convergedFlag = .TRUE.
      ELSE
        convergedFlag = .FALSE.
      ENDIF
    ENDIF

    !Distribute any updated fields
    IF(ASSOCIATED(equationsSetField3D)) THEN
      CALL FieldVariable_ParameterSetUpdateStart(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(uEquationsSetVariable3D,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDIF
    IF(ASSOCIATED(dependentField3D)) THEN
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable3D,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
    ENDIF

    EXITS("NavierStokes_CalculateBoundaryFlux3D0D")
    RETURN
999 ERRORSEXITS("NavierStokes_CalculateBoundaryFlux3D0D",err,error)
    RETURN 1
    
  END SUBROUTINE NavierStokes_CalculateBoundaryFlux3D0D

  !
  !================================================================================================================================
  !

  !>Calculates the wall shear stress for fluid problems at all the boundary nodes
  SUBROUTINE NavierStokes_WallShearStressCalculate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The equations set to calculate the wall shear stress for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryIdx,coordinateIdx,elementNumber,esSpecification(3),faceIdx,lineIdx,localNodeIdx, &
      & meshComponentNumber,nodeNumber,nodeIdx,numberOfBoundaries,numberOfDimensions,numberOfNodes,numberOfXi, &
      & startNode,stopNode,xiIdx
    REAL(DP) :: boundaryXi(2),deludelxi(3),fullXi(3),gradu,mu,normal(3),position(3),tangents(3,2),wss
    LOGICAL :: boundaryNode,found
    TYPE(BasisType), POINTER :: boundaryBasis,velocityBasis
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DecompositionTopologyType), POINTER :: decompositionTopology
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainFacesType), POINTER :: domainFaces
    TYPE(DomainFaceType), POINTER :: face
    TYPE(DomainLinesType), POINTER :: domainLines
    TYPE(DomainLineType), POINTER :: line
    TYPE(DomainMappingType), POINTER :: nodesMappings
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("NavierStokes_WallShearStressCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ALE_NAVIER_STOKES_SUBTYPE, &
      &  EQUATIONS_SET_ALE_RBS_NAVIER_STOKES_SUBTYPE)

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentVariable)
      CALL Field_VariableExists(dependentField,FIELD_W_VARIABLE_TYPE,dependentVariable,err,error,*999)
      IF(ASSOCIATED(dependentVariable)) THEN
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        NULLIFY(materialsField)
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        NULLIFY(equations)
        CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        NULLIFY(vectorMapping)
        CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
        NULLIFY(residualMapping)
        CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
        NULLIFY(dependentVariable)
        CALL EquationsMappingResidual_VariableGet(residualMapping,1,dependentVariable,err,error,*999)
        NULLIFY(decomposition)
        CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
        NULLIFY(decompositionTopology)
        CALL Decomposition_DecompositionTopologyGet(decomposition,decompositionTopology,err,error,*999)
        CALL FieldVariable_ComponentMeshComponentGet(dependentVariable,1,meshComponentNumber,err,error,*999)
        NULLIFY(domain)
        CALL Decomposition_DomainGet(decomposition,meshComponentNumber,domain,err,error,*999)
        NULLIFY(domainMappings)
        CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainElements)
        CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
        NULLIFY(nodesMappings)
        CALL DomainMappings_NodesMappingGet(domainMappings,nodesMappings,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)

        IF(numberOfDimensions==2) THEN
          CALL Decomposition_AssertCalculateLines(decomposition,err,error,*999)
          NULLIFY(domainLines)
          CALL DomainTopology_DomainLinesGet(domainTopology,domainLines,err,error,*999)
        ELSE IF(numberOfDimensions==3) THEN
          CALL Decomposition_AssertCalculateFaces(decomposition,err,error,*999)
          NULLIFY(domainFaces)
          CALL DomainTopology_DomainFacesGet(domainTopology,domainFaces,err,error,*999)
        ELSE
          localError="The number of dimensions of "//TRIM(NumberToVString(numberofDimensions,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        ENDIF

        NULLIFY(equationsInterpolation)
        CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
        NULLIFY(geometricInterpParameters)
        CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpParameters, &
          & err,error,*999)
        NULLIFY(geometricInterpPoint)
        CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
          & err,error,*999)
        NULLIFY(geometricInterpPointMetrics)
        CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpPointMetrics,err,error,*999)
        NULLIFY(dependentInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpParameters, &
          & err,error,*999)
        NULLIFY(dependentInterpPoint)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
          & err,error,*999)
        NULLIFY(materialsInterpParameters)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpParameters, &
          & err,error,*999)
        NULLIFY(materialsInterpPoint)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
          & err,error,*999)

        !Loop over internal and boundary nodes
        CALL DomainMapping_InternalStartGet(nodesMappings,startNode,err,error,*999)
        CALL DomainMapping_BoundaryFinishGet(nodesMappings,stopNode,err,error,*999)
        !Loop over internal and boundary nodes
        nodes: DO nodeIdx=startNode,stopNode
          CALL DomainMapping_NumberGet(nodesMappings,nodeIdx,nodeNumber,err,error,*999)
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeNumber,boundaryNode,err,error,*999)
          IF(.NOT.boundaryNode) CYCLE nodes
          !Node is on the boundary. Loop over surrounding faces/lines
          IF(numberOfDimensions == 2) THEN
            CALL DomainNodes_NodeNumberOfLinesGet(domainNodes,nodeNumber,numberOfBoundaries,err,error,*999)
          ELSE
            CALL DomainNodes_NodeNumberOfFacesGet(domainNodes,nodeNumber,numberOfBoundaries,err,error,*999)
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
!!TODO: This bit needs to be redone. A boundary line could be involved in multiple elements. The code should look at the surrounding elements via the decomposition lines?
            IF(numberOfDimensions==2) THEN
              CALL DomainNodes_NodeLineNumberGet(domainNodes,boundaryIdx,nodeNumber,lineIdx,err,error,*999)
              NULLIFY(line)
              CALL DomainLines_LineGet(domainLines,lineIdx,line,err,error,*999)
              NULLIFY(boundaryBasis)
              CALL DomainLines_LineBasisGet(domainLines,lineIdx,boundaryBasis,err,error,*999)
              elementNumber=line%elementNumber
              NULLIFY(velocityBasis)
              CALL DomainElements_ElementBasisGet(domainElements,elementNumber,velocityBasis,err,error,*999)
              !Find node position
              found=.FALSE.
              CALL Basis_NumberOfLocalNodesGet(boundaryBasis,numberOfNodes,err,error,*999)
              DO localNodeIdx=1,numberOfNodes
                IF(line%nodesInLine(localNodeIdx)==nodeNumber) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !localNodeIdX
              IF(.NOT.found) THEN
                localError="Could not find node number "//TRIM(NumberToVString(nodeNumber,"*",err,error))// &
                  & " in line number "//TRIM(NumberToVString(lineIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
                & err,error,*999)
            ELSE
              CALL DomainNodes_NodeFaceNumberGet(domainNodes,boundaryIdx,nodeNumber,faceIdx,err,error,*999)
              NULLIFY(face)
              CALL DomainFaces_FaceGet(domainFaces,faceIdx,face,err,error,*999)
              NULLIFY(boundaryBasis)
              CALL DomainFaces_FaceBasisGet(domainFaces,faceIdx,boundaryBasis,err,error,*999)
              elementNumber=face%elementNumber
              NULLIFY(velocityBasis)
              CALL DomainElements_ElementBasisGet(domainElements,elementNumber,velocityBasis,err,error,*999)
              !Find node position
              found=.FALSE.
              CALL Basis_NumberOfLocalNodesGet(boundaryBasis,numberOfNodes,err,error,*999)
              DO localNodeIdx=1,numberOfNodes
                IF(face%nodesInFace(localNodeIdx)==nodeNumber) THEN
                  found=.TRUE.
                  EXIT
                ENDIF
              ENDDO !localNodeIdX
              IF(.NOT.found) THEN
                localError="Could not find node number "//TRIM(NumberToVString(nodeNumber,"*",err,error))// &
                  & " in face number "//TRIM(NumberToVString(faceIdx,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
                & err,error,*999)
            ENDIF
            CALL Field_InterpolateXi(FIRST_PART_DERIV,boundaryXi,geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(boundaryBasis%numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_PositionNormalTangentsCalculateIntPtMetric(geometricInterpPointMetrics,.FALSE.,position,normal,tangents, &
              & err,error,*999)
            CALL Basis_NumberOfXiGet(velocityBasis,numberOfXi,err,error,*999)
            CALL Basis_LocalNodeXiCalculate(boundaryBasis,localNodeIdx,boundaryXi,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters, &
              & err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters, &
              & err,error,*999)
            CALL Basis_BoundaryXiToXi(velocityBasis,boundaryIdx,boundaryXi,fullXi,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,fullXi,geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(velocityBasis%numberOfXi,geometricInterpPointMetrics,err,error,*999)
            CALL Field_InterpolateXi(FIRST_PART_DERIV,fullXi,dependentInterpPoint,err,error,*999)
            CALL Field_InterpolateXi(NO_PART_DERIV,fullXi,materialsInterpPoint,err,error,*999)
            mu=materialsInterpPoint%values(1,NO_PART_DERIV)
            DO coordinateIdx=1,numberOfDimensions
              DO xiIdx=1,numberOfXi
                deludelxi(xiIdx)=dependentInterpPoint%values(coordinateIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
              ENDDO !xiIdx
              gradu=DOT_PRODUCT(deludelxi(1:numberOfXi),geometricInterpPointMetrics%dXidX(1:numberOfXi,coordinateIdx))
              wss=wss+mu*gradu*normal(coordinateIdx)
            ENDDO !coordinateIdx
          ENDDO !boundaryIdx
          wss=wss/REAL(numberOfBoundaries)
          CALL Field_ParameterSetUpdateLocalNode(dependentField,FIELD_W_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,nodeNumber,1, &
            & wss,err,error,*999)
        ENDDO nodes !nodeIdx
      ENDIF
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
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

END MODULE NavierStokesEquationsRoutines
