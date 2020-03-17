!> \file
!> \author Chris Bradley
!> \brief This module handles all diffusion equation routines.
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

!>This module handles all advection-diffusion equation routines.
MODULE AdvectionDiffusionEquationsRoutines

  USE AnalyticAnalysisRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types
! temporary input for setting velocity field
  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE 

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC AdvectionDiffusion_EquationsSetSetup

  PUBLIC AdvectionDiffusion_EquationsSetSolnMethodSet
  
  PUBLIC AdvectionDiffusion_BoundaryConditionsAnalyticCalculate

  PUBLIC AdvectionDiffusion_EquationsSetSpecificationSet
  
  PUBLIC AdvectionDiffusion_FiniteElementCalculate

  PUBLIC AdvectionDiffusion_ProblemSpecificationSet
  
  PUBLIC AdvectionDiffusion_ProblemSetup

  PUBLIC AdvectionDiffusion_PreSolveUpdateInputData

  PUBLIC AdvectionDiffusion_PreSolveGetSourceValue
  
  PUBLIC AdvectionDiffusion_PreSolveStoreCurrentSoln

  PUBLIC AdvectionDiffusion_PreSolve

  PUBLIC AdvectionDiffusion_PostSolve
  
CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !>For the advection-diffusion analytic example it is required that the advective velocity
  !>and the source field are set to a particular analytic value, which is performed within this subroutine.
  SUBROUTINE AdvectionDiffusion_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,derivativeIdx,dimensionIdx,globalDerivativeIndex,localDOFIdx, &
      & nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodeDerivatives,numberOfVariables,numberOfNodes,variableIdx, &
      & variableType
    REAL(DP) :: alpha,dependentValue,x(3),sourceValue,independentValue,materialValue,Peclet,phi,tanphi
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField,sourceField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,independentVariable,materialVariable,sourceVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate",err,error,*999)

!!TODO: do this properly with analytic field values.    
    alpha = 1.0_DP
    phi = 0.2_DP
    tanphi = TAN(phi)
    Peclet= 10.0_DP

    !>Set the analytic boundary conditions
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
          DO dimensionIdx=1,numberOfDimensions
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
              !This is a steady-state solution of advection-diffusion equation
              !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
              !Solution is u(x,y)=tanh(1 - alpha.(x.tan(Phi) - y))
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=TANH(1.0-alpha*(x(1)*tanphi-x(2)))
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implmented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=0.0_DP
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
              & dependentValue,err,error,*999)
            IF(variableType==FIELD_U_VARIABLE_TYPE.AND.boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
              CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentField,variableType,localDOFIdx, &
                & BOUNDARY_CONDITION_FIXED,dependentValue,err,error,*999)
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    
    !>Set the independent field (i.e. the advective velocity) to a specified analytical function
    NULLIFY(indendentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    CALL Field_NumberOfVariablesGet(indendentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(independentVariable)
      CALL Field_VariableIndexGet(independentField,variableIdx,independentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(independentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,umberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(independentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(independentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,domainNodes%numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
          DO dimensionIdx=1,numberOfDimensions
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
              !Velocity field takes form v(x,y)=(sin(6y),cos(6x))
              IF(componentIdx==1) THEN
                independentValue=SIN(6.0_DP*x(1))
              ELSE
                independentValue=COS(6.0_DP*x(2))           
              ENDIF
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(independentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(independentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
              & independentValue,err,error,*999)
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL Field_ParameterSetUpdateStart(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx

    !>Set the source field to a specified analytical function
    NULLIFY(sourceField)
    CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
    CALL Field_NumberOfVariablesGet(sourceField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(sourceVariable)
      CALL Field_VariableIndexGet(sourceField,variableIdx,sourceVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(sourceVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(independentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(independentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
          DO dimensionIdx=1,numberOfDimensions
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
              sourceValue=(1.0_DP/Peclet)*(2.0_DP*TANH(-0.1E1_DP+alpha*(tanphi*x(1)-x(2)))*(1.0_DP-(TANH(-0.1E1_DP+ &
                & alpha*(tanphi*x(1)-x(2)))**2))*alpha*alpha*tanphi*tanphi+2.0_DP*TANH(-0.1E1_DP+ &
                & alpha*(tanphi*x(1)-x(2)))*(1.0_DP-(TANH(-0.1E1_DP+alpha*(tanphi*x(1)-x(2)))**2))*alpha*alpha- &
                & Peclet*(-SIN(6.0_DP*x(2))*(1.0_DP-(TANH(-0.1E1_DP+alpha*(tanphi*x(1)-x(2)))**2))*alpha*tanphi+ &
                & COS(6.0_DP*x(1))*(1.0_DP-(TANH(-0.1E1_DP+alpha*(tanphi*x(1)-x(2)))**2))*alpha))
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(independentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(sourceVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,sourceValue, &
              & err,error,*999)
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx

    !>Set the material field to a specified analytical value
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    CALL Field_NumberOfVariablesGet(materialsField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(materialsVariable)
      CALL Field_VariableIndexGet(materialsField,variableIdx,materialsVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(materialsVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        SELECT CASE(analyticFunctionType)
        CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
          materialValue= (1.0/Peclet)
        CASE DEFAULT
          localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        CALL FieldVariable_ParameterSetUpdateConstant(materialsVariable,FIELD_VALUES_SET_TYPE,componentIdx,materialValue, &
          & err,error,*999)
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(materialsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(materialsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)

    EXITS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("AdvectionDiffusion_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_BoundaryConditionsAnalyticCalculate


  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSetup",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(equationsSet%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_equationsSet_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_equationsSet_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL AdvectionDiffusion_EquationsSetLinearSetup(equationsSet,equationsSetSetup,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_EquationsSetSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSetup",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSolnMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, & 
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)        
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
        & " is not valid for an advection-diffusion equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("AdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSolnMethodSet",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSolnMethodSet")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetSolnMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a diffusion equation type of a classical field equations set class.
  SUBROUTINE AdvectionDiffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for an advection-diffusion type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TYPE,subtype]

    EXITS("AdvectionDiffusion_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetSpecificationSet",err,error)
    EXITS("AdvectionDiffusion_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the advection-diffusion equation.
  SUBROUTINE AdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3),numberOfVariables,numberOfVariablesCount,numberOfMaterialsCouplingComponents    
    INTEGER(INTG) :: esFieldNumberOfVariables,esFieldNumberOfComponents    
    INTEGER(INTG), POINTER :: esFieldData(:)
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:),variableUTypes(:),couplingMatrixStorageTypes(:), &
      & couplingMatrixStructureTypes(:)
    INTEGER(INTG) :: equationsSetSubtype
    INTEGER(INTG) :: componentIdx,geometricMeshComponent,geometricScalingType,numberOfDimensions, &
      & numberOfMaterialsComponents,numberOfSourceComponents,numberOfIndependentComponents,myMatrixIdx,numberOfCompartments,&
      & geometricComponentNumber,numberOfIndependentUComponents,numberOfIndependentVComponents,variableIdx
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetEquationsSetFieldType), POINTER :: equationsEquationsSetField
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetSourceType), POINTER :: equationsSource
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,equationsSetFieldField
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    equationsSetSubtype=esSpecification(3)
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid for an advection-diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
      
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL AdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          esFieldNumberOfVariables = 1
          esFieldNumberOfComponents = 2
          NULLIFY(equationsEquationsSetField)
          CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsEquationsSetField,err,error,*999)
          IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
            NULLIFY(region)
            CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
            !Create the auto created equations set field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsEquationsSetField%equationsSetFieldField, &
              & err,error,*999)
            CALL Field_LabelSet(equationsEquationsSetField%equationsSetFieldField,"Equations Set Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsEquationsSetField%equationsSetFieldField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsEquationsSetField%equationsSetFieldField,FIELD_INDEPENDENT_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsEquationsSetField%equationsSetFieldField,esFieldNumberOfVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsEquationsSetField%equationsSetFieldField,[FIELD_U_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsEquationsSetField%equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsEquationsSetField%equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsEquationsSetField%equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
              & esFieldNumberOfComponents,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,esFieldNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,esFieldNumberOfComponents, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          IF(equationsSet%equationsSetField%equationsSetFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSet%equationsSetField%equationsSetFieldField,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsSet%equationsSetField%equationsSetFieldField,&
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1_INTG,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsSet%equationsSetField%equationsSetFieldField,&
              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,1_INTG,err,error,*999)
          ENDIF
        ENDIF
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for an advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          esFieldNumberOfComponents = 2
          NULLIFY(equationsEquationsSetField)
          CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsEquationsSetField,err,error,*999)
          IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsEquationsSetField%equationsSetFieldField,geometricDecomposition, &
              & err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsEquationsSetField%equationsSetFieldField,geometricField,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
            DO componentIdx=1,esFieldNumberOfComponents
              CALL Field_ComponentMeshComponentSetAndLock(equationsSet%equationsSetField%equationsSetFieldField, &
                & FIELD_U_VARIABLE_TYPE,componentIdx,geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%equationsSetField%equationsSetFieldField, &
                & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            END DO
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%equationsSetField%equationsSetFieldField,geometricScalingType,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          ! do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Do nothing???
      CASE(EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
        !ALE types
        CALL Field_ParameterSetEnsureCreate(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetEnsureCreate(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
      CASE DEFAULT
        !Do nothing
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & err,error,*999)
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            !Default the scaling to the geometric field scaling
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
          SELECT CASE(equationsSetSubtype)
          CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
            & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
            !uses number of compartments to check that appropriate number and type of variables have been set on the dependent field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)                 
            equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
            CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
            numberOfCompartments=esFieldData(2)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2*numberOfCompartments,err,error,*999)
            !Create & populate array storing all of the relevant variable types against which to check the field variables
            ALLOCATE(variableTypes(2*numberOfCompartments))
            DO variableIdx=1,numberOfCompartments
              variableTypes(2*variableIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
              variableTypes(2*variableIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(variableIdx-1))
            ENDDO !variableIdx
            CALL Field_VariableTypesCheck(equationsSetSetup%field,variableTypes,err,error,*999)
            DO variableIdx=1,2*numberOfCompartments
              CALL Field_DimensionCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,variableTypes(variableIdx),FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,variableTypes(variableIdx),1,err,error,*999)
            ENDDO !variableIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)            
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              componentIdx=1
              DO variableIdx=1,2*numberOfCompartments
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,variableTypes(variableIdx),componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !variableIdx
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
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfDimensions, &
              & err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)            
            SELECT CASE(equationsSet%solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !                   DO componentIdx=1,numberOfDimensions
              componentIdx=1
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              !                   ENDDO !componentIdx
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
          END SELECT
        ENDIF
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsSet%dependent% &
            & dependentField,err,error,*999)
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition, &
            & err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,equationsSet%geometry% &
            & geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,4,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField, &
            & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & 1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & 1,err,error,*999)
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_DELVDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                   !Default the scaling to the geometric field scaling
            CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
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
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE,&
            & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
          
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,1,err,error,*999)
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !                  DO componentIdx=1,numberOfDimensions
            componentIdx=1
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
              & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE, & 
              & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            !                   ENDDO !componentIdx
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
            localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      END SELECT
    CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
      IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
        CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
          & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
        & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a linear advection-diffusion equation"
      CALL FlagError(localError,err,error,*999)
    END SELECT
  CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
    !-----------------------------------------------------------------
    ! M a t e r i a l s   f i e l d
    !-----------------------------------------------------------------
    SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            equationsMaterials=>equationsSet%MATERIALS
            IF(ASSOCIATED(equationsMaterials)) THEN
              IF(equationsMaterials%materialsFieldAutoCreated) THEN
               IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsMaterials% &
                    & materialsField,err,error,*999)
                  CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
                  CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
                  CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition, &
                    & err,error,*999)
                  CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,equationsSet%geometry% &
                    & geometricField,err,error,*999)
                  CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,2,err,error,*999)
                  CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                    numberOfMaterialsComponents=numberOfDimensions
                   !Set the number of materials components
                  CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfMaterialsComponents,err,error,*999)
                  equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                  CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                    numberOfCompartments=esFieldData(2)
                    numberOfMaterialsCouplingComponents=numberOfCompartments
                  CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                    & numberOfMaterialsCouplingComponents,err,error,*999)
                  !Default the k materials components to the geometric interpolation setup with constant interpolation
                  DO componentIdx=1,numberOfDimensions
                    CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,geometricMeshComponent,err,error,*999)
                  ENDDO !componentIdx
                    CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                      & 1,geometricMeshComponent,err,error,*999)
                  DO componentIdx=1,numberOfMaterialsCouplingComponents
                    CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                      & componentIdx,geometricMeshComponent,err,error,*999)
                  ENDDO 
                    !Default the field scaling to that of the geometric field
                  CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
               ELSE !standard materials field
                !Create the auto created materials field
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsMaterials% &
                  & materialsField,err,error,*999)
                CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,equationsSet%geometry% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  numberOfMaterialsComponents=numberOfDimensions
                ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  numberOfMaterialsComponents=numberOfDimensions+1
                ELSE
                  numberOfMaterialsComponents=numberOfDimensions
                ENDIF
                 !Set the number of materials components
                CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfMaterialsComponents,err,error,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO componentIdx=1,numberOfDimensions
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricMeshComponent,err,error,*999)   
                  DO componentIdx=numberOfDimensions+1,numberOfMaterialsComponents
                    CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,geometricMeshComponent,err,error,*999)
                  ENDDO !componentIdx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
               ENDIF
              ELSE
               IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                  CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
                  CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
                     & FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                    & err,error,*999)
                  equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                  CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                    numberOfCompartments=esFieldData(2)
                    numberOfMaterialsCouplingComponents=numberOfCompartments
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE, &
                    & numberOfMaterialsCouplingComponents,err,error,*999)
               ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                    & err,error,*999)
                ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1, &
                    & err,error,*999)
                ENDIF
               ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            equationsMaterials=>equationsSet%MATERIALS
            IF(ASSOCIATED(equationsMaterials)) THEN
              IF(equationsMaterials%materialsFieldAutoCreated) THEN
                !Finish creating the materials field
                CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
                !Set the default values for the materials field
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  numberOfMaterialsComponents=numberOfDimensions             
                ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                  !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                  numberOfMaterialsComponents=numberOfDimensions+1
                ELSE
                  numberOfMaterialsComponents=numberOfDimensions
                ENDIF
                !First set the k values to 1.0
                DO componentIdx=1,numberOfDimensions
                  !WRITE(*,'("Setting materials components values :")')
                  CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                ENDDO !componentIdx
                IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                  CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                    numberOfCompartments=esFieldData(2)
                  DO componentIdx=1,numberOfCompartments
                   CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                     & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
                   ENDDO !componentIdx
                ENDIF
                IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  !Now set the linear source values to 1.0
                  DO componentIdx=numberOfDimensions+1,numberOfMaterialsComponents
                    CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                  ENDDO !componentIdx
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsSource=>equationsSet%SOURCE
            IF(ASSOCIATED(equationsSource)) THEN
              IF(equationsSource%sourceFieldAutoCreated) THEN
                !Create the auto created source field
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsSource% &
                  & sourceField,err,error,*999)
                CALL Field_LabelSet(equationsSource%sourceField,"Source Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,equationsSet%geometry% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                numberOfSourceComponents=1
                !Set the number of source components
                CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfSourceComponents,err,error,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &  
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  DO componentIdx=1,numberOfSourceComponents
                    CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                     & componentIdx,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                      & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                ENDIF
                  !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            equationsSource=>equationsSet%SOURCE
            IF(ASSOCIATED(equationsSource)) THEN
              IF(equationsSource%sourceFieldAutoCreated) THEN
                !Finish creating the source field
                CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
                !Set the default values for the source field
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  numberOfSourceComponents=1
                ELSE
                  numberOfSourceComponents=0
                ENDIF
                !Now set the source values to 1.0
                DO componentIdx=1,numberOfSourceComponents
                  CALL Field_ComponentValuesInitialise(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          !Setup the equations set for the advective velocity field
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            equationsIndependent=>equationsSet%INDEPENDENT
            IF(ASSOCIATED(equationsIndependent)) THEN
             IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              IF(equationsIndependent%independentFieldAutoCreated) THEN
                !Create the auto created independent field
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsIndependent% &
                  & independentField,err,error,*999)
                CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField, &
                      & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,equationsSet%geometry% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                  numberOfIndependentUComponents=numberOfDimensions
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfIndependentUComponents,err,error,*999)
                equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                  numberOfCompartments=esFieldData(2)
                  numberOfIndependentVComponents=numberOfCompartments-1
                 !Set the number of independent components
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                  & numberOfIndependentVComponents,err,error,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO componentIdx=1,numberOfIndependentUComponents
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                DO componentIdx=1,numberOfIndependentVComponents
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                  !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],&
                   & err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                   & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                  & err,error,*999)
                equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                   & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                  numberOfCompartments=esFieldData(2)
                  numberOfIndependentVComponents=numberOfCompartments-1
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                  & numberOfIndependentVComponents,err,error,*999)
                DO componentIdx=1,numberOfIndependentUComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO
                DO componentIdx=1,numberOfIndependentVComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO
              ENDIF
             ELSE
              IF(equationsIndependent%independentFieldAutoCreated) THEN
                !Create the auto created independent field
                CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsIndependent% &
                  & independentField,err,error,*999)
                CALL Field_LabelSet(equationsIndependent%independentField,"Independent Field",err,error,*999)
                CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField, &
                      & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
                CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,equationsSet%geometry% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                  numberOfIndependentComponents=numberOfDimensions
                 !Set the number of independent components
                CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfIndependentComponents,err,error,*999)
                !Default the k independent components to the geometric interpolation setup with constant interpolation
                DO componentIdx=1,numberOfDimensions
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                  !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
                  & err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDIF
             ENDIF
            ELSE
              CALL FlagError("Equations set independent is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            equationsIndependent=>equationsSet%INDEPENDENT
            IF(ASSOCIATED(equationsIndependent)) THEN
              IF(equationsIndependent%independentFieldAutoCreated) THEN
                !Finish creating the independent field
                CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
                !Set the default values for the independent field
                CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)

                 IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                  & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
!                   CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
!                      & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
                 ENDIF
!                 CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
!                   & numberOfDimensions,err,error,*999)
!                   numberOfIndependentComponents=numberOfDimensions             
!                 !First set the k values to 1.0
!                 DO componentIdx=1,numberOfIndependentComponents
!                   CALL Field_ComponentValuesInitialise(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
!                ENDDO !componentIdx
              ENDIF
            ELSE
              CALL FlagError("Equations set independent is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)              
            dependentField=>equationsSet%dependent%dependentField
            IF(ASSOCIATED(dependentField)) THEN
              geometricField=>equationsSet%geometry%geometricField
              IF(ASSOCIATED(geometricField)) THEN
                CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
                SELECT CASE(equationsSetSetup%analyticFunctionType)
                CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
                  IF(numberOfDimensions/=2) THEN
                    localError="The number of geometric dimensions of "// &
                      & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                      & " is invalid. The analytic function type of "// &
                      & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                      & " requires that there be 2 geometric dimensions."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                  equationsSet%ANALYTIC%analyticFunctionType= & 
                    & EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1
                CASE DEFAULT
                  localError="The specified analytic function type of "// &
                    & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                    & " is invalid for a linear advection-diffusion equation."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
              analyticField=>equationsSet%ANALYTIC%analyticField
              IF(ASSOCIATED(analyticField)) THEN
                IF(equationsSet%ANALYTIC%analyticFieldAutoCreated) THEN
                  CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSetSetup%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
            CALL Equations_CreateStart(equationsSet,EQUATIONS,err,error,*999)
            CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
            IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
              & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            ELSEIF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_STATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set subtype not valid.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(equationsSet%solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. & 
                & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                !Finish the equations
                CALL EquationsSet_EquationsGet(equationsSet,EQUATIONS,err,error,*999)
                CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                SELECT CASE(equationsSetSubtype)
                CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
                  & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
                  & EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE,&
                  & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                  CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                  CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                  CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                  CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                  IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                    & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                    CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
                    CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
                  ENDIF
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                  & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)              
                  equationsSetFieldField=>equationsSet%equationsSetField%equationsSetFieldField
                  CALL Field_ParameterSetDataGet(equationsSetFieldField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
                  myMatrixIdx = esFieldData(1)
                  numberOfCompartments = esFieldData(2)    
                  CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                  CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                  CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,numberOfCompartments-1,err,error,*999)
                  
                  ALLOCATE(variableTypes(2*numberOfCompartments))
                  ALLOCATE(variableUTypes(numberOfCompartments-1))
                  DO numberOfVariables=1,numberOfCompartments
                    variableTypes(2*numberOfVariables-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(numberOfVariables-1))
                    variableTypes(2*numberOfVariables)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(numberOfVariables-1))
                  ENDDO
                  numberOfVariablesCount=0
                  DO numberOfVariables=1,numberOfCompartments
                    IF(numberOfVariables/=myMatrixIdx)THEN
                      numberOfVariablesCount=numberOfVariablesCount+1
                      variableUTypes(numberOfVariablesCount)=variableTypes(2*numberOfVariables-1)
                    ENDIF
                  ENDDO
                  CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx-1),err,error,*999)
                  CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,variableUTypes,err,error,*999)
                  CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx),err,error,*999)
                  CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
                  CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
                END SELECT
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                !Set up matrix storage and structure
                IF(EQUATIONS%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                  !Set up lumping
                  CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                ELSE
                  SELECT CASE(EQUATIONS%sparsityType)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                    CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices, &
                      & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                      & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                      & err,error,*999)
                    CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                      [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)    
                    IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
                      & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)THEN
                      ALLOCATE(couplingMatrixStorageTypes(numberOfCompartments-1))
                      ALLOCATE(couplingMatrixStructureTypes(numberOfCompartments-1))
                      DO numberOfVariables=1,numberOfCompartments-1
                        couplingMatrixStorageTypes(numberOfVariables)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                        couplingMatrixStructureTypes(numberOfVariables)=EQUATIONS_MATRIX_FEM_STRUCTURE
                      ENDDO
                      CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,couplingMatrixStorageTypes,err,error,*999)
                      CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,couplingMatrixStructureTypes, &
                        & err,error,*999)
                    ENDIF
                  CASE DEFAULT
                    localError="The equations matrices sparsity type of "// &
                      & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
              ELSE IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                !Finish the creation of the equations
                NULLIFY(equations)
                CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
                CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                  & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                  & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                  & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                  CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
                  CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
                ENDIF
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                    & err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
          & " is not a linear advection-diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_EquationsSetLinearSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_EquationsSetLinearSetup",err,error)
    EXITS("AdvectionDiffusion_EquationsSetLinearSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_EquationsSetLinearSetup

  !
  !===============================================================================================================================
  !
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP(equationsSet,equationsSetSetup,err,error,*)
!     !Argument variables
!     TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
!     TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
!     INTEGER(INTG), INTENT(OUT) :: err !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
! 
!     CALL FlagError("Not implemented.",err,error,*999)
!     EXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP")
!     RETURN
! 999 ERRORSEXITS("ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP",err,error)
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_EQUATIONS_SET_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE AdvectionDiffusion_ProblemSetup(PROBLEM,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_ProblemSetup",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
!        CALL ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,problemSetup,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_ProblemSetup")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE AdvectionDiffusion_FiniteElementCalculate(equationsSet,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) mh,mhs,ms,ng,nh,nhs,ni,nj,ns,FIELD_VAR_TYPE,my_compartment,numberOfCompartments,imatrix,numberOfVariablesCount
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),ADVEC_VEL,A_PARAM,COUPLING_PARAM,PGM,PGN
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField,sourceField,independentField,EQUATIONS_SET_FIELD
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE,geometricVariable
    TYPE(FieldVariablePtrType) :: FIELD_VARIABLES(99)
    TYPE(EquationsMatrixPtrType) :: COUPLING_MATRICES(99) 
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    TYPE(QuadratureSchemeType), POINTER :: QUADRATURE_SCHEME
    TYPE(QuadratureSchemeType), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(VARYING_STRING) :: localError
    TYPE(FieldInterpolationParametersType), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS
    TYPE(FieldInterpolatedPointType), POINTER :: DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT
    INTEGER(INTG), POINTER :: esFieldData(:)
    LOGICAL :: updateDampingMatrix,updateStiffnessMatrix,updateRHSVector,updateSourceVector
    INTEGER(INTG) :: equationsSetSubtype

    updateDampingMatrix = .FALSE.
    updateStiffnessMatrix = .FALSE.
    updateRHSVector = .FALSE.
    updateSourceVector = .FALSE.

    ENTERS("AdvectionDiffusion_FiniteElementCalculate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>equationsSet%EQUATIONS
      
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(equationsSet%SPECIFICATION,1)<3) THEN
          CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
        END IF
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        equationsSetSubtype=equationsSet%SPECIFICATION(3)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          dependentField=>equations%interpolation%dependentField
          geometricField=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          independentField=>equations%interpolation%independentField !Stores the advective velocity field
          IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             sourceField=>equations%interpolation%sourceField
          ENDIF
          vectorMatrices=>vectorEquations%vectorMatrices
          IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             dynamicMatrices=>vectorMatrices%dynamicMatrices
             stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
             dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          ELSEIF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             linearMatrices=>vectorMatrices%linearMatrices
             stiffnessMatrix=>linearMatrices%matrices(1)%ptr
          ELSE

          ENDIF
          rhsVector=>vectorMatrices%rhsVector
          IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN            
             sourceVector=>vectorMatrices%sourceVector
          ENDIF
          IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. & 
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             IF(ASSOCIATED(dampingMatrix)) updateDampingMatrix=dampingMatrix%updateMatrix
          ENDIF
          IF(ASSOCIATED(stiffnessMatrix)) updateStiffnessMatrix=stiffnessMatrix%updateMatrix
          IF(ASSOCIATED(rhsVector)) updateRHSVector=rhsVector%updateVector
          IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN             
             IF(ASSOCIATED(sourceVector)) updateSourceVector=sourceVector%updateVector
          ENDIF
          vectorMapping=>vectorEquations%vectorMapping
          IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
           & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
           EQUATIONS_SET_FIELD=>equationsSet%equationsSetField%equationsSetFieldField
           CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
             & FIELD_VALUES_SET_TYPE,esFieldData,err,error,*999)
            my_compartment = esFieldData(1)
            numberOfCompartments  = esFieldData(2)
            linearMatrices=>vectorMatrices%linearMatrices
            linearMapping=>vectorMapping%linearMapping
            numberOfVariablesCount=0
            DO imatrix = 1,numberOfCompartments
             IF(imatrix/=my_compartment)THEN
              numberOfVariablesCount=numberOfVariablesCount+1
              COUPLING_MATRICES(numberOfVariablesCount)%ptr=>linearMatrices%matrices(numberOfVariablesCount)%ptr
              FIELD_VARIABLES(numberOfVariablesCount)%ptr=>linearMapping%equationsMatrixToVarMaps(numberOfVariablesCount)%VARIABLE
              FIELD_VAR_TYPES(numberOfVariablesCount)=FIELD_VARIABLES(numberOfVariablesCount)%ptr%variableType
              COUPLING_MATRICES(numberOfVariablesCount)%ptr%elementMatrix%MATRIX=0.0_DP
             ENDIF
            END DO
          ENDIF
          IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &  
             & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             dynamicMapping=>vectorMapping%dynamicMapping
             FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
          ELSEIF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
             linearMapping=>vectorMapping%linearMapping
             FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
          ELSE

          ENDIF
          FIELD_VAR_TYPE=FIELD_VARIABLE%variableType
          geometricVariable=>geometricField%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%decomposition%meshComponentNumber)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%decomposition%meshComponentNumber)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
               & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF

          IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF
          !the following line has been changed to use fieldinputdata1settype
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN  
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
            & equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
           CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,err,error,*999)
           DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
            & equations%interpolation%dependentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
          ENDIF


          !Select whether using standard Galerkin scheme, or the stabilised streamwise-upwinding Petrov-Galerkin scheme
          SELECT CASE(equationsSetSubtype)
          CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,& 
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE , & 
               & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%numberOfGauss
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%numberOfXi,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Interpolate to get the advective velocity
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)            
            IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
             & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)   
            ENDIF
            IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                 & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
              & QUADRATURE_SCHEME%gaussWeights(ng)
            !Loop over field components

            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,geometricVariable%numberOfComponents
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%numberOfXi                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%gaussBasisFunctions(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%gaussBasisFunctions(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                        ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. & 
                            & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                          A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(geometricVariable%numberOfComponents,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ELSEIF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                        !for multi-compartment model must include additional terms into the
                          COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                            & VALUES(my_compartment,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+ & 
                            & SUM*RWG + QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG*COUPLING_PARAM
                        ENDIF
                      ENDIF
                    IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN 
                IF(updateSourceVector) THEN
                    C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%numberOfComponents
                    !DO mh=1,DEPENDENT_VARIABLE%numberOfComponents
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ELSEIF(equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                IF(updateSourceVector) THEN
                    CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                          & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,err,error,*999)
                    C_PARAM=DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%numberOfComponents
                    !DO mh=1,DEPENDENT_VARIABLE%numberOfComponents
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=0.0_DP
  
            IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
            !Calculate the coupling matrices

              !Loop over element rows
              mhs=0
              DO mh=1,FIELD_VARIABLE%numberOfComponents !field_variable is the variable associated with the equations set under consideration

                MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%meshComponentNumber
                DEPENDENT_BASIS_1 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                  & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                  & quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian * &
                  & QUADRATURE_SCHEME_1%gaussWeights(ng)

                DO ms=1,DEPENDENT_BASIS_1%numberOfElementParameters
                  mhs=mhs+1

                  numberOfVariablesCount=0
                  DO imatrix = 1,numberOfCompartments
                  IF(imatrix/=my_compartment)THEN
                    numberOfVariablesCount=numberOfVariablesCount+1

!need to test for the case where imatrix==mycompartment
!the coupling terms then needs to be added into the stiffness matrix
                    IF(COUPLING_MATRICES(numberOfVariablesCount)%ptr%updateMatrix) THEN

!                       !Loop over element columns
                      nhs=0
! !                       DO nh=1,FIELD_VARIABLE%numberOfComponents
                      DO nh=1,FIELD_VARIABLES(numberOfVariablesCount)%ptr%numberOfComponents

                        MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%meshComponentNumber
                        DEPENDENT_BASIS_2 => dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        !--- We cannot use two different quadrature schemes here !!!
                        QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                         & quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        !RWG = equations%interpolation%geometricInterpPointMetrics%jacobian * &
                        !  & QUADRATURE_SCHEME_2%gaussWeights(ng)

                        DO ns=1,DEPENDENT_BASIS_2%numberOfElementParameters
                          nhs=nhs+1

!                           !-------------------------------------------------------------------------------------------------------------
!                           !concentration test function, concentration trial function
!                           !For now, this is only a dummy implementation - this still has to be properly set up.
!                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN ! don't need this for diffusion equation

!                             SUM = 0.0_DP

                            PGM=QUADRATURE_SCHEME_1%gaussBasisFunctions(ms,NO_PART_DERIV,ng)
                            PGN=QUADRATURE_SCHEME_2%gaussBasisFunctions(ns,NO_PART_DERIV,ng)

                            !Get the coupling coefficients 
                              COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                                & VALUES(imatrix,NO_PART_DERIV)

!                              SUM = SUM + COUPLING_PARAM * PGM * PGN
 
                             COUPLING_MATRICES(numberOfVariablesCount)%ptr%elementMatrix%matrix(mhs,nhs) = &
                               & COUPLING_MATRICES(numberOfVariablesCount)%ptr%elementMatrix%matrix(mhs,nhs) + & 
                               & COUPLING_PARAM * PGM * PGN * RWG
!                           ENDIF
 
                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                   ENDIF
                  ENDDO !imatrix
                ENDDO !ms
              ENDDO !mh

            ENDIF

          ENDDO !ng
          
          !Scale factor adjustment
          IF(dependentField%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1                    
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF
                    IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. & 
                       & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                IF(updateSourceVector) sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
               & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999) 
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%numberOfGauss
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%numberOfXi,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Interpolate to get the advective velocity
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)            
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
              & QUADRATURE_SCHEME%gaussWeights(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,geometricVariable%numberOfComponents
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%numberOfXi                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%gaussBasisFunctions(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%gaussBasisFunctions(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          !Advection term is constructed here and then added to SUM for updating the stiffness matrix outside of this loop
                          ADVEC_VEL=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+ADVEC_VEL*QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*PGNJ(nj)   
                        ENDDO !nj
                        IF (equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                        ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                          A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(geometricVariable%numberOfComponents,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG
                            ! A_PARAM is the material parameter that multiplies the linear source u
                        ENDIF
                      ENDIF
                    IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
               ENDDO !ms
              ENDDO !mh
              IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
                IF(updateSourceVector) THEN
                    C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
                    mhs=0
                    DO mh=1,FIELD_VARIABLE%numberOfComponents
                    !DO mh=1,DEPENDENT_VARIABLE%numberOfComponents
                     !Loop over element rows
                      DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                       mhs=mhs+1
                       sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                      ENDDO !ms
                    ENDDO !mh
                ENDIF
              ENDIF
            IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=0.0_DP
          ENDDO !ng
          
          !Scale factor adjustment
          IF(dependentField%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1                    
                nhs=0
                IF(updateStiffnessMatrix .OR. updateDampingMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(updateStiffnessMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF
                    IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                       & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE) THEN
                      IF(updateDampingMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF
                    ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(updateRHSVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              IF(equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_CONSTANT_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                 & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                IF(updateSourceVector) sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
              ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
              & " is not valid for an advection-diffusion equation type of a classical field equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
             & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
          CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
            & " is not valid for an advection-diffusion equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("AdvectionDiffusion_FiniteElementCalculate")
    RETURN
999 ERRORS("AdvectionDiffusion_FiniteElementCalculate",err,error)
    EXITS("AdvectionDiffusion_FiniteElementCalculate")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !INSERT CODE HERE TO DEAL WITH THE NON-LINEAR SOURCE TERMS - DECIDE HOW TO SOLVE THEM, USE JACOBIAN AS FOR NON-LINEAR POISSON EXAMPLE?

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for an advection diffusion problem type.
  SUBROUTINE AdvectionDiffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("AdvectionDiffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
          & PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for an advection-diffusion type of a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Advection-diffusion problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("AdvectionDiffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("AdvectionDiffusion_ProblemSpecificationSet",err,error)
    EXITS("AdvectionDiffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE AdvectionDiffusion_ProblemLinearSetup(PROBLEM,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS
    TYPE(SolversType), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: PROBLEM_SUBTYPE
    
    ENTERS("AdvectionDiffusion_ProblemLinearSetup",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
      END IF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      IF(PROBLEM_SUBTYPE==PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
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
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_NumberOfSolversSet(SOLVERS,1,err,error,*999)
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
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
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SolverEquationsType)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSEIF(PROBLEM_SUBTYPE==PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
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
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_NumberOfSolversSet(SOLVERS,1,err,error,*999)
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
             !Start the linear solver creation
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SolverEquationsType)
          SELECT CASE(problemSetup%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
              & " is invalid for a linear static advection-diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
            & " is invalid for a linear static advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
          & " does not equal a linear advection-diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("AdvectionDiffusion_ProblemLinearSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_ProblemLinearSetup",err,error)
    EXITS("AdvectionDiffusion_ProblemLinearSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemLinearSetup
  
  !
  !================================================================================================================================
  !
    !>Sets up the diffusion equations.
!   SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,problemSetup,err,error,*)
! 
!     !Argument variables
!     TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to setup
!     TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
!     INTEGER(INTG), INTENT(OUT) :: err !<The error code
!     TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
! 
!     CALL FlagError("Not implemented.",err,error,*999) 
! 
!     EXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
!     RETURN
! 999 ERRORSEXITS("ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",err,error)
!     RETURN 1
! 
!   END SUBROUTINE ADVECTION_DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP

  !
  !================================================================================================================================
  !

  !>Performs pre-solve operations for advection-diffusion problems.
  SUBROUTINE AdvectionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubType
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an advection-diffusion problem.",err,error,*999)
    
    problemSubType=problem%specification(3)
    IF(problemSubType==PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
      & problemSubType==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE) THEN
       CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
       !Update independent data fields
       CALL AdvectionDiffusion_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
       !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(controlLoop,solver,err,error,*999)
     ELSE IF(problemSubType==PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
       & problemSubType==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
       CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
       !Update independent data fields
       CALL AdvectionDiffusion_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
       !CALL ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(controlLoop,solver,err,error,*999)
       CALL AdvectionDiffusion_PreSolveALEUpdateMesh(controlLoop,solver,err,error,*999)
     ELSE IF(problemSubType==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
       & problemSubType==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN         
       CALL FlagError("Not implemented.",err,error,*999)
     ELSE IF(problemSubType==PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
       & problemSubType==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
       !do nothing
     ELSE
       localError="The third problem specification of "// &
         & TRIM(NumberToVString(problemSubType,"*",err,error))// &
         & " is not valid for a advection-diffusion type of a classical field problem."
       CALL FlagError(localError,err,error,*999)
     ENDIF

    EXITS("AdvectionDiffusion_PreSolve")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolve
  !   
  !================================================================================================================================
  !
  !>Update mesh position and velocity for ALE advection-diffusion problem
  SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh(SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldType), POINTER :: geometricField
    TYPE(SolverType), POINTER :: SOLVER_ALE_DIFFUSION !<A pointer to the solvers
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)

    INTEGER(INTG) :: dof_number,totalNumberOfDofs,NDOFS_TO_PRINT

    INTEGER(INTG) :: numberOfDimensions
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)

    ENTERS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(SOLVER_ALE_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%loopType==CONTROL_TIME_LOOP_TYPE) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ELSE IF(CONTROL_LOOP%controlLoopLevel>1) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP%parentLoop,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ENDIF
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%solverEquations
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)<3) THEN
                        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
                      END IF
                      SELECT CASE(equationsSet%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
                          ! do nothing ???
                        CASE(EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_GENERALISED_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_QUAD_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE, &
                          & EQUATIONS_SET_EXP_SOURCE_ALE_ADVECTION_DIFF_SUPG_SUBTYPE)
                          CALL WriteString(GENERAL_OUTPUT_TYPE,"Advection-diffusion update mesh ... ",err,error,*999)
                          geometricField=>equationsSet%geometry%geometricField
                          IF(ASSOCIATED(geometricField)) THEN
                            !--- First, read mesh displacement values from file

                           CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                            & numberOfDimensions,err,error,*999)

                           INPUT_TYPE=42
                           INPUT_OPTION=2
                           NULLIFY(INPUT_DATA1)
                           !CALL Field_ParameterSetDataGet(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, & 
                            !& FIELD_VALUES_SET_TYPE,INPUT_DATA1,err,error,*999)
                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & numberOfDimensions,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%timeLoop%iterationNumber,1.0_DP, &
                            & err,error,*999)

                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & err,error,*999)
                            ENDIF

                           CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                            & numberOfDimensions,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%timeLoop%iterationNumber,1.0_DP, &
                            & err,error,*999)

                            totalNumberOfDofs = geometricField%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & totalNumberOfDofs

                            !--- Second, update geometric field
                            DO dof_number=1,totalNumberOfDofs
                              CALL Field_ParameterSetAddLocalDOF(geometricField, & 
                                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, & 
                                & MESH_DISPLACEMENT_VALUES(dof_number), &
                                & err,error,*999)
                            END DO
                            CALL Field_ParameterSetUpdateStart(geometricField, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetUpdateFinish(geometricField, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)

                            !--- Third, use displacement values to calculate velocity values
                            ALPHA=1.0_DP/TIME_INCREMENT
                            CALL Field_ParameterSetsCopy(geometricField,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
                            CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE, & 
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                          ELSE
                            CALL FlagError("Geometric field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(equationsSet%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for an advection-diffusion equation type of a classical field problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error)
    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh
  !   
  !================================================================================================================================
  !


  SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln(SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SolverType), POINTER :: SOLVER_ADVECTION_DIFFUSION !<A pointer to the solvers
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD_ADVECTION_DIFFUSION
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I

    ENTERS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%globalNumber==1) THEN
                !--- Get the dependent field of the advection-diffusion equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value of advection-diffusion & 
                  & (dependent field - U variable_type) at time, t ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%solvers,1,SOLVER_ADVECTION_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%solverEquations
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%solverMapping
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%equationsSets(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      DEPENDENT_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%dependent%dependentField
                      IF(ASSOCIATED(DEPENDENT_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL Field_NumberOfComponentsGet(DEPENDENT_FIELD_ADVECTION_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_ADVECTION_DIFFUSIONE is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Advection-diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Advection-diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Advection-diffusion solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_ADVECTION_DIFFUSION
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,err,error,*999)
                  END DO

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL Field_ParameterSetDataGet(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL Field_ParameterSetDataRestore(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error)
    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln    
  !
  !================================================================================================================================
  !
  SUBROUTINE AdvectionDiffusion_PreSolveGetSourceValue(SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SolverType), POINTER :: SOLVER_ADVECTION_DIFFUSION, SOLVER_DIFFUSION  !<A pointer to the solvers
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD_DIFFUSION, SOURCE_FIELD_ADVECTION_DIFFUSION
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS_ADVECTION_DIFFUSION, SOLVER_EQUATIONS_DIFFUSION  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING_ADVECTION_DIFFUSION, SOLVER_MAPPING_DIFFUSION !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET_ADVECTION_DIFFUSION, EQUATIONS_SET_DIFFUSION !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
    INTEGER(INTG) :: I


    ENTERS("AdvectionDiffusion_PreSolveGetSourceValue",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_ADVECTION_DIFFUSION)
      NULLIFY(SOLVER_DIFFUSION)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
               & PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%globalNumber==1) THEN
                !--- Get the dependent field of the diffusion equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Update advection-diffusion source field ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%solvers,2,SOLVER_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION=>SOLVER_DIFFUSION%solverEquations
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION)) THEN
                  SOLVER_MAPPING_DIFFUSION=>SOLVER_EQUATIONS_DIFFUSION%solverMapping
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION)) THEN
                    EQUATIONS_SET_DIFFUSION=>SOLVER_MAPPING_DIFFUSION%equationsSets(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION)) THEN
                      DEPENDENT_FIELD_DIFFUSION=>EQUATIONS_SET_DIFFUSION%dependent%dependentField
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION)) THEN
                        CALL Field_NumberOfComponentsGet(DEPENDENT_FIELD_DIFFUSION, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_DIFFUSION is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion solver equations are not associated.",err,error,*999)
                END IF


                !--- Get the source field for the advection-diffusion equations
                CALL SOLVERS_SOLVER_GET(SOLVER%solvers,1,SOLVER_ADVECTION_DIFFUSION,err,error,*999)
                SOLVER_EQUATIONS_ADVECTION_DIFFUSION=>SOLVER_ADVECTION_DIFFUSION%solverEquations
                IF(ASSOCIATED(SOLVER_EQUATIONS_ADVECTION_DIFFUSION)) THEN
                  SOLVER_MAPPING_ADVECTION_DIFFUSION=>SOLVER_EQUATIONS_ADVECTION_DIFFUSION%solverMapping
                  IF(ASSOCIATED(SOLVER_MAPPING_ADVECTION_DIFFUSION)) THEN
                    EQUATIONS_SET_ADVECTION_DIFFUSION=>SOLVER_MAPPING_ADVECTION_DIFFUSION%equationsSets(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_ADVECTION_DIFFUSION)) THEN
                      SOURCE_FIELD_ADVECTION_DIFFUSION=>EQUATIONS_SET_ADVECTION_DIFFUSION%SOURCE%sourceField
                      IF(ASSOCIATED(SOURCE_FIELD_ADVECTION_DIFFUSION)) THEN
                        CALL Field_NumberOfComponentsGet(SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION,err,error,*999)
                      ELSE
                        CALL FlagError("SOURCE_FIELD_ADVECTION_DIFFUSION is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Advection-diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Advection-diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Advection-diffusion solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the result from diffusion's dependent field to advection-diffusion's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_ADVECTION_DIFFUSION
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_ADVECTION_DIFFUSION, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,err,error,*999)
                  END DO
                ELSE
                  localError="Number of components of diffusion dependent field "// &
                    & "is not consistent with advection-diffusion equation source field."
                  CALL FlagError(localError,err,error,*999)
                END IF

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL Field_ParameterSetDataGet(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL Field_ParameterSetDataRestore(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("AdvectionDiffusion_PreSolveGetSourceValue")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveGetSourceValue",err,error)
    EXITS("AdvectionDiffusion_PreSolveGetSourceValue")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveGetSourceValue
  !   
  !================================================================================================================================
  !
  !>Update independent field (velocity) for advection-diffusion pre solve
  SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData(SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS

    INTEGER(INTG) :: numberOfDimensions
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)
    INTEGER(INTG) :: PROBLEM_SUBTYPE

    ENTERS("AdvectionDiffusion_PreSolveUpdateInputData",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          PROBLEM_SUBTYPE=CONTROL_LOOP%PROBLEM%SPECIFICATION(3)
          IF(PROBLEM_SUBTYPE==PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
            SOLVER_EQUATIONS=>SOLVER%solverEquations
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
              EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
              IF(ASSOCIATED(EQUATIONS)) THEN
                EQUATIONS_SET=>equations%equationsSet
                IF(ASSOCIATED(EQUATIONS_SET)) THEN
                  
                  CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDimensions,err,error,*999)
                  
                  INPUT_TYPE=1
                  INPUT_OPTION=1
                  NULLIFY(INPUT_DATA1)
                  CALL Field_ParameterSetDataGet(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_VALUES_SET_TYPE,INPUT_DATA1,err,error,*999)
                  CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, & 
                    & numberOfDimensions,INPUT_TYPE,INPUT_OPTION,CONTROL_LOOP%timeLoop%iterationNumber,1.0_DP, &
                    & err,error,*999)
                  
                ELSE
                  CALL FlagError("Equations set is not associated.",err,error,*999)
                END IF
                
              ELSE
                CALL FlagError("Equations are not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("Solver equations are not associated.",err,error,*999)
            END IF
          ELSE IF(PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE .OR. &
            & PROBLEM_SUBTYPE==PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE) THEN
            CALL FlagError("Not implemented.",err,error,*999)
          ENDIF
          
          CALL Field_ParameterSetUpdateStart(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveUpdateInputData",err,error)
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !
  !Update the boundary conditions
  SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC(SOLVER,err,error,*)
    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!\todo: Reduce number of variable used
    INTEGER(INTG) :: BOUNDARY_CONDITION_CHECK_VARIABLE,globalDOFIdx,localDOFIdx,nodeIdx,NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: 
    INTEGER(INTG), POINTER :: boundaryNodeS(:)
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
    TYPE(BoundaryConditionVariableType), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      WRITE (*,*) CONTROL_LOOP%timeLoop%iterationNumber
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE (PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,&
            & PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
          SOLVER_EQUATIONS=>SOLVER%solverEquations
          IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
            SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
            EQUATIONS=>SOLVER_MAPPING%equationsSetToSolverMatricesMap(1)%EQUATIONS
            IF(ASSOCIATED(EQUATIONS)) THEN
              EQUATIONS_SET=>equations%equationsSet
              IF(ASSOCIATED(EQUATIONS_SET)) THEN
                boundaryConditions=>SOLVER_EQUATIONS%boundaryConditions
                IF(ASSOCIATED(boundaryConditions)) THEN
                  DEPENDENT_FIELD=>equationsSet%dependent%dependentField
                  IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                    CALL Field_VariableGet(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VARIABLE,err,error,*999)
                    IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                      CALL BoundaryConditions_VariableGet(boundaryConditions,FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                        & err,error,*999)
                      IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                        CALL Field_NumberOfComponentsGet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                         & NUMBER_OF_COMPONENTS,err,error,*999)
                        NULLIFY(BOUNDARY_VALUES)
                        NULLIFY(boundaryNodes)
                        CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION(SOLVER_LINEAR_TYPE,BOUNDARY_VALUES, &
                          & boundaryNodes,NUMBER_OF_COMPONENTS,BOUNDARY_CONDITION_FIXED,CONTROL_LOOP%timeLoop%inputNumber, &
                          & CONTROL_LOOP%timeLoop%iterationNumber,err,error,*999)
                        WRITE(*,*) SIZE(BOUNDARY_VALUES)
                        DO nodeIdx=1,SIZE(BOUNDARY_VALUES)
                          !Default to version 1 of each node derivative
                          CALL Field_ComponentDOFGetUserNode(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                            & 1,NO_GLOBAL_DERIV,boundaryNodes(nodeIdx), &
                            & 1,localDOFIdx,globalDOFIdx,err,error,*999)
                          BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                            & conditionTypes(localDOFIdx)
                          IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                            CALL Field_ParameterSetUpdateLocalDOF(equationsSet%dependent%dependentField, &
                              & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                              & BOUNDARY_VALUES(nodeIdx),err,error,*999)
                          END IF
                        ENDDO
                      ELSE
                        CALL FlagError("Boundary condition variable is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Dependent field variable is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Equations set dependent variable is not associated.",err,error,*999)
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
          CALL Field_ParameterSetUpdateStart(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_VALUES_SET_TYPE,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for an advection-diffusion equation of a classical field problem class."
          CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC")
    RETURN
999 ERRORSEXITS("ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC",err,error)
    RETURN 1

  END SUBROUTINE ADVECTION_DIFFUSION_PRE_SOLVE_UPDATE_BC
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE AdvectionDiffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PostSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for an advection-diffusion problem.",err,error,*999)
 
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE)
      CALL ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(controlLoop,solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_NONLINEAR_SOURCE_ALE_ADVECTION_DIFFUSION_SUBTYPE, &
      & PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE,PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
      ! do nothing ???
      !CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "// &
        & TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PostSolve")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PostSolve
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(SolverEquationsType), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    ENTERS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_GENERALISED_ADVECTION_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%solverEquations
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%timeLoop%iterationNumber
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%timeLoop%outputNumber

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%timeLoop%currentTime<=CONTROL_LOOP%timeLoop%stopTime) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
!                        FILE="TRANSIENT_OUTPUT"
!                         METHOD="FORTRAN"
!                         EXPORT_FIELD=.TRUE.
!                         IF(EXPORT_FIELD) THEN          
!                          IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                          CALL FLUID_MECHANICS_IO_WRITE_CMGUI(equationsSet%REGION,equationsSet%globalNumber,FILE, &
                              & err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
!                           ENDIF
!                         ENDIF 
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_NONLINEAR_SOURCE_ADVECTION_DIFFUSION_SUBTYPE)
              ! do nothing ???
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for an advection-diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
    
  END SUBROUTINE ADVECTION_DIFFUSION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !
  
END MODULE AdvectionDiffusionEquationsRoutines

