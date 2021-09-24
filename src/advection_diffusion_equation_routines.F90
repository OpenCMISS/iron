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
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
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
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
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
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,independentVariable,materialsVariable,sourceVariable
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
              CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
                & dependentValue,err,error,*999)
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    
    !>Set the independent field (i.e. the advective velocity) to a specified analytical function
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    CALL Field_NumberOfVariablesGet(independentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(independentVariable)
      CALL Field_VariableIndexGet(independentField,variableIdx,independentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(independentVariable,numberOfComponents,err,error,*999)
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
      CALL FieldVariable_ParameterSetUpdateStart(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
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

  !>Sets up the advection-diffusion equation.
  SUBROUTINE AdvectionDiffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,equationsSetSubtype,esFieldNumberOfComponents,esFieldNumberOfVariables,esSpecification(3), &
      & geometricComponentNumber,geometricMeshComponent,geometricScalingType,lumpingType,myMatrixIdx, &
      & numberOfCompartments,numberOfDimensions,numberOfIndependentComponents,numberOfIndependentUComponents, &
      & numberOfIndependentVComponents,numberOfMaterialsComponents,numberOfSourceComponents,numberOfVariables, &
      & numberOfVariablesCount,numberOfMaterialsCouplingComponents,solutionMethod,sparsityType,variableIdx
    INTEGER(INTG), POINTER :: esFieldData(:)
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:),variableUTypes(:),couplingMatrixStorageTypes(:), &
      & couplingMatrixStructureTypes(:)
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetSourceType), POINTER :: equationsSource
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(FieldType), POINTER :: analyticField,dependentField,equationsSetField,geometricField
    TYPE(FieldVariableType), POINTER :: geometricVariable
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    equationsSetSubtype=esSpecification(3)
    SELECT CASE(equationsSetSubtype)
    CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " has not been implemented."
      CALL FlagError(localError,err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid for an advection-diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
      
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   s e t u p
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL AdvectionDiffusion_EquationsSetSolnMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          esFieldNumberOfVariables = 1
          esFieldNumberOfComponents = 2
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            !Create the auto created equations set field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField,err,error,*999)
            CALL Field_LabelSet(equationsField%equationsSetField,"Equations Set Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsField%equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsField%equationsSetField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsField%equationsSetField,esFieldNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
              & esFieldNumberOfComponents,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,esFieldNumberOfVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,esFieldNumberOfComponents, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 1,1_INTG,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & 2,1_INTG,err,error,*999)
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
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      SELECT CASE(equationsSetSubtype)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          esFieldNumberOfComponents = 2
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsField%equationsSetField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsField%equationsSetField,geometricField,err,error,*999)
            CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,1,geometricComponentNumber,err,error,*999)
            DO componentIdx=1,esFieldNumberOfComponents
              CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsField%equationsSetField,geometricScalingType,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          !Do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for an advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Do nothing???
      CASE(EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE)
        !ALE types
        CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetEnsureCreated(geometricVariable,FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
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
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
        IF(esSpecification(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
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
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
              & err,error,*999)
            !Default to the geometric interpolation setup
            CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,1,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
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
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !DO componentIdx=1,numberOfDimensions
              componentIdx=1
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE, & 
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              !ENDDO !componentIdx
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
        ELSE
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
            CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
              & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & err,error,*999)
            !Default to the geometric interpolation setup
            CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,1,geometricMeshComponent,err,error,*999)
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
              !Uses number of compartments to check that appropriate number and type of variables have been set on the
              !dependent field
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
              CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)                 
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
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
              SELECT CASE(solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !DO componentIdx=1,numberOfDimensions
                componentIdx=1
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
                  & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !ENDDO !componentIdx
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
        ENDIF
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
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
            CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,2,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            numberOfMaterialsComponents=numberOfDimensions
            !Set the number of materials components
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
            numberOfCompartments=esFieldData(2)
            numberOfMaterialsCouplingComponents=numberOfCompartments
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
              & numberOfMaterialsCouplingComponents,err,error,*999)
            !Default the k materials components to the geometric interpolation setup with constant interpolation
            DO componentIdx=1,numberOfDimensions
              CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,componentIdx,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfMaterialsCouplingComponents
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
          ELSE !standard materials field
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%REGION,equationsMaterials% &
              & materialsField,err,error,*999)
            CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              numberOfMaterialsComponents=numberOfDimensions
            ELSEIF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
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
              CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,componentIdx,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
              CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,1,geometricMeshComponent,err,error,*999)   
              DO componentIdx=numberOfDimensions+1,numberOfMaterialsComponents
                CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                  & geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
            ENDIF
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
          ENDIF
        ELSE
          !Not auto-created
          IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
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
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            ELSE IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions+1,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            numberOfMaterialsComponents=numberOfDimensions             
          ELSE IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
            !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
            numberOfMaterialsComponents=numberOfDimensions+1
          ELSE
            numberOfMaterialsComponents=numberOfDimensions
          ENDIF
          !First set the k values to 1.0
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
          IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
            numberOfCompartments=esFieldData(2)
            DO componentIdx=1,numberOfCompartments
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
          IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            !Now set the linear source values to 1.0
            DO componentIdx=numberOfDimensions+1,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsSource)
      CALL EquationsSet_SourceGet(equationsSet,equationsSource,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          !Create the auto created source field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSource%sourceField,err,error,*999)
          CALL Field_LabelSet(equationsSource%sourceField,"Source Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          numberOfSourceComponents=1
          !Set the number of source components
          CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,numberOfSourceComponents, &
            & err,error,*999)
          !Default the source components to the geometric interpolation setup with constant interpolation
          IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. &  
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            DO componentIdx=1,numberOfSourceComponents
              CALL FieldVariable_ComponentMeshComponentGet(geometricVariable,componentIdx,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
          ENDIF
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          !Finish creating the source field
          CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
          !Set the default values for the source field
          IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            numberOfSourceComponents=1
          ELSE
            numberOfSourceComponents=0
          ENDIF
          !Now set the source values to 1.0
          DO componentIdx=1,numberOfSourceComponents
            CALL Field_ComponentValuesInitialise(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsIndependent)
      CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      !Setup the equations set for the advective velocity field
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
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
            numberOfIndependentUComponents=numberOfDimensions
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfIndependentUComponents,err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
            numberOfCompartments=esFieldData(2)
            numberOfIndependentVComponents=numberOfCompartments-1
            !Set the number of independent components
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE, &
              & numberOfIndependentVComponents,err,error,*999)
            !Default the k independent components to the geometric interpolation setup with constant interpolation
            DO componentIdx=1,numberOfIndependentUComponents
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                & err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx=1,numberOfIndependentVComponents
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            NULLIFY(esFieldData)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
              & err,error,*999)
            numberOfCompartments=esFieldData(2)
            numberOfIndependentVComponents=numberOfCompartments-1
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfIndependentVComponents, &
              & err,error,*999)
            DO componentIdx=1,numberOfIndependentUComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx=1,numberOfIndependentVComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        ELSE
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
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            numberOfIndependentComponents=numberOfDimensions
            !Set the number of independent components
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfIndependentComponents,err,error,*999)
            !Default the k independent components to the geometric interpolation setup with constant interpolation
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
                & err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsIndependent%independentFieldAutoCreated) THEN
          !Finish creating the independent field
          CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
          !Set the default values for the independent field
          CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
            
          !IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
          !  & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          !  CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
          !    & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
          !ENDIF
          !numberOfIndependentComponents=numberOfDimensions             
          !!First set the k values to 1.0
          !DO componentIdx=1,numberOfIndependentComponents
          !  CALL Field_ComponentValuesInitialise(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
          !    & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
          !ENDDO !componentIdx
          
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(geometricVariable)
      CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        SELECT CASE(equationsSetSetup%analyticFunctionType)
        CASE(EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1)
          IF(numberOfDimensions/=2) THEN
            localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 2 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_ADVECTION_DIFFUSION_EQUATION_TWO_DIM_1
        CASE DEFAULT
          localError="The specified analytic function type of "// &
            & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
            & " is invalid for a linear advection-diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsAnalytic%analyticFieldAutoCreated) THEN
          !CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) 
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CASE(EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) 
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
        CASE DEFAULT
          localError="The equations set subtype of "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))//" is not valid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE .OR. & 
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN 
            !Finish the equations
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            CALL Equations_CreateFinish(equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            SELECT CASE(equationsSetSubtype)
            CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
              & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
              & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE,&
              & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE,&
              & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE,&
              & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE,&
              & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE,&
              & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE)
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
                & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE) THEN
                CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
                CALL EquationsMappingVector_SourcesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
              & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
              NULLIFY(equationsSetField)
              CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
              NULLIFY(esFieldData)
              CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,esFieldData, &
                & err,error,*999)
              myMatrixIdx = esFieldData(1)
              numberOfCompartments = esFieldData(2)    
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,numberOfCompartments-1,err,error,*999)
              
              ALLOCATE(variableTypes(2*numberOfCompartments))
              ALLOCATE(variableUTypes(numberOfCompartments-1))
              DO numberOfVariables=1,numberOfCompartments
                variableTypes(2*numberOfVariables-1)=FIELD_U_VARIABLE_TYPE+ &
                  & (FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(numberOfVariables-1))
                variableTypes(2*numberOfVariables)=FIELD_DELUDELN_VARIABLE_TYPE+ &
                  & (FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(numberOfVariables-1))
              ENDDO !numberOfVariables
              numberOfVariablesCount=0
              DO numberOfVariables=1,numberOfCompartments
                IF(numberOfVariables/=myMatrixIdx)THEN
                  numberOfVariablesCount=numberOfVariablesCount+1
                  variableUTypes(numberOfVariablesCount)=variableTypes(2*numberOfVariables-1)
                ENDIF
              ENDDO !numberOfVariables
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx-1),err,error,*999)
              CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,variableUTypes,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,variableTypes(2*myMatrixIdx),err,error,*999)
              CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_SourcesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            END SELECT
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            !Set up matrix storage and structure
            CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
            IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_UNLUMPED, &
                & EQUATIONS_MATRIX_LUMPED],err,error,*999)
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                & EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            ELSE
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                  & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)    
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
                  & TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
          ELSE IF(equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE .OR. &
            & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
            !Finish the creation of the equations
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            CALL Equations_CreateFinish(equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            NULLIFY(vectorMapping)
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
            IF(equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE .OR. &
              & equationsSetSubtype==EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
              CALL EquationsMappingVector_SourcesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            ENDIF
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
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
    CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE,&
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE,&
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE,&
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, & 
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
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
    CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
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

  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE AdvectionDiffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,componentIdx, &
      & equationsSetSubtype,esSpecification(3),gaussPointIdx,matrixIdx,myCompartment,numberOfColumnElementParameters, &
      & numberOfColsComponents,numberOfCompartments,numberOfCouplingComponents,numberOfDimensions,numberOfGauss, &
      & numberOfRowElementParameters,numberOfRowsComponents,numberOfVariablesCount,numberOfXi,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx, rowsVariableType,rowXiIdx,scalingType,sourceVariableType,xiIdx
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COUPLING_MATRICES=99
    INTEGER(INTG) :: couplingVariableTypes(MAX_NUMBER_OF_COUPLING_MATRICES)
    INTEGER(INTG), POINTER :: esVariableData(:)
    REAL(DP) :: aParam,bParam,columnPhi,columndPhidXi(3),columndPhidX(3),conductivity(3,3),cParam,coupledSourceParam, &
      & couplingParam,dXidX(3,3),gaussWeight,jacobian,jacobianGaussWeight,rowPhi,rowdPhidXi(3),sourceParam,sum,w(3),v(3)
    LOGICAL :: update,updateCoupling,updateCouplingMatrices(MAX_NUMBER_OF_COUPLING_MATRICES),updateDamping,updateMatrices, &
      & updateStiffness,updateRHS,updateSource
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsMatrixPtrType) :: couplingMatrices(MAX_NUMBER_OF_COUPLING_MATRICES) 
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,fibreField,geometricField,independentField,materialsField, &
      & sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,couplingInterpParameters,fibreInterpParameters, &
      & geometricInterpParameters,independentInterpParameters,rowsInterpParameters,sourceInterpParameters, &
      & uDependentInterpParameters,uMaterialsInterpParameters,vDependentInterpParameters,vMaterialsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint,geometricInterpPoint,independentInterpPoint,sourceInterpPoint, &
      & uDependentInterpPoint,uMaterialsInterpPoint,vDependentInterpPoint,vMaterialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,equationsSetVariable,fibreVariable,geometricVariable, &
      & independentVariable,rowsVariable,sourceVariable,uMaterialsVariable,vDependentVariable,vMaterialsVariable
    TYPE(FieldVariablePtrType) :: couplingVariables(MAX_NUMBER_OF_COUPLING_MATRICES)
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    equationsSetSubtype=esSpecification(3)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(equationsInterpolation)
    CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(equationsSetField)
    NULLIFY(equationsSetVariable)
    NULLIFY(esVariableData)
    myCompartment=0
    numberOfCompartments=0
    NULLIFY(dynamicMapping)
    NULLIFY(dynamicMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dampingMatrix)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(colsVariable)
    updateStiffness=.FALSE.
    updateDamping=.FALSE.
    updateCoupling=.FALSE.
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_StiffnessMatrixGet(dynamicMatrices,stiffnessMatrix,err,error,*999)
      CALL EquationsMatricesDynamic_DampingMatrixGet(dynamicMatrices,dampingMatrix,err,error,*999)      
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
      CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_StiffnessMatrixGet(dynamicMatrices,stiffnessMatrix,err,error,*999)
      CALL EquationsMatricesDynamic_DampingMatrixGet(dynamicMatrices,dampingMatrix,err,error,*999)      
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      CALL Field_VariableGet(equationsSetField,FIELD_U_VARIABLE_TYPE,equationsSetVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(equationsSetVariable,FIELD_VALUES_SET_TYPE,esVariableData,err,error,*999)
      myCompartment=esVariableData(1)
      numberOfCompartments=esVariableData(2)
      IF(numberOfCompartments>MAX_NUMBER_OF_COUPLING_MATRICES) THEN
        localError="Insufficient number of coupling matrices. The maximum number of allowed coupling matrices is "// &
          & TRIM(NumberToVString(MAX_NUMBER_OF_COUPLING_MATRICES,"*",err,error))//" and "// &
          & TRIM(NumberToVString(numberOfCompartments,"*",err,error))//" are required."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      numberOfVariablesCount=0
      DO matrixIdx=1,numberOfCompartments
        IF(matrixIdx/=myCompartment) THEN
          numberOfVariablesCount=numberOfVariablesCount+1
          NULLIFY(couplingMatrices(numberOfVariablesCount)%ptr)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,couplingMatrices(numberOfVariablesCount)%ptr, &
            & err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(couplingMatrices(numberOfVariablesCount)%ptr, &
            & updateCouplingMatrices(numberOfVariablesCount),err,error,*999)
          updateCoupling=(updateCoupling.OR.updateCouplingMatrices(numberOfVariablesCount))
        ENDIF
      ENDDO !matrixIdx
    CASE(EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_STATIC_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUPG_SUBTYPE)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)      
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
      & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE)
      CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    NULLIFY(sourceVariable)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
      CALL EquationsMappingSource_SourceVariableGet(sourceMapping,sourceVariable,err,error,*999)
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateCoupling.OR.updateRHS.OR.updateSource)

    IF(update) THEN

      IF(ASSOCIATED(equationsSetVariable)) THEN
        numberOfVariablesCount=0
        DO matrixIdx=1,numberOfCompartments
          IF(matrixIdx/=myCompartment)THEN
            numberOfVariablesCount=numberOfVariablesCount+1
            NULLIFY(couplingVariables(numberOfVariablesCount)%ptr)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,matrixIdx, &
              & couplingVariables(numberOfVariablesCount)%ptr,err,error,*999)
            CALL EquationsMappingLinear_LinearMatrixVariableTypeGet(linearMapping,matrixIdx, &
              & couplingVariableTypes(numberOfVariablesCount),err,error,*999)
            couplingMatrices(numberOfVariablesCount)%ptr%elementMatrix%matrix=0.0_DP
          ENDIF
        ENDDO !matrixIdx
      ENDIF

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
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
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,geometricInterpPoint, &
        & err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)

      NULLIFY(fibreField)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      NULLIFY(fibreVariable)
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      IF(ASSOCIATED(fibreField)) THEN
        CALL Field_VariableGet(fibreField,FIELD_U_VARIABLE_TYPE,fibreVariable,err,error,*999)
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF

      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      NULLIFY(vDependentVariable)
      IF(equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
        CALL Field_VariableGet(dependentField,FIELD_V_VARIABLE_TYPE,vDependentVariable,err,error,*999)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & vDependentInterpParameters,err,error,*999)
        NULLIFY(uDependentInterpPoint)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vDependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vDependentInterpParameters,err,error,*999)
      ENDIF

      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsVariable)
      NULLIFY(vMaterialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,uMaterialsVariable,err,error,*999)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & uMaterialsInterpParameters,err,error,*999)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uMaterialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,uMaterialsInterpParameters,err,error,*999)
      IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
        CALL Field_VariableGet(materialsField,FIELD_V_VARIABLE_TYPE,vMaterialsVariable,err,error,*999)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & vMaterialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vMaterialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vMaterialsInterpParameters,err,error,*999)
      ENDIF

      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(independentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
      NULLIFY(independentInterpParameters)
      CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & independentInterpParameters,err,error,*999)
      NULLIFY(independentInterpPoint)
      CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,independentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters,err,error,*999)

      NULLIFY(sourceVariable)
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      IF(ASSOCIATED(sourceVariable)) THEN
        CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
        CALL FieldVariable_VariableTypeGet(sourceVariable,sourceVariableType,err,error,*999)
        NULLIFY(sourceInterpParameters)
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,sourceVariableType,sourceInterpParameters, &
          & err,error,*999)
        NULLIFY(sourceInterpPoint)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,sourceVariableType,sourceInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
      ENDIF

      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)

      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)     

      !Select whether using standard Galerkin scheme, or the stabilised streamwise-upwinding Petrov-Galerkin scheme

!!TODO: NEED TO ADD IN SUPG

      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss

        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        IF(ASSOCIATED(fibreVariable)) &
          & CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,uMaterialsInterpPoint, &
          & err,error,*999)
        !Get coupling source parameter
        IF(ASSOCIATED(vDependentVariable)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,vDependentInterpPoint, &
            & err,error,*999)
          coupledSourceParam=vDependentInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        !Interpolate to get the advective velocity
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
          & err,error,*999)
        IF(ASSOCIATED(sourceVariable)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
            & err,error,*999)            
          sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        IF(ASSOCIATED(vMaterialsVariable)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,vMaterialsInterpPoint, &
            & err,error,*999)
          couplingParam=vMaterialsInterpPoint%values(myCompartment,NO_PART_DERIV)
        ENDIF

        aParam=uMaterialsInterpPoint%values(1,NO_PART_DERIV)
        bParam=uMaterialsInterpPoint%values(2,NO_PART_DERIV)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE , & 
          & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
          cParam=uMaterialsInterpPoint%values(3,NO_PART_DERIV)
          !Calculate conductivity tensor
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & uMaterialsInterpPoint%values(4:4+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)
        CASE DEFAULT
          !Calculate conductivity tensor
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & uMaterialsInterpPoint%values(3:3+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)  
        END SELECT

        !Get the advection velocity
        v(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)

        !Get the ALE mesh velocity
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE,& 
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
          w(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        CASE DEFAULT
          w=0.0_DP
        END SELECT

        !Calculate Jacobian and Gauss weight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        DO columnComponentIdx=1,numberOfColsComponents
          DO rowXiIdx=1,numberOfXi
            dXidX(rowXiIdx,columnComponentIdx)=geometricInterpPointMetrics%dXidX(rowXiIdx,columnComponentIdx)
          ENDDO !rowXiIdx
        ENDDO !columnXiIdx

        !Loop over row components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowPhi,err,error,*999)
            DO rowXiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,rowdPhidXi(rowXiIdx),err,error,*999)
            ENDDO !rowXiIIdx
            IF(updateMatrices) THEN
              columnElementDOFIdx=0
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                  & err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                    & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  IF(updateStiffness) THEN
                    columndPhidX=0.0_DP
                    DO columnXiIdx=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,columndPhidXi(columnXiIdx), &
                        & err,error,*999)
                      DO componentIdx=1,numberOfDimensions
                        columndPhidX(componentIdx)=columndPhidX(componentIdx)+ &
                          & columndPhidXi(columnXiIdx)*dXidX(columnXiIdx,componentIdx)
                      ENDDO !componentIdx
                    ENDDO !columnXiiIdx
                    sum=0.0_DP
                    DO rowXiIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    !Add in the advection term
                    DO componentIdx=1,numberOfDimensions
                      sum=sum+bParam*v(componentIdx)*rowPhi*columndPhidX(componentIdx)
                    ENDDO !componentIdx
                    !Add in the linear source term
                    SELECT CASE(equationsSetSubtype)
                    CASE(EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
                      & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
                      & EQUATIONS_SET_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
                      sum=cParam*rowPhi*columnPhi
                    END SELECT
                    IF(equationsSetSubtype==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE) THEN
                      !For multi-compartment model must include additional terms into the
                      sum=sum+couplingParam*rowPhi*columnPhi
                    ENDIF
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                  IF(updateDamping) THEN
                    sum=aParam*rowPhi*columnPhi
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateCoupling) THEN
              numberOfVariablesCount=0
              DO matrixIdx=1,numberOfCompartments
                IF(matrixIdx/=myCompartment) THEN
                  numberOfVariablesCount=numberOfVariablesCount+1
                  couplingParam=vMaterialsInterpPoint%values(matrixIdx,NO_PART_DERIV)
                  !Need to test for the case where imatrix==mycompartment the coupling terms then needs to be added into the
                  !stiffness matrix                  
                  IF(updateCouplingMatrices(numberOfVariablesCount)) THEN
                    CALL FieldVariable_NumberOfComponentsGet(couplingVariables(numberOfVariablesCount)%ptr, &
                      & numberOfCouplingComponents,err,error,*999)
                    !Loop over element columns
                    columnElementDOFIdx=0
                    DO columnComponentIdx=1,numberOfCouplingComponents
                      NULLIFY(columnDomain)
                      CALL FieldVariable_ComponentDomainGet(couplingVariables(numberOfVariablesCount)%ptr,columnComponentIdx, &
                        & columnDomain,err,error,*999)
                      NULLIFY(columnDomainTopology)
                      CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                      NULLIFY(columnDomainElements)
                      CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                      NULLIFY(columnBasis)
                      CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                      CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme, &
                        & err,error,*999)
                      CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                          & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                        
                        !-------------------------------------------------------------------------------------------------
                        !concentration test function, concentration trial function
                        !For now, this is only a dummy implementation - this still has to be properly set up.
                        
                        couplingMatrices(numberOfVariablesCount)%ptr%elementMatrix% &
                          & matrix(rowElementDOFIdx,columnElementDOFIdx) = couplingMatrices(numberOfVariablesCount)%ptr% &
                          & elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + & 
                          & couplingParam*rowPhi*columnPhi*jacobianGaussWeight
                        
                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                  ENDIF !updateCouplingMatrix
                ENDIF !myCompartment
              ENDDO !matrixIdx
            ENDIF !updateCoupling
            IF(updateSource) THEN
              IF(equationsSetSubtype==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
                sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                  & coupledSourceParam*rowPhi*jacobianGaussWeight
              ELSE
                sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                  & sourceParam*rowPhi*jacobianGaussWeight
              ENDIF
            ENDIF !updateSource
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
            ENDIF !updateRHS
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDDO !gausssPointIdx
      
      !Scale factor adjustment
      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        NULLIFY(rowsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsInterpParameters, &
          & err,error,*999)
        NULLIFY(colsInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
          & err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,colsInterpParameters,err,error,*999)
        
        !Loop over element rows
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  IF(updateStiffness) THEN
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update stiffness
                  IF(updateDamping) THEN
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF !update damping
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrices
            IF(updateCoupling) THEN
              numberOfVariablesCount=0
              DO matrixIdx=1,numberOfCompartments
                IF(matrixIdx/=myCompartment) THEN
                  numberOfVariablesCount=numberOfVariablesCount+1
                  couplingParam=vMaterialsInterpPoint%values(matrixIdx,NO_PART_DERIV)
                  !Need to test for the case where imatrix==mycompartment the coupling terms then needs to be added into the
                  !stiffness matrix                  
                  IF(updateCouplingMatrices(numberOfVariablesCount)) THEN
                    NULLIFY(couplingInterpParameters)
                    CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation, &
                      & couplingVariableTypes(numberOfVariablesCount),couplingInterpParameters, &
                      & err,error,*999)
                    CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,couplingInterpParameters,err,error,*999)
                    CALL FieldVariable_NumberOfComponentsGet(couplingVariables(numberOfVariablesCount)%ptr, &
                      & numberOfCouplingComponents,err,error,*999)
                    !Loop over element columns
                    columnElementDOFIdx=0
                    DO columnComponentIdx=1,numberOfCouplingComponents
                      NULLIFY(columnDomain)
                      CALL FieldVariable_ComponentDomainGet(couplingVariables(numberOfVariablesCount)%ptr,columnComponentIdx, &
                        & columnDomain,err,error,*999)
                      NULLIFY(columnDomainTopology)
                      CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                      NULLIFY(columnDomainElements)
                      CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                      NULLIFY(columnBasis)
                      CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                      CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                      DO columnElementParameterIdx=1,numberOfColumnElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        couplingMatrices(numberOfVariablesCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                          & couplingMatrices(numberOfVariablesCount)%ptr% &
                          & elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* & 
                          & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                          & couplingInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                  ENDIF !updateCouplingMatrix
                ENDIF !myCompartment
              ENDDO !matrixIdx
            ENDIF
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update
    
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
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for an advection-diffusion type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_ADVECTION_DIFFUSION_EQUATION_TYPE,problemSubtype]

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
  SUBROUTINE AdvectionDiffusion_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables    
    INTEGER(INTG) :: problemSubtype,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("AdvectionDiffusion_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    problemSubtype=pSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(problemSubtype,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      !
      ! Initial Setup
      !
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
      !
      ! Control loop setup
      !
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        IF(problemSubtype==PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE.OR. &
          & problemSubtype==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
          !Set up a simple control loop
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_SIMPLE_TYPE,err,error,*999)
        ELSE
          !Set up a time control loop
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
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
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !
      ! Solvers setup
      !
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
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        IF(problemSubtype==PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE.OR. &
          & problemSubtype==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
          !Start the linear solver creation
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        ELSE
          !Set the solver to be a first order dynamic solver 
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)        
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
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
          & " is invalid for a linear advection-diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !
      ! Solver equations setup
      !
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      !Get the control loop
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
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        IF(problemSubtype==PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE.OR. &
          & problemSubtype==PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE) THEN
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        ELSE
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        ENDIF
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
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
       
    EXITS("AdvectionDiffusion_ProblemSetup")
    RETURN
999 ERRORS("AdvectionDiffusion_ProblemSetup",err,error)
    EXITS("AdvectionDiffusion_ProblemSetup")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Performs pre-solve operations for advection-diffusion problems.
  SUBROUTINE AdvectionDiffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: outputType,problemSubType,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solver2
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    
    problemSubType=problem%specification(3)
    SELECT CASE(problemSubType)
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE)
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
      !Update independent data fields
      CALL AdvectionDiffusion_PreSolveUpdateInputData(solver,err,error,*999)
      !CALL AdvectionDiffusion_PreSolveUpdateBC(solver,err,error,*999)
    CASE(PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE) 
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
      !Update independent data fields
      CALL AdvectionDiffusion_PreSolveUpdateInputData(solver,err,error,*999)
      !CALL AdvectionDiffusion_PreSolveUpdateBC(solver,err,error,*999)
      CALL AdvectionDiffusion_PreSolveALEUpdateMesh(solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)         
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
      !do nothing
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(problemSubType,"*",err,error))// &
        & " is not valid for a advection-diffusion type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PreSolve")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolve
  
  !   
  !================================================================================================================================
  !
  
  !>Update mesh position and velocity for ALE advection-diffusion problem
  SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,dofIdx,dofNumber,esSpecification(3),inputIteration,inputOption,inputType,numberOfDimensions, &
      &numberDOFsToPrint,outputIteration,outputTYpe,pSpecification(3),totalNumberOfDofs
    REAL(DP) :: alpha,currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: inputData1(:),meshDisplacementValues(:)
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: geometricField
    TYPE(FieldVariableType), POINTER :: geometricVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverALEDiffusion
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
      
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_EXP_SOURCE_ADVEC_DIFF_SUBTYPE)
        !do nothing ???
      CASE(EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE, &
        & EQUATIONS_SET_EXP_SOURCE_ALE_ADVEC_DIFF_SUPG_SUBTYPE)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Advection-diffusion update mesh ...",err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        !--- First, read mesh displacement values from file        
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)

        inputType=42
        inputOption=2
        NULLIFY(inputData1)
        !CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,inputData1,err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputData1,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)

        NULLIFY(meshDisplacementValues)
        CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,meshDisplacementValues, &
          & err,error,*999)
        
        IF(diagnostics1) THEN
          numberDOFsToPrint = SIZE(meshDisplacementValues,1)
          CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,&
            & meshDisplacementValues,'(" meshDisplacementValues = ",3(X,E13.6))','3(3(X,E13.6))', err,error,*999)
        ENDIF

        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputData1,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)

        CALL FieldVariable_TotalNumberOfDOFsGet(geometricVariable,totalNumberOfDOFs,err,error,*999)

        !--- Second, update geometric field
        DO dofIdx=1,totalNumberOfDofs
          CALL FieldVariable_ParameterSetAddLocalDOF(geometricVariable,FIELD_VALUES_SET_TYPE,dofIdx, &
            & meshDisplacementValues(dofIdx),err,error,*999)
        ENDDO !dofIdx
        CALL FieldVariable_ParameterSetUpdateStart(geometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(geometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)

        !--- Third, use displacement values to calculate velocity values
        alpha=1.0_DP/timeIncrement
        CALL FieldVariable_ParameterSetsCopy(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
          & alpha,err,error,*999)
        CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,meshDisplacementValues, &
          & err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for an advection-diffusion equation type of a classical field problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveALEUpdateMesh",err,error)
    EXITS("AdvectionDiffusion_PreSolveALEUpdateMesh")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveALEUpdateMesh
  
  !   
  !================================================================================================================================
  !

  SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln(solver,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,numberOfComponents,outputType,pSpecification(3),solverGlobalNumber
    TYPE(EquationsSetType), POINTER :: advectionDiffusionEquationsSet
    TYPE(FieldType), POINTER :: advectionDiffusionDependentField
    TYPE(FieldVariableType), POINTER :: advectionDiffusionDependentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: advectionDiffusionSolver
    TYPE(SolverEquationsType), POINTER :: advectionDiffusionSolverEquations
    TYPE(SolverMappingType), POINTER :: advectionDiffusonSolverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==1) THEN
        !--- Get the dependent field of the advection-diffusion equations
        CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value of advection-diffusion dependent field at time, t ... ", &
          & err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(advectionDiffusionSolver)
        CALL Solvers_SolverGet(solvers,1,advectionDiffusionSolver,err,error,*999)
        NULLIFY(advectionDiffusionSolverEquations)
        CALL Solver_SolverEquationsGet(advectionDiffusionSolver,advectionDiffusionSolverEquations,err,error,*999)
        NULLIFY(advectionDiffusonSolverMapping)
        CALL SolverEquations_SolverMappingGet(advectionDiffusionSolverEquations,advectionDiffusonSolverMapping,err,error,*999)
        NULLIFY(advectionDiffusionEquationsSet)
        CALL SolverMapping_EquationsSetGet(advectionDiffusonSolverMapping,1,advectionDiffusionEquationsSet,err,error,*999)
        NULLIFY(advectionDiffusionDependentField)
        CALL EquationsSet_DependentFieldGet(advectionDiffusionEquationsSet,advectionDiffusionDependentField,err,error,*999)
        NULLIFY(advectionDiffusionDependentVariable)
        CALL Field_VariableGet(advectionDiffusionDependentField,FIELD_U_VARIABLE_TYPE,advectionDiffusionDependentVariable, &
          & err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(advectionDiffusionDependentVariable,numberOfComponents,err,error,*999)

        !--- Copy the current time value parameters set from diffusion-one's dependent field 
        DO componentIdx=1,numberOfComponents
          CALL FieldVariable_ParametersToFieldVariableParametersCopy(advectionDiffusionDependentVariable, & 
            & FIELD_VALUES_SET_TYPE,componentIdx,advectionDiffusionDependentVariable,FIELD_PREVIOUS_VALUES_SET_TYPE, &
            & componentIdx,err,error,*999)
        ENDDO !componentIdx

        !IF(diagnostics3) THEN
        !  NULLIFY(dummyValues2)
        !  CALL FieldVariable_ParameterSetDataGet(finiteElasticityDependentVariable,FIELD_VALUES_SET_TYPE,dummyValues2, &
        !    & err,error,*999)
        !  numberDOFsToPrint = SIZE(dummyValues2,1)
        !  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
        !    & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
        !    & '4(4(X,E13.6))',err,error,*999)
        !  CALL FieldVariable_ParameterSetDataRestore(finiteElasticityDependentVariable,FIELD_VALUES_SET_TYPE,dummyValuess2, &
        !    & err,error,*999)
        !ENDIF

      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveStoreCurrentSoln",err,error)
    EXITS("AdvectionDiffusion_PreSolveStoreCurrentSoln")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveStoreCurrentSoln
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE AdvectionDiffusion_PreSolveGetSourceValue(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,numberOfDependentComponents,numberOfSourceComponents,outputType,pSpecification(3), &
      & solverGlobalNumber
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsSetType), POINTER :: advectionDiffusionEquationsSet,diffusionEquationsSet
    TYPE(FieldType), POINTER :: diffusionDependentField,advectionDiffusionSourceField
    TYPE(FieldVariableType), POINTER :: diffusionDependentVariable,advectionDiffusionSourceVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: advectionDiffusionSolver,diffusionSolver
    TYPE(SolverEquationsType), POINTER :: advectionDiffusionSolverEquations,diffusionSolverEquations
    TYPE(SolverMappingType), POINTER :: advectionDiffusonSolverMapping,diffusionSolverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolveGetSourceValue",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(SolverGlobalNumber==1) THEN
        !--- Get the dependent field of the diffusion equations
        CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Update advection-diffusion source field ... ",err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(diffusionSolver)
        CALL Solvers_SolverGet(solvers,2,diffusionSolver,err,error,*999)
        NULLIFY(diffusionSolverEquations)
        CALL Solver_SolverEquationsGet(diffusionSolver,diffusionSolverEquations,err,error,*999)
        NULLIFY(diffusionSolverMapping)
        CALL SolverEquations_SolverMappingGet(diffusionSolverEquations,diffusionSolverMapping,err,error,*999)
        NULLIFY(diffusionEquationsSet)
        CALL SolverMapping_EquationsSetGet(diffusionSolverMapping,1,diffusionEquationsSet,err,error,*999)
        NULLIFY(diffusionDependentField)
        CALL EquationsSet_DependentFieldGet(diffusionEquationsSet,diffusionDependentField,err,error,*999)
        NULLIFY(diffusionDependentVariable)
        CALL Field_VariableGet(diffusionDependentField,FIELD_U_VARIABLE_TYPE,diffusionDependentVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(diffusionDependentVariable,numberOfDependentComponents,err,error,*999)
       
        !--- Get the source field for the advection-diffusion equations
        NULLIFY(advectionDiffusionSolver)
        CALL Solvers_SolverGet(solvers,1,advectionDiffusionSolver,err,error,*999)
        NULLIFY(advectionDiffusionSolverEquations)
        CALL Solver_SolverEquationsGet(advectionDiffusionSolver,advectionDiffusionSolverEquations,err,error,*999)
        NULLIFY(advectionDiffusonSolverMapping)
        CALL SolverEquations_SolverMappingGet(advectionDiffusionSolverEquations,advectionDiffusonSolverMapping,err,error,*999)
        NULLIFY(advectionDiffusionEquationsSet)
        CALL SolverMapping_EquationsSetGet(advectionDiffusonSolverMapping,1,advectionDiffusionEquationsSet,err,error,*999)
        NULLIFY(advectionDiffusionSourceField)
        CALL EquationsSet_SourceFieldGet(advectionDiffusionEquationsSet,advectionDiffusionSourceField,err,error,*999)
        NULLIFY(advectionDiffusionSourceVariable)
        CALL Field_VariableGet(advectionDiffusionSourceField,FIELD_U_VARIABLE_TYPE,advectionDiffusionSourceVariable,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(advectionDiffusionSourceVariable,numberOfSourceComponents,err,error,*999)

        !--- Copy the result from diffusion's dependent field to advection-diffusion's source field
        IF(numberOfSourceComponents/=numberOfDependentComponents) THEN
          localError="The number of advection diffusion source field components of "// &
            & TRIM(NumberToVString(numberOfSourceComponents,"*",err,error))// &
            & " does not match the number of diffusion dependent field components of "// &
            & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
        ENDIF
        
        DO componentIdx=1,numberOfSourceComponents
          CALL FieldVariable_ParametersToFieldVariableParametersCopy(diffusionDependentVariable,FIELD_VALUES_SET_TYPE, &
            & componentIdx,advectionDiffusionSourceVariable,FIELD_VALUES_SET_TYPE,componentIdx,err,error,*999)
        ENDDO !componentIdx
       
        !IF(diagnostics3) THEN
        !  NULLIFY(dummyValues2)
        !  CALL FieldVariable_ParameterSetDataGet(finiteElasticityDependentVariable,FIELD_VALUES_SET_TYPE,dummyValues2, &
        !    & err,error,*999)
        !  numberDOFsToPrint=SIZE(dummyValues2,1)
        !  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDOFsToPrint,numberDOFsToPrint,numberDOFsToPrint,dummyValues2, &
        !    & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
        !    & '4(4(X,E13.6))',err,error,*999)
        !  CALL FieldVariable_ParameterSetDataRestore(finiteElasticityDependentVariable,FIELD_VALUES_SET_TYPE,dummyValues2, &
        !    & err,error,*999)
        !ENDIF

      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

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
  SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,inputIteration,inputOption,inputType,numberOfDimensions,outputIteration,outputType, &
      & problemSubtype,pSpecification(3)
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: inputData1(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: geometricField,independentField
    TYPE(FieldVariableType), POINTER :: geometricVariable,independentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolveUpdateInputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
        & CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      NULLIFY(independentVariable)
      CALL Field_VariableGet(independentField,FIELD_U_VARIABLE_TYPE,independentVariable,err,error,*999)
            
      inputType=1
      inputOption=1
      NULLIFY(inputData1)
      CALL FieldVariable_ParameterSetDataGet(independentVariable,FIELD_VALUES_SET_TYPE,inputData1,err,error,*999)
      
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputData1,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      
      CALL FieldVariable_ParameterSetDataRestore(independentVariable,FIELD_VALUES_SET_TYPE,inputData1,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(independentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN
999 ERRORS("AdvectionDiffusion_PreSolveUpdateInputData",err,error)
    EXITS("AdvectionDiffusion_PreSolveUpdateInputData")
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !
  
  !>Update the boundary conditions
  SUBROUTINE AdvectionDiffusion_PreSolveUpdateBC(solver,err,error,*)
    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!\todo: Reduce number of variable used
    INTEGER(INTG) :: boundaryConditionCheckVariable,currentIteration,inputIteration,globalDOFIdx,localDOFIdx,nodeIdx, &
      & numberOfComponents,outputIteration,pSpecification(3)
    INTEGER(INTG), POINTER :: boundaryNodes(:)    
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: boundaryValues(:)
    LOGICAL :: ghostDOF
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PreSolveUpdateBC",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE (PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(boundaryConditions)
      CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(dependentVariable)
      CALL Field_VariableGet(dependentField,FIELD_U_VARIABLE_TYPE,dependentVariable,err,error,*999)
      CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      NULLIFY(boundaryValues)
      NULLIFY(boundaryNodes)
      CALL FLUID_MECHANICS_IO_READ_BOUNDARY_CONDITIONS_ITERATION(SOLVER_LINEAR_TYPE,boundaryValues,boundaryNodes,&
        & numberOfComponents,BOUNDARY_CONDITION_FIXED,inputIteration,currentIteration,err,error,*999)
      DO nodeIdx=1,SIZE(boundaryValues)
        !Default to version 1 of each node derivative
        CALL FieldVariable_ComponentDOFGetUserNode(dependentVariable,1,NO_GLOBAL_DERIV,boundaryNodes(nodeIdx),1,localDOFIdx, &
          & globalDOFIdx,err,error,*999)
        CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDOFIdx, &
          & boundaryConditionCheckVariable,err,error,*999)
        IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) &
          & CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
          & boundaryValues(nodeIdx),err,error,*999)
      ENDDO !nodeIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("AdvectionDiffusion_PreSolveUpdateBC")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PreSolveUpdateBC",err,error)
    RETURN 1

  END SUBROUTINE AdvectionDiffusion_PreSolveUpdateBC
  
  !
  !================================================================================================================================
  !
  
  SUBROUTINE AdvectionDiffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PostSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE)
      CALL AdvectionDiffusion_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_GENERALISED_STATIC_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_STATIC_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
      !CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  
  SUBROUTINE AdvectionDiffusion_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: currentIteration,equationsSetGlobalNumber,equationsSetIdx,inputIteration,numberOfEquationsSets, &
      & outputIteration,outputType,pSpecification(3)
    
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    CHARACTER(14) :: file,outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("AdvectionDiffusion_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_GENERALISED_ADVEC_DIFF_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_ADVEC_DIFF_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      !Make sure the equations sets are up to date
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        CALL EquationsSet_GlobalNumberGet(equationsSet,equationsSetGlobalNumber,err,error,*999)
        
        IF(outputIteration/=0) THEN
          IF(currentTime<=stopTime) THEN
            IF(currentIteration<10) THEN
              WRITE(outputFile,'("TimeStep_000",I0)') currentIteration
            ELSE IF(currentIteration<100) THEN
              WRITE(outputFile,'("TimeStep_00",I0)') currentIteration
            ELSE IF(currentIteration<1000) THEN
              WRITE(outputFile,'("TimeStep_0",I0)') currentIteration
            ELSE IF(currentIteration<10000) THEN
              WRITE(outputFile,'("TimeStep_",I0)') currentIteration
            END IF
            file=outputFile
            !file="TRANSIENT_OUTPUT"
            !method="FORTRAN"
            !exportField=.TRUE.
            !IF(exportField) THEN          
            IF(MOD(currentIteration,outputIteration)==0)  THEN   
              IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
              ENDIF
              CALL FLUID_MECHANICS_IO_WRITE_CMGUI(region,equationsSetGlobalNumber,file,err,error,*999)
              IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
              ENDIF
            ENDIF
            !ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE(PROBLEM_NONLINEAR_SOURCE_ADVEC_DIFF_SUBTYPE)
      !Do nothing ???
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for an advection-diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("AdvectionDiffusion_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("AdvectionDiffusion_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE AdvectionDiffusion_PostSolveOutputData

  !
  !================================================================================================================================
  !
  
END MODULE AdvectionDiffusionEquationsRoutines

