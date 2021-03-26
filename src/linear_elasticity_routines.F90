!> \file
!> \author Chris Bradley
!> \brief This module handles all linear elasticity routines.
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

!>This module handles all linear elasticity routines.
MODULE LinearElasticityRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
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
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE ProblemAccessRoutines
  USE Strings
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC LinearElasticity_BoundaryConditionsAnalyticCalculate
  
  PUBLIC LinearElasticity_EquationsSetSetup
  
  PUBLIC LinearElasticity_EquationsSetSolutionMethodSet
  
  PUBLIC LinearElasticity_EquationsSetSpecificationSet
  
  PUBLIC LinearElasticity_FiniteElementCalculate

  PUBLIC LinearElasticity_ProblemSpecificationSet
  
  PUBLIC LinearElasticity_ProblemSetup

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE LinearElasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the analytic boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<The boundary conditions to set the analytic boundary conditions for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,bcXCounter,bcXForceCounter,bcXNodes,bcYCounter,bcZNodes,componentIdx,derivativeIdx, &
      & dimensionIdx,globalDerivativeIndex,localDOFIdx,nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodes, &
      & numberOfNodeDerivatives,numberOfVariables,variableIdx,variableType
    REAL(DP) :: analyticValue,bcValue,E,geometricTolerance,forceX,forceXArea,forceY,forceYArea,forceZ,height,Iyy,length,width, &
      & vX,x(3)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: setBC
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(VARYING_STRING) :: localError    
    
    ENTERS("LinearElasticity_BoundaryConditionsAnalyticCalculate",err,error,*999)
    
    geometricTolerance=1.0E-6_DP
    
!!TODO: Use Geometric/Material Field values to prescribe values in analytic solution, currently hardcodded geometry & material properties

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)

    NULLIFY(equationsAnalytic)
    CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
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
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    !
    ! IDENTIFY BOUNDARY CONDITION NODES
    !
    bcXNodes = 0
    bcZNodes = 0
    bcXCounter = 0
    bcYCounter = 0
    bcXForceCounter = 0
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
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
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            setBC = .FALSE.
            SELECT CASE(analyticFunctionType)
              !
              ! ONE DIMENSIONAL LINEAR ELASTICITY
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1)
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  !pass
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  !pass
                END SELECT
              END SELECT              
              !
              ! TWO DIMENSIONAL LINEAR ELASTICITY
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1)
              SELECT CASE(componentIdx)
              CASE(1) !u component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              CASE(2) !v component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                      bcXCounter = bcXCounter + 1
                    ENDIF
                    IF(ABS(x(2)-0.0_DP) < geometricTolerance) THEN
                      bcYCounter = bcYCounter + 1
                    ENDIF
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              END SELECT
              !
              ! THREE DIMENSIONAL LINEAR ELASTICITY
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1)
              length=20.0_DP
              width=20.0_DP
              height=5.0_DP
              SELECT CASE(componentIdx)
              CASE(1) !u component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              CASE(2) !v component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    IF(ABS(x(2)-width) < geometricTolerance) THEN
                      IF(ABS(x(1)-length) < geometricTolerance) THEN
                        bcZNodes = bcZNodes + 1
                      ENDIF
                      IF(ABS(x(3)-height) < geometricTolerance) THEN
                        bcXNodes = bcXNodes + 1 
                      ENDIF
                    ENDIF
                  END SELECT
                END SELECT
              CASE(3) !w component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              END SELECT
              !
              ! THREE DIMENSIONAL LINEAR ELASTICITY FLEXURE
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_2)
              !length = 1.0_DP
              !width=2.0_DP
              !height=5.0_DP
              SELECT CASE(componentIdx)
              CASE(1) !u component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              CASE(2) !v component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              CASE(3) !w component
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !pass
                  END SELECT
                END SELECT
              END SELECT              
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
    ENDDO !variableIdx
    !
    ! SET BOUNDARY CONDITIONS & ANALYTIC SOLUTION VALUES
    !
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
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
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            setBC = .FALSE.
            SELECT CASE(analyticFunctionType)
              !
              ! ONE DIMENSIONAL LINEAR ELASTICITY
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1)
              forceX=50.0_DP
              length=20.0_DP
              width=20.0_DP
              height=5.0_DP
              forceXArea=width*height
              E=10.0E3_DP
              SELECT CASE(variableType)
!!TODO set material parameters from material field
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=(x(1)*(forceX/forceXArea))/E
                  IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                    setBC = .TRUE.
                    bcValue=0.0_DP
                  ENDIF
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=1.0_DP
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                    analyticValue=-forceX
                  ELSE
                    analyticValue=0.0_DP
                  ENDIF
                  IF(ABS(x(1)-length) < geometricTolerance) THEN
                    setBC = .TRUE.
                    bcValue=-forceX
                  ENDIF
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=1.0_DP
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
              !
              ! TWO DIMENSIONAL LINEAR ELASTICITY
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1)
              length=20.0_DP
              width=20.0_DP
              height=5.0_DP
              forceY=50.0_DP
              forceYArea=width*height
              E=10.0E3_DP
              vX=0.3_DP
              SELECT CASE(componentIdx)
              CASE(1) !u component
                !u=Sigmax*x/E
                SELECT CASE(variableType)
!!TODO set material parameters from material field
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=(-vX*x(1)*(forceY/forceYArea))/E
                    IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=1.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=1.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=1.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=1.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=1.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=1.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(2) !v component
                !v=Sigmay*y/E
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=(x(2)*(forceY/forceYArea))/E
                    IF(ABS(x(2)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                    !If node located on a line edge of mesh
                    IF(ABS(x(2)-width) < geometricTolerance) THEN !Apply Force BC
                      setBC = .TRUE.
                      IF((ABS(x(1)-0.0_DP) < geometricTolerance).OR.(ABS(x(1)-length) < geometricTolerance)) THEN
                        !lies on a corner of the mesh
                        bcValue=-forceY/(bcYCounter-1.0_DP)*0.5_DP
                      ELSE
                         !does not lie on a corner of the mesh
                        bcValue=-forceY/(bcYCounter-1.0_DP)
                      ENDIF
                    ELSE IF(ABS(x(2)-0) < geometricTolerance) THEN
                      !Provide Analytic reaction force, node located on fixed displacment edge
                      IF((ABS(x(1)-0.0_DP) < geometricTolerance).OR.(ABS(x(1)-length) < geometricTolerance)) THEN
                         !lies on a corner of the mesh
                        analyticValue=-forceY/(bcYCounter-1.0_DP)*0.5_DP
                      ELSE
                         !does not lie on a corner of the mesh
                        analyticValue=-forceY/(bcYCounter-1.0_DP)
                      ENDIF
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
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
                localError="The component index of "//TRIM(NumberToVString(componentIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT              
              !
              ! THREE DIMENSIONAL LINEAR ELASTICITY
              !              
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1)
              length=20.0_DP
              width=20.0_DP
              height=5.0_DP
              forceY=50.0_DP
              forceYArea=width*height
              E=10.0E3_DP
              vX=0.3_DP
              SELECT CASE(componentIdx)
              CASE(1) !u component
                !u=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=(-vX*x(1)*(forceY/forceYArea))/E
                    IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(2) !v component
                !v=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=(x(2)*(forceY/forceYArea))/E
                    IF(ABS(x(2)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                    IF(ABS(x(2)-width) < geometricTolerance) THEN !Apply Force BC
                      setBC = .TRUE.
                      IF(((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance))) THEN
                        !lies on a corner of the mesh
                        bcValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.25_DP
                      ELSE IF((ABS(x(1)-0.0_DP) < geometricTolerance).OR.(ABS(x(1)-length) < geometricTolerance) .OR. &
                        & (ABS(x(3)-0.0_DP) < geometricTolerance).OR.(ABS(x(3)-height) < geometricTolerance)) THEN
                        !lies on an edge of the mesh
                        bcValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.5_DP
                      ELSE !lies on xi2=1 face
                        bcValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))
                      ENDIF
                    ELSE IF(ABS(x(2)-0.0_DP) < geometricTolerance) THEN
                      !Provide Analytic reaction force, node located on fixed displacment edge
                      IF (((ABS(x(1)-0.0_DP) < geometricTolerance) .AND. (ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-length) < geometricTolerance) .AND. (ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-0.0_DP) < geometricTolerance) .AND. (ABS(x(3)-height) < geometricTolerance)) .OR. &
                        & ((ABS(x(1)-length) < geometricTolerance) .AND. (ABS(x(3)-height) < geometricTolerance))) THEN
                         !lies on a corner of the mesh
                        analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.25_DP
                      ELSEIF ((ABS(x(1)-0.0_DP) < geometricTolerance) .OR. (ABS(x(1)-length) < geometricTolerance) .OR. &
                        & (ABS(x(3)-0.0_DP) < geometricTolerance) .OR. (ABS(x(3)-height) < geometricTolerance)) THEN
                         !lies on an edge of the mesh
                        analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.5_DP
                      ELSE
                         !lies on xi2=1 face
                        analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))
                      ENDIF
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*", err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(3) !w component
                !w=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=(-vX*x(3)*(forceY/forceYArea))/E
                    IF(ABS(x(3)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
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
                localError="The component index of "//TRIM(NumberToVString(componentIdx,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT              
              !
              ! THREE DIMENSIONAL LINEAR ELASTICITY FLEXURE
              !
            CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_2)
              length=5.0_DP
              width=2.0_DP
              height=2.0_DP
              forceZ=-1.0_DP
              !forceZ_AREA=width*height
              E=10.0E3_DP
              vX=0.3_DP
              Iyy=(height*(width**3.0_DP))/12.0_DP
              SELECT CASE(componentIdx)
              CASE(1) !u component
                !u=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !analyticValue=(-vX*x(1)*(forceY/forceYArea))/E
                    IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(2) !v component
                !v=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    !analyticValue=(x(2)*(forceY/forceYArea))/E
                    !IF (ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                    !  setBC = .TRUE.
                    !  bcValue=0.0_DP
                    !ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                    !IF(ABS(x(2)-width) < geometricTolerance) THEN !Apply Force BC
                    !  setBC = .TRUE.
                    !  IF(((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance))) THEN
                    !    !lies on a corner of the mesh
                    !    bcValue=forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.25_DP
                    !  ELSE IF((ABS(x(1)-0.0_DP) < geometricTolerance).OR.(ABS(x(1)-length) < geometricTolerance) .OR. &
                    !    & (ABS(x(3)-0.0_DP) < geometricTolerance).OR.(ABS(x(3)-height) < geometricTolerance)) THEN
                    !    !lies on an edge of the mesh
                    !    bcValue=forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.5_DP
                    !  ELSE
                    !    !lies on xi2=1 face
                    !    bcValue=forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))
                    !  ENDIF
                    !ELSE IF(ABS(x(2)-0.0_DP) < geometricTolerance) THEN
                    !  !Provide Analytic reaction force, node located on fixed displacment edge
                    !  IF(((ABS(x(1)-0.0_DP) < geometricTolerance) .AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-0.0_DP) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance)) .OR. &
                    !    & ((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(3)-height) < geometricTolerance))) THEN
                    !    !lies on a corner of the mesh
                    !    analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.25_DP
                    !  ELSE IF((ABS(x(1)-0.0_DP) < geometricTolerance).OR.(ABS(x(1)-length) < geometricTolerance) .OR. &
                    !    & (ABS(x(3)-0.0_DP) < geometricTolerance).OR.(ABS(x(3)-height) < geometricTolerance)) THEN
                    !    !lies on an edge of the mesh
                    !    analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))*0.5_DP
                    !  ELSE
                    !    !lies on xi2=1 face
                    !    analyticValue=-forceY/((bcXNodes-1.0_DP)*(bcZNodes-1.0_DP))
                    !  ENDIF
                    !ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE DEFAULT
                  localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(3) !w component
                !w=
                SELECT CASE(variableType)
                CASE(FIELD_U_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    analyticValue=0.0_DP
                    IF((ABS(x(1)-0.0_DP) < geometricTolerance).AND.(ABS(x(2)-0.0_DP) < geometricTolerance)) THEN
                      analyticValue=(forceZ*(3.0_DP*length-x(3))*x(3)**2.0_DP)/(6.0_DP*Iyy*E)
                    ENDIF
                    IF(ABS(x(1)-0.0_DP) < geometricTolerance) THEN
                      setBC = .TRUE.
                      bcValue=0.0_DP
                    ENDIF
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
                  CASE DEFAULT
                    localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                  SELECT CASE(globalDerivativeIndex)
                  CASE(NO_GLOBAL_DERIV)
                    IF((ABS(x(1)-length) < geometricTolerance).AND.(ABS(x(2)-0.0_DP) < geometricTolerance) .AND. &
                      (ABS(x(3)-0.0_DP) < geometricTolerance)) THEN
                      setBC = .TRUE.
                      bcValue=-forceZ
                    ENDIF
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S2_S3)
                    analyticValue=0.0_DP
                  CASE(GLOBAL_DERIV_S1_S2_S3)
                    analyticValue=0.0_DP
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
                localError="The component index of "//TRIM(NumberToVString(componentIdx,"*",err,error))//" is invalid."
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
              & analyticValue,err,error,*999)
            IF(setBC) THEN
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
              CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
                & bcValue,err,error,*999)
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    
    EXITS("LinearElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("LinearElasticity_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("LinearElasticity_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_BoundaryConditionsAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a linear elasticity finite element equations set.
  SUBROUTINE LinearElasticity_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element 
    !<calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: colsComponentIdx,colsVariableType,columnXiIdx,columnElementDOFIdx,columnElementParameterIdx, &
      & esSpecification(3),fieldVarType,gaussPointIdx,numberOfColumnComponents,numberDependentElementParameters(3), &
      & numberOfColsComponents,numberOfDimensions,numberOfGauss,numberOfRowsComponents,numberOfXi,rowElementParameterIdx, &
      & rowElementDOFIdx,rowXiIdx,rowsVariableType,scalingType,totalNumberDependentElementParameters,xiIdx
    INTEGER(INTG) :: offDiagonalComponents(3),offDiagonalDependentVariable(2,2,3),diagonalSubMatrixLocation(3), &
      & offDiagonalSubMatLocation(2,3)
    REAL(DP) :: colsdPhidXi,gaussWeight,jacobian,jacobianGaussWeight,C(6,6),jacobianGaussWeightDiagC(3,3), &
      & jacobianGaussWeightOffDiagC(2,3),sf(64*3),sourceParam
    LOGICAL :: update,updateMatrix,updateRHS,updateSource
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(BasisPtrType) :: dependentBases(3)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,sourceInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,materialsInterpPoint,sourceInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: quadratureSchemes(3)
    TYPE(VARYING_STRING) :: localError
    TYPE DPHI_DX_COMP_TYPE !A type to store dPhidX for each mesh component
    REAL(DP) :: dPhidX(64,3)
    END TYPE DPHI_DX_COMP_TYPE
    TYPE(DPHI_DX_COMP_TYPE) :: dPhidXComponent(3)

    ENTERS("LinearElasticity_FiniteElementCalculate",err,error,*999)
    
!!Have a look at XPES40.f in the old CMISS code.
!!Q - CPB: Need to think about anisotropic materials with fibre fields.
!!Q - CPB: why store this dPhidX(columnElementParameterIdx,xiIdx) as opposed to just using it directly? A - to minimize operations - Otherwise it would be calculated many more times than
!!         necessary within the loops below 
!!Q - TPBG: Need to be able to use different Quadrature schemes with different bases? A - No use highest quadrature scheme for all directions
!!TODO:: Check whether quadrature scheme being used is suffient to interpolate highest order basis function    
!!Q - TPBG: Need to be able to use different Interpolation for Geometric & Dependent field?

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of a elasticty equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
    ENDIF
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,equationsMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(equationsMatrix,updateMatrix,err,error,*999)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    update=(updateMatrix.OR.updateSource.OR.updateRHS)

    IF(update) THEN
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
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
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      IF(ASSOCIATED(fibreField)) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpParameters,err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpPoint,err,error,*999)
      ENDIF
      
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
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)

      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
            
      NULLIFY(sourceField)
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
      IF(ASSOCIATED(sourceField)) THEN
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpParameters,err,error,*999)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpPoint,err,error,*999)
      ENDIF
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)      
     
!!TODO:: Use highest interpolation scheme's guass points. Warn if Gauss Points insufficient
      !Create an array of Bases with each component 
      DO colsComponentIdx=1,numberOfColsComponents
        NULLIFY(columnDomain)
        CALL FieldVariable_ComponentDomainGet(colsVariable,colsComponentIdx,columnDomain,err,error,*999)
        NULLIFY(columnDomainTopology)
        CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
        NULLIFY(columnDomainElements)
        CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
        NULLIFY(dependentBases(colsComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,dependentBases(colsComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(dependentBases(colsComponentIdx)%ptr, &
          & numberDependentElementParameters(colsComponentIdx),err,error,*999)
        NULLIFY(quadratureSchemes(colsComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(dependentBases(colsComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & quadratureSchemes(colsComponentIdx)%ptr,err,error,*999)
      ENDDO !colsComponentIdx
      totalNumberDependentElementParameters = SUM(numberDependentElementParameters(1:numberOfColsComponents))

      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
        & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
        & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
        & EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
        !
        !ONE, TWO & THREE DIMENSIONAL LINEAR ELASTICITY
        !
        !Loop over gauss points & integrate upper triangular portion of Stiffness matrix
        DO gaussPointIdx=1,numberOfGauss !Gauss point index

          !Parameters for number of off diagonal stress/strain terms for a given number of xi directions and order of
          !calculation for shear terms
          !These parameters do not change for 1D,2D,3D Linear Elasticity
          offDiagonalComponents = [0,1,3]
          offDiagonalDependentVariable(1,1,:) = [1,1,2]
          offDiagonalDependentVariable(1,2,:) = [2,3,3]
          offDiagonalDependentVariable(2,1,:) = offDiagonalDependentVariable(1,2,:)
          offDiagonalDependentVariable(2,2,:) = offDiagonalDependentVariable(1,1,:)
          !
          diagonalSubMatrixLocation(:) = [0,numberDependentElementParameters(1),SUM(numberDependentElementParameters(1:2))]
          offDiagonalSubMatLocation(1,:) = [0,0,numberDependentElementParameters(1)]
          offDiagonalSubMatLocation(2,:) = [numberDependentElementParameters(1),diagonalSubMatrixLocation(3), &
            & diagonalSubMatrixLocation(3)]

          !Interpolate geometric, fibre and material fields at gauss points & calculate geometric field metrics
          CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
            & err,error,*999)
!!TODO:: Add option to only evaluate required metrics
          CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
          IF(ASSOCIATED(fibreField)) &
            & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
            & err,error,*999)
          IF(ASSOCIATED(sourceField)) THEN
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
              & err,error,*999)
            sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
          ENDIF

          !Calculate jacobianGaussWeight.
          CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
          CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
          jacobianGaussWeight=jacobian*gaussWeight
          
          DO colsComponentIdx=1,numberOfColsComponents
            dPhidXComponent(colsComponentIdx)%dPhidX = 0.0_DP
            DO columnElementParameterIdx=1,numberDependentElementParameters(colsComponentIdx)
              DO columnXiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureSchemes(colsComponentIdx)%ptr, &
                  & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx, &
                  & colsdPhidXi,err,error,*999)
                DO rowXiIdx=1,numberOfXi
                  !!TODO: second index in dXidX should be component not xi?
                  dPhidXComponent(colsComponentIdx)%dPhidX(columnElementParameterIdx,rowXiIdx) = &
                    & dPhidXComponent(xiIdx)%dPhidX(columnElementParameterIdx,rowXiIdx)+ &
                    & geometricInterpPointMetrics%dXidX(rowXiIdx,columnXiIdx)*colsdPhidXi
                ENDDO !rowXiIdx
              ENDDO !columnXiIdx
            ENDDO !columnElementParameterIdx
          ENDDO !xiIdx
          !TODO:: what about fibres?
          !Create Linear Elasticity Tensor C
          CALL LinearElasticity_ElasticityTensor(esSpecification(3),materialsInterpPoint,C,err,error,*999)
          !Store Elasticity Tensor diagonal & off diagonal stress coefficients
          jacobianGaussWeightDiagC(3,:) = [jacobianGaussWeight*C(4,4),jacobianGaussWeight*C(5,5),jacobianGaussWeight*C(3,3)]
          jacobianGaussWeightDiagC(2,:) = [jacobianGaussWeight*C(6,6),jacobianGaussWeight*C(2,2),jacobianGaussWeightDiagC(3,2)]
          jacobianGaussWeightDiagC(1,:) = [jacobianGaussWeight*C(1,1),jacobianGaussWeightDiagC(2,1),jacobianGaussWeightDiagC(3,1)]
          jacobianGaussWeightOffDiagC(1,:) = [jacobianGaussWeight*C(1,2),jacobianGaussWeight*C(1,3),jacobianGaussWeight*C(2,3)]
          jacobianGaussWeightOffDiagC(2,:) = [jacobianGaussWeight*C(6,6),jacobianGaussWeight*C(4,4),jacobianGaussWeight*C(5,5)]
          !Construct Element Matrix Diagonal Terms
          DO xiIdx=1,numberOfXi
            DO columnElementParameterIdx=1,numberDependentElementParameters(xiIdx)
              DO rowElementParameterIdx=columnElementParameterIdx,numberDependentElementParameters(xiIdx)
                equationsMatrix%elementMatrix%matrix(diagonalSubMatrixLocation(xiIdx)+ &
                  & columnElementParameterIdx,diagonalSubMatrixLocation(xiIdx)+rowElementParameterIdx) = &
                  & equationsMatrix%elementMatrix%matrix(diagonalSubMatrixLocation(xiIdx)+ &
                  & columnElementParameterIdx,diagonalSubMatrixLocation(xiIdx)+rowElementParameterIdx) + &
                  & DOT_PRODUCT(dPhidXComponent(xiIdx)%dPhidX(columnElementParameterIdx,1:numberOfXi)* &
                  & dPhidXComponent(xiIdx)%dPhidX(rowElementParameterIdx,1:numberOfXi), &
                  & jacobianGaussWeightDiagC(xiIdx,1:numberOfXi))
              ENDDO !rowElementParameterIdx
            ENDDO !columnElementParameterIdx
          ENDDO !xiIdx
          !Construct Element Matrix Off-Diagonal Terms
          DO xiIdx=1,offDiagonalComponents(numberOfXi)
            DO columnElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,1,xiIdx))
              DO rowElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,2,xiIdx))
                equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
                  & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx) = &
                  & equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
                  & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx)+ &
                  & DOT_PRODUCT(dPhidXComponent(offDiagonalDependentVariable(1,1,xiIdx))% &
                  & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(1,:,xiIdx))* &
                  & dPhidXComponent(offDiagonalDependentVariable(1,2,xiIdx))% &
                  & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(2,:,xiIdx)), &
                  & jacobianGaussWeightOffDiagC(:,xiIdx))
              ENDDO !rowElementParameterIdx
            ENDDO !columnElementParameterIdx
          ENDDO !xiIdx
          
          !Below is the full form of constructing the off-Diagonal terms. This will be documented in the linear elasticity
          !equation set page on doxygen for clarity
                        
          !Expanding the DOT_PRODUCT terms
          
          ! offDiagonalDependentVariable(1,:) = [1,1,2]
          ! offDiagonalDependentVariable(2,:) = [2,3,3]
          ! DO xiIdx=1,offDiagonalComponents(numberOfXi)
          !   DO columnElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,1,xiIdx))
          !     DO rowElementParameterIdx=1,numberDependentElementParameters(offDiagonalDependentVariable(1,2,xiIdx))
          !       equatiocolumnElementParameterIdxMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
          !         & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx) = &
          !         & equationsMatrix%elementMatrix%matrix(offDiagonalSubMatLocation(1,xiIdx)+ &
          !         & columnElementParameterIdx,offDiagonalSubMatLocation(2,xiIdx)+rowElementParameterIdx)+ &
          !         & jacobianGaussWeightOffDiagC(1,xiIdx)*dPhidXComponent(offDiagonalDependentVariable(1,xiIdx))% &
          !         & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(1,xiIdx))* &
          !         & dPhidXComponent(offDiagonalDependentVariable(2,xiIdx))% &
          !         & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(2,xiIdx))+ &
          !         & jacobianGaussWeightOffDiagC(2,xiIdx)*dPhidXComponent(offDiagonalDependentVariable(1,xiIdx))% &
          !         & dPhidX(columnElementParameterIdx,offDiagonalDependentVariable(2,xiIdx))* &
          !         & dPhidXComponent(offDiagonalDependentVariable(2,xiIdx))% &
          !         & dPhidX(rowElementParameterIdx,offDiagonalDependentVariable(1,xiIdx))
          !     ENDDO !rowElementParameterIdx
          !   ENDDO !columnElementParameterIdx
          ! ENDDO !xiIdx
          
          ! !Expanding the xiIdx loop above
          
          ! DO columnElementParameterIdx=1,numberDependentElementParameters(1)
          !   DO rowElementParameterIdx=1,numberDependentElementParameters(2)
          !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
          !       & numberDependentElementParameters(1))= &
          !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
          !       & numberDependentElementParameters(1))+ &
          !       & jacobianGaussWeight*C(1,2)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,1)* &
          !       & dPhidXComponent(2)%dPhidX(rowElementParameterIdx,2)+ &
          !       & jacobianGaussWeight*C(6,6)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,2)* &
          !       & dPhidXComponent(2)%dPhidX(rowElementParameterIdx,1)
          !   ENDDO !columnElementParameterIdx
          ! ENDDO !rowElementParameterIdx
          ! DO columnElementParameterIdx=1,numberDependentElementParameters(1)
          !   DO rowElementParameterIdx=1,numberDependentElementParameters(3)
          !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
          !       & numberDependentElementParameters(1)+numberDependentElementParameters(2))  &
          !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx,rowElementParameterIdx+ &
          !       & numberDependentElementParameters(1)+numberDependentElementParameters(2)) + &
          !       & jacobianGaussWeight*C(1,3)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,1)* &
          !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,3)+ &
          !       & jacobianGaussWeight*C(4,4)*dPhidXComponent(1)%dPhidX(columnElementParameterIdx,3)* &
          !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,1)
          !   ENDDO !columnElementParameterIdx
          ! ENDDO !rowElementParameterIdx
          ! DO columnElementParameterIdx=1,numberDependentElementParameters(2)
          !   DO rowElementParameterIdx=1,numberDependentElementParameters(3)
          !     equationsMatrix%elementMatrix%matrix(columnElementParameterIdx+ &
          !       & numberDependentElementParameters(1),rowElementParameterIdx+numberDependentElementParameters(1)+ &
          !       & numberDependentElementParameters(2))= &
          !       & equationsMatrix%elementMatrix%matrix(columnElementParameterIdx+ &
          !       & numberDependentElementParameters(1),rowElementParameterIdx+numberDependentElementParameters(1)+ &
          !       & numberDependentElementParameters(2))+ &
          !       & jacobianGaussWeight*C(2,3)*dPhidXComponent(2)%dPhidX(columnElementParameterIdx,2)* &
          !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,3)+ &
          !       & jacobianGaussWeight*C(5,5)*dPhidXComponent(2)%dPhidX(columnElementParameterIdx,3)* &
          !       & dPhidXComponent(3)%dPhidX(rowElementParameterIdx,2)
          !   ENDDO !columnElementParameterIdx
          ! ENDDO !rowElementParameterIdx
          
        ENDDO !gaussPointIdx
        
        !If Plane Stress/Strain problem multiply equation matrix by thickness
        IF(esSpecification(3) == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE .OR. &
          & esSpecification(3) == EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE .OR. & 
          & esSpecification(3) == EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE) THEN
          DO rowElementDOFIdx=1,totalNumberDependentElementParameters
            DO columnElementDOFIdx=rowElementDOFIdx,totalNumberDependentElementParameters
!!TODO::Bring 2D plane stress/strain element thickness in through a field - element constant when it can be exported by field i/o.
!!      Currently brought in through material field (Temporary)
              
              equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                & materialsInterpPoint%values(1,NO_PART_DERIV)
            ENDDO !columnElementDOFIdx
          ENDDO !rowElementDOFIdx
        ENDIF
      
!!TODO:: Is this RHS Vector update required? find out/check - RHS not used - BC are prescribed during assembling
!!       eg update RHS only when BC change - stiffness matrix should be the same
        IF(updateRHS) THEN
          rhsVector%elementVector%vector=0.0_DP
        ENDIF
        
        !Scale factor adjustment, Application of scale factors is symmetric
        CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
        IF(scalingType/=FIELD_NO_SCALING) THEN
          NULLIFY(colsInterpParameters)
          CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,colsInterpParameters, &
            & err,error,*999)
          DO xiIdx=1,numberOfXi
            sf(diagonalSubMatrixLocation(xiIdx)+1:SUM(numberDependentElementParameters(1:xiIdx)))= &
              & colsInterpParameters%scaleFactors(:,xiIdx)
          ENDDO !xiIdx
          DO rowElementDOFIdx=1,totalNumberDependentElementParameters
            IF(updateMatrix) THEN
              DO columnElementDOFIdx=rowElementDOFIdx,totalNumberDependentElementParameters
                equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                & SF(rowElementDOFIdx)*SF(columnElementDOFIdx)
              ENDDO !columnElementDOFIdx
            ENDIF
!!TODO:: Check if RHS update required for Linear Elasticity ie is the RHS the force terms but they are set during assembling and not here?
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)= &
                & rhsVector%elementVector%vector(rowElementDOFIdx)*SF(rowElementDOFIdx)
            ENDIF
          ENDDO !rowElementDOFIdx
          
          IF(updateMatrix) THEN
            !Transpose upper triangular portion of Stiffness matrix to give lower triangular portion.
            !Has to be done after scale factors are applied
!!TODO:: Use symmetric linear equation solver or alternatively traspose to give full matrix when asemmbling or when
!!       creating solver matrices        
!!TODO:: Better to use SIZE(equationsMatrix%elementMatrix%matrix,1) as apposed to totalNumberDependentElementParameters?
!!       Is the size re-calculated at end of every loop?
            DO rowElementDOFIdx=2,totalNumberDependentElementParameters
              DO columnElementDOFIdx=1,rowElementDOFIdx-1
                equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  &  equationsMatrix%elementMatrix%matrix(columnElementDOFIdx,rowElementDOFIdx)
              ENDDO !columnElementDOFIdx
            ENDDO !rowElementDOFIdx
          ENDIF
        ENDIF !scaling
        
      CASE(EQUATIONS_SET_PLATE_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_SHELL_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Linear Elasticity equation type of a Elasticty equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    ENDIF !update

    EXITS("LinearElasticity_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("LinearElasticity_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the linear elasticity tensor
  SUBROUTINE LinearElasticity_ElasticityTensor(equationsSetSubtype,materialsInterpolatedPoint,elasticityTensor,err,error,*)

    !Argument variables    
    INTEGER(INTG), INTENT(IN) :: equationsSetSubtype !<The subtype of the particular equation set being used
    TYPE(FieldInterpolatedPointType), POINTER :: materialsInterpolatedPoint
    REAL(DP), INTENT(OUT) :: elasticityTensor(:,:) !<The Linear Elasticity Tensor C
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    REAL(DP) :: E1,E2,E3,v13,v23,v12,v31,v32,v21,gama
    REAL(DP) :: C11,C22,C33,C12,C13,C23,C21,C31,C32,C44,C55,C66
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_ElasticityTensor",err,error,*999)
    elasticityTensor=0.0_DP
    SELECT CASE(equationsSetSubtype)
    !Note: Fortran uses column major format for arrays.
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
      !General Orthotropic 3D Linear Elasticity Tensor
      E1 = materialsInterpolatedPoint%values(1,1)
      E2 = materialsInterpolatedPoint%values(2,1)
      E3 = materialsInterpolatedPoint%values(3,1)
      v13 = materialsInterpolatedPoint%values(4,1)
      v23 = materialsInterpolatedPoint%values(5,1)
      v12 = materialsInterpolatedPoint%values(6,1)
      v31 = v13
      v32 = v23
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12*v21-v23*v32-v31*v13-2.0_DP*v21*v32*v13)
      C11 = E1*(1.0_DP-v23*v32)*gama
      C22 = E2*(1.0_DP-v13*v31)*gama
      C33 = E3*(1.0_DP-v12*v21)*gama
      C12 = E1*(v21+v31*v23)*gama ! = E2*(v12+v32*v13)*gama
      C13 = E1*(v31+v21*v32)*gama ! = E3*(v13+v12*v23)*gama
      C23 = E2*(v32+v12*v31)*gama ! = E3*(v23+v21*v13)*gama
      C21 = C12
      C31 = C13
      C32 = C23
      C44 = E2/(2.0_DP*(1.0_DP+v23)) != G23
      C55 = E1/(2.0_DP*(1.0_DP+v13)) != G13
      C66 = E3/(2.0_DP*(1.0_DP+v12)) != G12
      elasticityTensor(1:6,1)=[C11,C21,C31,0.0_DP,0.0_DP,0.0_DP]
      elasticityTensor(1:6,2)=[C12,C22,C32,0.0_DP,0.0_DP,0.0_DP]
      elasticityTensor(1:6,3)=[C13,C23,C33,0.0_DP,0.0_DP,0.0_DP]
      elasticityTensor(4,4)=C44
      elasticityTensor(5,5)=C55
      elasticityTensor(6,6)=C66
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE)
      !Plane Stress Isotropic Elasticity Tensor
      E1 = materialsInterpolatedPoint%values(2,1)
      v12 = materialsInterpolatedPoint%values(3,1)
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12*v21)
      C11 = E1*gama
      C22 = C11
      C12 = C11*v21
      C21 = C12
      C66 = E1/(2.0_DP*(1.0_DP+v12)) != G12
      elasticityTensor(1,1)=C11
      elasticityTensor(1,2)=C21
      elasticityTensor(2,1)=C21
      elasticityTensor(2,2)=C22
      elasticityTensor(6,6)=C66
    CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
      !Plane Strain Isotropic Linear Elasticity Tensor
      E1 = materialsInterpolatedPoint%values(2,1)
      E2 = E1
      v12 = materialsInterpolatedPoint%values(3,1)
      v21 = v12
      gama = 1.0_DP/(1.0_DP-v12-v21)
      C11 = E1*gama*(1.0_DP-v12)/(1.0_DP+v12)
      C22 = E2*gama*(1.0_DP-v21)/(1.0_DP+v21)
      C12 = C22*v12
      C21 = C12
      C66 = E1/(2.0_DP*(1.0_DP+v12)) != G12
      elasticityTensor(1,1)=C11
      elasticityTensor(1,2)=C21
      elasticityTensor(2,1)=C21
      elasticityTensor(2,2)=C22
      elasticityTensor(6,6)=C66
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
      !Plane Strain Isotropic Linear Elasticity Tensor
      E1 = materialsInterpolatedPoint%values(2,1)
      C11 = E1
      elasticityTensor(1,1)=C11
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a Linear Elasticity equation type of a Elasticty equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("LinearElasticity_ElasticityTensor")
    RETURN
999 ERRORSEXITS("LinearElasticity_ElasticityTensor",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ElasticityTensor

  !
  !================================================================================================================================
  !

  !>Sets up the Linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a linear elasticity equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricComponentNumber,geometricMeshComponent, &
      & geometricScalingType,numberOfComponents,numberOfDimensions,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of an elasticity equation set class."
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
        !Default to FEM solution
        CALL LinearElasticity_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !Do nothing???
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          numberOfComponents=numberOfDimensions
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
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
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          numberOfComponents=numberOfDimensions
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
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
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation"
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
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
        numberOfComponents=2
      CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE,EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
        numberOfComponents=3
      CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
        numberOfComponents=6
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a linear elasticity equation type of an elasticity equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Default to the general 3D orthotropic material
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,numberOfComponents
            !Default to to the first geometric component with constant interpolation
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
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
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
            !Set the default values for the materials field
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,err,error,*999) !Young's Modulus
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,0.25_DP,err,error,*999) !Poisson's Ratio
          CASE(EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE,EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE)
            !2 Components for 2D Isotropic Linear Elasticity
            !TODO:: Temporarily set to 3 to allow thickness to passed in.
            !Remove once a thickness, element constant field is defined and can be exported/viewed by cmgui
            !Set the default values for the materials field
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,1,30.0E6_DP,err,error,*999) !Young's Modulus
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,2,0.25_DP,err,error,*999) !Poisson's Ratio
          CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
            !Set the default values for the materials field
            DO componentIdx=1,3
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,30.0E6_DP,err,error,*999)
            ENDDO !componentIdx
            DO componentIdx=4,6
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & componentIdx,0.25_DP,err,error,*999)
            ENDDO !componentIdx
          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is not valid for a linear elasticity equation type of an elasticity equation set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
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
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        !List 3 Dimensional Analytic function types currently implemented
        SELECT CASE(equationsSetSetup%analyticFunctionType)
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1)
          !Check that we are in 1D
          IF(numberOfDimensions/=1) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 1 geometric dimension."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_ELASTICITY_ONE_DIM_1
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1)
          !Check that we are in 2D
!!TODO:: This check may have been done before
          IF(numberOfDimensions/=2) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 2 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_ELASTICITY_TWO_DIM_1
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1)
          !Check that we are in 3D
          IF(numberOfDimensions/=3) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 3 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_1
        CASE(EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_2)
          !Check that we are in 3D
          IF(numberOfDimensions/=3) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 3 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_ELASTICITY_THREE_DIM_2
        CASE DEFAULT
          localError="The specified analytic function type of "// &
            & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
            & " is invalid for a standard Linear Elasticity equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(analyticField,err,error,*999)
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Linear Elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        !Create the equations
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
            & err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a linear elasticity equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("LinearElasticity_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("LinearElasticity_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE)
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
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equation type of an elasticity equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("LinearElasticity_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("LinearElasticity_EquationsSetSolutionMethodSet",err,error)
    EXITS("LinearElasticity_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a linear elasticity equation type of an elasticity equations set class.
  SUBROUTINE LinearElasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    SELECT CASE(specification(3))
    CASE(EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRESS_SUBTYPE, &
      & EQUATIONS_SET_TWO_DIMENSIONAL_PLANE_STRAIN_SUBTYPE, &
      & EQUATIONS_SET_ONE_DIMENSIONAL_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_PLATE_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_SHELL_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a linear elasticity equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_ELASTICITY_CLASS,EQUATIONS_SET_LINEAR_ELASTICITY_TYPE,specification(3)]
 
    EXITS("LinearElasticity_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("LinearElasticity_EquationsSetSpecificationSet",err,error)
    EXITS("LinearElasticity_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear elasticity problem.
  SUBROUTINE LinearElasticity_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a linear elasticity equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("LinearElasticity_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_NO_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a linear elasticity type of an elasticity problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
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
          & " is invalid for a linear elasticity problem."
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
          & " is invalid for a linear elasticity problem."
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
        !Set the solver to be a linear solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
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
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
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
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Create the solver equations
        NULLIFY(solverEquations)
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solver equations
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear elasticity problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a linear elasticity problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("LinearElasticity_ProblemSetup")
    RETURN
999 ERRORSEXITS("LinearElasticity_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a linear elasticity type problem.
  SUBROUTINE LinearElasticity_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specifiation to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("LinearElasticity_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_NO_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a linear elasticity problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_ELASTICITY_CLASS, PROBLEM_LINEAR_ELASTICITY_TYPE, problemSubtype]

    EXITS("LinearElasticity_ProblemSpecificationSet")
    RETURN
999 ERRORS("LinearElasticity_ProblemSpecificationSet",err,error)
    EXITS("LinearElasticity_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE LinearElasticity_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

END MODULE LinearElasticityRoutines
