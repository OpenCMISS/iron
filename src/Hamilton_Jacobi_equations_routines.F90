!!> \file
!> \author Chris Bradley
!> \brief This module handles all Hamilton-Jacobi equations routines.
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

!>This module handles all Hamilton-Jacobi equations routines.
MODULE HamiltonJacobiRoutines

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
  USE MatrixVector
  USE Maths
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

!!MERGE: move

  PUBLIC HamiltonJacobi_BoundaryConditionsAnalyticCalculate
  
  PUBLIC HamiltonJacobi_EquationsSetSolutionMethodSet
  
  PUBLIC HamiltonJacobi_EquationsSetSetup
  
  PUBLIC HamiltonJacobi_EquationsSetSpecificationSet

  PUBLIC HamiltonJacobi_FiniteElementCalculate
  
  PUBLIC HamiltonJacobi_ProblemSetup
  
  PUBLIC HamiltonJacobi_ProblemSpecificationSet
  
  PUBLIC NUMBER_OF_INPUT_NODES,PRE_PROCESS_INFORMATION,SOLVE_PROBLEM_FMM,SOLVE_PROBLEM_GEODESIC
  
  PUBLIC SOLVE_PROBLEM_GEODESIC_CONNECTIVITY,SOLVE_PROBLEM_FMM_CONNECTIVITY
  
  PUBLIC FIND_MINIMAX,POST_PROCESS_DATA

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE HamiltonJacobi_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,dependentVariableType,derivativeIdx,dimensionIdx,esSpecification(3), &
      & globalDerivativeIndex,localDOFIdx,nodeIdx,numberOfDependentComponents,numberOfDependentVariables,numberOfDimensions, &
      & numberOfNodes,numberOfNodeDerivatives,variableIdx,variableType
    REAL(DP) :: analyticValue,X(3)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(VARYING_STRING) :: localError    
    
    ENTERS("HamiltonJacobi_BoundaryConditionsAnalyticCalculate",err,error,*999)

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
    CALL Field_NumberOfVariablesGet(dependentField,numberOfDependentVariables,err,error,*999)
    DO variableIdx=1,numberOfDependentVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,dependentVariableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfDependentComponents,err,error,*999)
      DO componentIdx=1,numberOfDependentComponents
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
            CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1)
              !u=x^2+2.x.y-y^2
              SELECT CASE(dependentVariableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=x(1)*x(1)-2.0_DP*x(1)*x(2)-x(2)*x(2)
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=2.0_DP*x(1)+2.0_DP*x(2)
                CASE(GLOBAL_DERIV_S2)
                  analyticValue=2.0_DP*x(1)-2.0_DP*x(2)
                CASE(GLOBAL_DERIV_S1_S2)
                  analyticValue=2.0_DP
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=0.0_DP !!TODO
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
                localError="The dependent variable type of "//TRIM(NumberToVString(dependentVariableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2)
              !u=cos(x).cosh(y)
              SELECT CASE(dependentVariableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=COS(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=-SIN(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S2)
                  analyticValue=COS(x(1))*SINH(x(2))
                CASE(GLOBAL_DERIV_S1_S2)
                  analyticValue=-SIN(x(1))*SINH(x(2))
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=0.0_DP !!TODO
                CASE(GLOBAL_DERIV_S1)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The dependent variable type of "//TRIM(NumberToVString(dependentVariableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1)
              !u=x^2+y^2-2.z^2
              SELECT CASE(dependentVariableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=X(1)*X(1)+X(2)*X(2)-2.0_DP*X(3)*X(3)
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=2.0_DP*X(1)
                CASE(GLOBAL_DERIV_S2)
                  analyticValue=2.0_DP*X(2)
                CASE(GLOBAL_DERIV_S1_S2)
                  analyticValue=0.0_DP
                CASE(GLOBAL_DERIV_S3)
                  analyticValue=-4.0_DP*X(3)
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
                  analyticValue=0.0_DP !!TODO
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S3)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S3)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2_S3)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2_S3)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The dependent variable type of "//TRIM(NumberToVString(dependentVariableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2)
              !u=cos(x).cosh(y).z
              SELECT CASE(dependentVariableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=COS(X(1))*COSH(X(2))*X(3)
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=-SIN(X(1))*COSH(X(2))*X(3)
                CASE(GLOBAL_DERIV_S2)
                  analyticValue=COS(X(1))*SINH(X(2))*X(3)
                CASE(GLOBAL_DERIV_S1_S2)
                  analyticValue=-SIN(X(1))*SINH(X(2))*X(3)
                CASE(GLOBAL_DERIV_S3)
                  analyticValue=COS(X(1))*COSH(X(2))
                CASE(GLOBAL_DERIV_S1_S3)
                  analyticValue=-SIN(X(1))*COSH(X(2))
                CASE(GLOBAL_DERIV_S2_S3)
                  analyticValue=COS(X(1))*SINH(X(2))
                CASE(GLOBAL_DERIV_S1_S2_S3)
                  analyticValue=-SIN(X(1))*SINH(X(2))
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=0.0_DP !!TODO
                CASE(GLOBAL_DERIV_S1)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S3)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S3)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2_S3)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2_S3)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(dependentVariableType,"*",err,error))// &
                  & " is invalid."
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
            IF(dependentVariableType==FIELD_U_VARIABLE_TYPE) THEN
              IF(boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
                CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx,BOUNDARY_CONDITION_FIXED, &
                  & analyticValue,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    
    EXITS("HamiltonJacobi_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("HamiltonJacobi_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("HamiltonJacobi_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Hamilton-Jacobi equation finite element equations set.
  SUBROUTINE HamiltonJacobi_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx, &
      & esSpecification(3),gaussPointIdx,numberOfColComponents,numberOfColElementParameters,numberOfDependentXi, &
      & numberOfDimensions,numberOfGauss,numberOfGeometricXi,numberOfRowComponents,numberOfRowElementParameters, &
      & rowComponentIdx,rowElementDOFIdx,rowXiIdx,rowElementParameterIdx,rowsVariableType,scalingType
    REAL(DP) :: gaussWeight,jacobian,jacobianGaussWeight,sum,dRowPhidXi(3),dColPhidXi(3)
    LOGICAL :: update,updateMatrix,updateRHS
    TYPE(BasisType), POINTER :: colBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: colDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: colDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: colDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,geometricInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,fieldVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: colQuadratureScheme,dependentQuadratureScheme,rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("HamiltonJacobi_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Hamilton-Jacobi equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(linearMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,linearMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF
    
    update=(updateMatrix.OR.updateRHS)

    IF(update) THEN

      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      
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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfGeometricXi,err,error,*999)

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
      CALL Basis_NumberOfXiGet(dependentBasis,numberOfDependentXi,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
     
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColComponents,err,error,*999)
      
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      
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
   
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss

        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfDependentXi,geometricInterpPointMetrics,err,error,*999)
        !Calculate jacobianGaussWeight.

        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        
!!TODO: Think about symmetric problems. 
        jacobianGaussWeight=jacobian*gaussWeight
        
        !Loop over row components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          NULLIFY(rowQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            DO rowXiIdx=1,numberOfDependentXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,dRowPhidXi(rowXiIdx),err,error,*999)
            ENDDO !rowXiIdx
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColComponents
                NULLIFY(colDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colDomain,err,error,*999)
                NULLIFY(colDomainTopology)
                CALL Domain_DomainTopologyGet(colDomain,colDomainTopology,err,error,*999)
                NULLIFY(colDomainElements)
                CALL DomainTopology_DomainElementsGet(colDomainTopology,colDomainElements,err,error,*999)
                NULLIFY(colBasis)
                CALL DomainElements_ElementBasisGet(colDomainElements,elementNumber,colBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colBasis,numberOfColElementParameters,err,error,*999)
                NULLIFY(colQuadratureScheme)
                CALL Basis_QuadratureSchemeGet(colBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colQuadratureScheme,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  DO columnXiIdx=1,numberOfDependentXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,dColPhidXi(columnXiIdx),err,error,*999)
                  ENDDO !columnXiIdx                  
                  sum=0.0_DP
                  DO rowXiIdx=1,numberOfDependentXi
                    DO columnXiIdx=1,numberOfDependentXi
                      sum=sum+dRowPhidXi(rowXiIdx)*dColPhidXi(columnXiIdx)*geometricInterpPointMetrics%gu(rowXiIdx,columnXiIdx)
                    ENDDO !columnXiIdx
                  ENDDO !rowXiIdx
                  linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF
            IF(updateRHS) rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDDO !gaussPointIdx
          
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
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowComponents
          NULLIFY(rowDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
          NULLIFY(rowDomainTopology)
          CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
          NULLIFY(rowDomainElements)
          CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
          NULLIFY(rowBasis)
          CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColComponents
                NULLIFY(colDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colDomain,err,error,*999)
                NULLIFY(colDomainTopology)
                CALL Domain_DomainTopologyGet(colDomain,colDomainTopology,err,error,*999)
                NULLIFY(colDomainElements)
                CALL DomainTopology_DomainElementsGet(colDomainTopology,colDomainElements,err,error,*999)
                NULLIFY(colBasis)
                CALL DomainElements_ElementBasisGet(colDomainElements,elementNumber,colBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colBasis,numberOfColElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                    & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                    & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
            IF(updateRHS) rhsVector%elementVector%vector(rowElementDOFIdx)= &
              & rhsVector%elementVector%vector(rowElementDOFIdx)* &
              & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
      
    ENDIF !update
       
    EXITS("HamiltonJacobi_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("HamiltonJacobi_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the Hamilton-Jacobi equation type of a classical field equations set class.
  SUBROUTINE HamiltonJacobi_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a Hamilton-Jacobi equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricMeshComponent,geometricScalingType,numberOfDimensions, &
      & numberOfMaterialsComponents,solutionMethod,sparsityType
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
    
    ENTERS("HamiltonJacobi_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))      
    CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Hamilton-Jacobi equation type of a classical field equation set class."
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
        CALL HamiltonJacobi_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      !Do nothing
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
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
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
            & err,error,*999)
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
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,"Phi",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & "del Phi/del n",err,error,*999)
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
              & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation"
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
      numberOfMaterialsComponents=1
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          !Set the number of materials components
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          !Default the materials components to the geometric interpolation setup with constant interpolation
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
            & err,error,*999)
          DO componentIdx=1,numberOfMaterialsComponents
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
            & err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          DO componentIdx=1,numberOfMaterialsComponents
            CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
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
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(equationsSetSetup%analyticFunctionType)
        CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1)
          !Check that we are in 2D
          IF(numberOfDimensions/=2) THEN
            localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 2 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_HJ_EQUATION_TWO_DIM_1
        CASE(EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2)
          !Check that we are in 2D
          IF(numberOfDimensions/=2) THEN
            localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 2 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_HJ_EQUATION_TWO_DIM_2
        CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1)
          !Check that we are in 3D
          IF(numberOfDimensions/=3) THEN
            localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 3 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_HJ_EQUATION_THREE_DIM_1
        CASE(EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2)
          !Check that we are in 3D
          IF(numberOfDimensions/=3) THEN
            localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 3 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analtyic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_HJ_EQUATION_THREE_DIM_2
        CASE DEFAULT
          localError="The specified analytic function type of "// &
            & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
            & " is invalid for a standard Hamilton-Jacobi equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsSet%analytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        !Create the equations
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
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
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Hamilton-Jacobi equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("HamiltonJacobi_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("HamiltonJacobi_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Hamilton-Jacobi equation type of an classical field equations set class.
  SUBROUTINE HamiltonJacobi_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("HamiltonJacobi_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
   
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
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
        & " is not valid for a Hamilton-Jacobi equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("HamiltonJacobi_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("HamiltonJacobi_EquationsSetSolutionMethodSet",err,error)
    EXITS("HamiltonJacobi_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Hamilton-Jacobi equation type of a classical field equations set class.
  SUBROUTINE HamiltonJacobi_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("HamiltonJacobi_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_STANDARD_HJ_SUBTYPE)
      !ok
    CASE(EQUATIONS_SET_GENERALISED_HJ_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a Hamilton-Jacobi type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_HJ_EQUATION_TYPE,subtype]
 
    EXITS("HamiltonJacobi_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("HamiltonJacobi_EquationsSetSpecificationSet",err,error)
    EXITS("HamiltonJacobi_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !
 
  !>Sets up the Hamilton-Jacobi problem.
  SUBROUTINE HamiltonJacobi_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a Hamilton-Jacobi equation on.
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
    
    ENTERS("HamiltonJacobi_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_HJ_SUBTYPE)
      !OK
    CASE(PROBLEM_GENERALISED_HJ_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Hamilton-Jacobi equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
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
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Finish the control loops
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)            
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !Get the control loop
      NULLIFY(controlLoopRoot)
      CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
      CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Start the solvers creation
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
        !Set the solver to be a linear solver
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
        !Set solver defaults
        CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        !Create the solver equations
        CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        !Get the solver equations
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        !Finish the solver equations creation
        CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)             
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Hamilton-Jacobi equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Hamilton-Jacobi equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("HamiltonJacobi_ProblemSetup")
    RETURN
999 ERRORSEXITS("HamiltonJacobi_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Hamilton-Jacobi equation type.
  SUBROUTINE HamiltonJacobi_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("HamiltonJacobi_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_STANDARD_HJ_SUBTYPE)
      !ok
    CASE(PROBLEM_GENERALISED_HJ_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Hamilton-Jacobi type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_HJ_EQUATION_TYPE,problemSubtype]

    EXITS("HamiltonJacobi_ProblemSpecificationSet")
    RETURN
999 ERRORS("HamiltonJacobi_ProblemSpecificationSet",err,error)
    EXITS("HamiltonJacobi_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE HamiltonJacobi_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates to give back the number of nodes from input file.
  SUBROUTINE NUMBER_OF_INPUT_NODES(INPUT_FILE_NAME,INPUT_FILE_FORMAT,totalNumberOfNodes,TOTAL_NUMBER_OF_ELEMENTS,&
  &TOTAL_NUMBER_OF_CONNECTIVITY,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(OUT) :: totalNumberOfNodes
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_CONNECTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10) :: INPUT_FILE_FORMAT
    INTEGER(INTG) :: Err

    !Argument variables
!    TYPE(VARYING_STRING) :: localError !<The error string

    !Local variables
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: CONNECTIVITY_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:,:):: ELEMENT_LIST
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: CONNECTIVITY_NUMBER
    INTEGER(INTG) :: I,J,K,N,NUMBER_OF_NODES_PER_ELEMENT,THERE_IS_IN_CONNECTIVITY_LIST
    CHARACTER (LEN=300) :: STRING
        
!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)

    
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
!      PRINT *, STRING
      totalNumberOfNodes=-1
      DO WHILE (STRING .ne. "Connectivity") 
        READ(11,*) STRING
        totalNumberOfNodes=totalNumberOfNodes+1
      ENDDO
      
      TOTAL_NUMBER_OF_CONNECTIVITY=0
      DO I=1,totalNumberOfNodes
        READ(11,*) STRING,J
        TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+J
      ENDDO
      
      CLOSE (11)
      TOTAL_NUMBER_OF_ELEMENTS = TOTAL_NUMBER_OF_CONNECTIVITY ! SHOULD BE DEFINED
    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,totalNumberOfNodes
!      PRINT *, STRING,totalNumberOfNodes
      
      DO I=1,INT((totalNumberOfNodes/3.0_DP)+0.5_DP)
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
!      print*,I,INT(totalNumberOfNodes/3.0+0.5),STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS
    
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(totalNumberOfNodes,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(totalNumberOfNodes),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".vtk"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING,totalNumberOfNodes,STRING
!      PRINT *, STRING,totalNumberOfNodes,STRING
      
      DO I=1,totalNumberOfNodes
      
        READ(11,*) STRING

      ENDDO
!      READ(11,*) STRING
      READ(11,*) STRING,TOTAL_NUMBER_OF_ELEMENTS

      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(totalNumberOfNodes,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(totalNumberOfNodes),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".pts"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) totalNumberOfNodes
!      PRINT *, totalNumberOfNodes
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".elem"
!      PRINT *, STRING
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(totalNumberOfNodes,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(totalNumberOfNodes),STAT=ERR)

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

      CLOSE (11)

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".node"
      OPEN (11,FILE=STRING)
      READ(11,*) totalNumberOfNodes
      CLOSE (11)
      
      STRING = INPUT_FILE_NAME(1:I)//".ele"
      OPEN (11,FILE=STRING)
      READ(11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      ALLOCATE(ELEMENT_LIST(TOTAL_NUMBER_OF_ELEMENTS,20),STAT=ERR)
      ALLOCATE(CONNECTIVITY_LIST(totalNumberOfNodes,50),STAT=ERR)
      ALLOCATE(CONNECTIVITY_NUMBER(totalNumberOfNodes),STAT=ERR)
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        DO J=1,20
          ELEMENT_LIST(I,J)=0
        ENDDO
      ENDDO
      DO I=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(I)=0
        DO J=1,50
          CONNECTIVITY_LIST(I,J)=0
        ENDDO
      ENDDO

      TOTAL_NUMBER_OF_CONNECTIVITY=0
      
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

              TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY+1
              
            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE (11)
      
      TOTAL_NUMBER_OF_CONNECTIVITY=TOTAL_NUMBER_OF_CONNECTIVITY*2

    ENDIF


  END SUBROUTINE NUMBER_OF_INPUT_NODES

  !
  !================================================================================================================================
  !


  !>to READ input file.
  SUBROUTINE PRE_PROCESS_INFORMATION(MATERIAL_BEHAVIOUR,INPUT_FILE_NAME,INPUT_FILE_FORMAT,totalNumberOfNodes,&
&INPUT_TYPE_FOR_SEED_VALUE,INPUT_TYPE_FOR_SPEED_FUNCTION,SPEED_FUNCTION_ALONG_EIGEN_VECTOR,INPUT_TYPE_FOR_CONDUCTIVITY,&
&STATUS_MASK,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,&
&SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)

    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: COLUMN_INDEX
    INTEGER(INTG), ALLOCATABLE, DIMENSION(:)  :: RAW_INDEX
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    REAL(DP) :: SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
    INTEGER(INTG), INTENT(OUT) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    INTEGER(INTG), INTENT(OUT) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SEED_VALUE
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_SPEED_FUNCTION
    CHARACTER (LEN=10)  :: INPUT_TYPE_FOR_CONDUCTIVITY
    CHARACTER (LEN=300) :: INPUT_FILE_NAME
    CHARACTER (LEN=10)  :: INPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: I,J,K,N,FIRST_NODE_NUMBER
    INTEGER(INTG) :: TEXT_LENGTH, THERE_IS_IN_CONNECTIVITY_LIST
    REAL(DP) :: A(3),B(3),C(3)
    REAL(DP) :: DOT_PRODUCT_VALUE
        
!INITIALIZE PARAMETERS:
    DO I=1,totalNumberOfNodes
      CONNECTIVITY_NUMBER(I)=0
      DO J=1,3
        NODE_LIST(I,J) = 0.0
        SPEED_FUNCTION_TABLE(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR(I,J) = 0.0
      ENDDO
      RAW_INDEX(I)=0
    ENDDO
    RAW_INDEX(totalNumberOfNodes+1)=0

    DO I=1,TOTAL_NUMBER_OF_CONNECTIVITY
      COLUMN_INDEX(I)=0
      DO J=1,3
        SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
      DO J=1,9
        CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(I,J) = 0.0
      ENDDO
    ENDDO
     
! """""""""""""""""""""""""""""""""""INPUT OF TABC FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TABC") THEN 

      I = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:I)//".tabc"
      OPEN (11,FILE=STRING)
      READ(11,*) STRING

! SOSIOISOISOISOISOIS      load data for the case material behaves = ISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!      input type for *velocity function* = FILE and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,totalNumberOfNodes 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),SEED_VALUE(I)
            
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
           
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
          DO I=1,totalNumberOfNodes 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SEED_VALUE(I)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

          ENDDO
        ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,totalNumberOfNodes 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
          DO I=1,totalNumberOfNodes 

            READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP

            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)

            IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
              SEED_VALUE(I) = 1000.0_DP
            ENDIF

          ENDDO
        ENDIF

      ENDIF

! ANANANAIANSOANAIANI      load data for the case material behaves = ANISOTROPIC 
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 


!      if conductivity format is TENSOR type i.e. three EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,9)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

        ENDIF


!      if conductivity format is VECTOR type i.e. first EigenVectors
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN

!      input type for *velocity function* = FILE and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3),SEED_VALUE(I)

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = FILE
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SEED_VALUE(I)
 
            ENDDO
          ENDIF

!      input type for *velocity function* = FILE and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1),&
                           &SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2),SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

!      input type for *velocity function* = FIXED and *seed points* = LIST
          IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE" .AND. INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN
            DO I=1,totalNumberOfNodes 

              READ(11,*) STRING,(NODE_LIST(I,J),J=1,3),(CONDUCTIVITY_TENSOR(I,J),J=1,3)

              IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
                SEED_VALUE(I) = 1000.0_DP
              ENDIF

            ENDDO
          ENDIF

          DO I=1,totalNumberOfNodes
!            CALL CALCULATE_SECOND_EIGENVECTOR()
            A=[CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)]
            B=[0.0_DP,0.0_DP,1.0_DP]
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CrossProduct(A,B,C,Err,Error,*999)
            ELSE
              B=[0.0_DP,1.0_DP,0.0_DP]
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ELSE
                B=[1.0_DP,0.0_DP,0.0_DP]
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=[CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)]
            CALL CrossProduct(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            CALL CALCULATE_SECOND_EIGENVECTOR()
          ENDDO

        ENDIF

      ENDIF

! CONNSDONCOCNCNOSKCN      load data for the CONNECTIVITY list 
      READ(11,*) STRING

      DO I=1,totalNumberOfNodes

        READ(11,*) STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
        
        DO J=1,3
          SPEED_FUNCTION_TABLE(I,J)=SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
        ENDDO
          
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))

      ENDDO

      CLOSE(11)

    ENDIF


! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-1

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDDO

      I=INT((totalNumberOfNodes/3.0_DP)+0.5_DP)

      IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 0) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3),(NODE_LIST(3*(I-1)+3,J),J=1,3)

      ENDIF
      IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 1) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3),(NODE_LIST(3*(I-1)+2,J),J=1,3)


      ENDIF
      IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 2) THEN

        READ(11,*) (NODE_LIST(3*(I-1)+1,J),J=1,3)

      ENDIF
      

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO
        
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO
        
      ENDDO
      
      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        DO I=1,totalNumberOfNodes
          STATUS_MASK(I) = ""
        ENDDO

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
!          PRINT*,(NODE_LIST(J+1,K),K=1,3),SEED_VALUE(J+1)      
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,totalNumberOfNodes 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,totalNumberOfNodes 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((totalNumberOfNodes/3.0_DP)+0.5_DP)

            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,totalNumberOfNodes 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,totalNumberOfNodes 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=[CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)]

            B=[0.0_DP,0.0_DP,1.0_DP]
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CrossProduct(A,B,C,Err,Error,*999)
            ELSE
              B=[0.0_DP,1.0_DP,0.0_DP]
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ELSE
                B=[1.0_DP,0.0_DP,0.0_DP]
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=[CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)]
            CALL CrossProduct(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 
        
        DO I=1,totalNumberOfNodes
          DO J=1,3
            SPEED_FUNCTION_TABLE(I,J) = 0.0_DP
          ENDDO
        ENDDO
        
          DO I=1,totalNumberOfNodes
!            DO J=1,CONNECTIVITY_NUMBER(I)
!              IF (CONNECTIVITY_LIST(I,J) > totalNumberOfNodes) THEN
!                PRINT*, I,J,CONNECTIVITY_NUMBER(I),CONNECTIVITY_LIST(I,J),totalNumberOfNodes
!              ENDIF
              SPEED_FUNCTION_TABLE(I,1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
              SPEED_FUNCTION_TABLE(I,2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
              SPEED_FUNCTION_TABLE(I,3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),1) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),2) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(2)
!              SPEED_FUNCTION_TABLE(I,CONNECTIVITY_LIST(I,J),3) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(3)
!            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
! """""""""""""""""""""""""""""""""""INPUT OF VTK FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "VTKTET1NPL") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"

      OPEN (11,FILE=STRING)
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING
      READ(11,*) STRING

      DO I=1,totalNumberOfNodes

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 

      READ(11,*) STRING

      DO N=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

! we have nodes starting number of 0
        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) N

        DO I=1,N 

          READ(11,*) J,SEED_VALUE(J+1)
       
          STATUS_MASK(J+1) = "SEED POINT"

        ENDDO

        CLOSE(11)

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,totalNumberOfNodes 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING
!          print *,STRING
        
        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,totalNumberOfNodes 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          IF (STRING .EQ. '#') THEN! begin if

            DO WHILE (STRING .NE. 'fiber')
      
              READ(11,*) STRING
        
            ENDDO
      
!      PRINT *,STRING
      
            DO I=1,INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-1

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDDO

            I=INT((totalNumberOfNodes/3.0_DP)+0.5_DP)

            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 0) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3), &
               & (CONDUCTIVITY_TENSOR(3*(I-1)+3,J),J=1,3)

            ENDIF
            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 1) THEN
    
              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3),(CONDUCTIVITY_TENSOR(3*(I-1)+2,J),J=1,3)

            ENDIF
            IF (3*INT((totalNumberOfNodes/3.0_DP)+0.5_DP)-totalNumberOfNodes .EQ. 2) THEN

              READ(11,*) (CONDUCTIVITY_TENSOR(3*(I-1)+1,J),J=1,3)

            ENDIF
      
          
          ELSE
            DO I=1,totalNumberOfNodes 
            
              READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3) 
              
            ENDDO               
          ENDIF ! end if
          
          DO I=1,totalNumberOfNodes 
          
!            READ(11,*) (CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=[CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)]
            B=[0.0_DP,0.0_DP,1.0_DP]
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CrossProduct(A,B,C,Err,Error,*999)
            ELSE
              B=[0.0_DP,1.0_DP,0.0_DP]
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ELSE
                B=[1.0_DP,0.0_DP,0.0_DP]
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=[CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)]
            CALL CrossProduct(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

! """""""""""""""""""""""""""""""""""INPUT OF CARP FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "CARP") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".pts"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      DO I=1,totalNumberOfNodes 

        READ(11,*) (NODE_LIST(I,J),J=1,3)

      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".elem"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT
          ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
        ENDDO

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,totalNumberOfNodes 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,totalNumberOfNodes 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".lon"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,totalNumberOfNodes 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR" .OR. INPUT_TYPE_FOR_CONDUCTIVITY .EQ. " ") THEN
          DO I=1,TOTAL_NUMBER_OF_ELEMENTS
            READ(11,*) (CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),J),J=1,3)

            A=[CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),1),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),2), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),3)]
            B=[0.0_DP,0.0_DP,1.0_DP]
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CrossProduct(A,B,C,Err,Error,*999)
            ELSE
              B=[0.0_DP,1.0_DP,0.0_DP]
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ELSE
                B=[1.0_DP,0.0_DP,0.0_DP]
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=[CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),4),CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),5), &
             & CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),6)]
            CALL CrossProduct(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>=ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF
            
            DO J=2,NUMBER_OF_NODES_PER_ELEMENT
              DO K=1,9
                CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,J),K)=CONDUCTIVITY_TENSOR(ELEMENT_LIST(I,1),K)
              ENDDO
            ENDDO

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF
    
! """""""""""""""""""""""""""""""""""INPUT OF TETGEN FORMAT""""""""""""""""""""""""""""""
    IF (INPUT_FILE_FORMAT .EQ. "TETGEN") THEN 

      NUMBER_OF_NODES_PER_ELEMENT = 4

! NDNAPDNDONOEEENODED      load data for NODAL POSITIONING info
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".node"

      OPEN (11,FILE=STRING)
      READ (11,*) STRING

      READ(11,*) FIRST_NODE_NUMBER,(NODE_LIST(1,J),J=1,3)

      DO I=2,totalNumberOfNodes 


        READ(11,*) STRING,(NODE_LIST(I,J),J=1,3)


      ENDDO

      CLOSE(11)

! EMMNLEMNTTMENLMNTMT      load data for ELEMENT CONNECTIVITY info 
      TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
      STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".ele"

      OPEN (11,FILE=STRING)
      READ (11,*) TOTAL_NUMBER_OF_ELEMENTS
      
      DO N=1,totalNumberOfNodes
        CONNECTIVITY_NUMBER(N)=0
      ENDDO

      DO I=1,TOTAL_NUMBER_OF_ELEMENTS

        READ(11,*) STRING,(ELEMENT_LIST(I,J),J=1,NUMBER_OF_NODES_PER_ELEMENT)

        IF (FIRST_NODE_NUMBER .EQ. 0) THEN        
          DO J=1,NUMBER_OF_NODES_PER_ELEMENT
            ELEMENT_LIST(I,J)=ELEMENT_LIST(I,J)+1
          ENDDO
        ENDIF

        DO J=1,NUMBER_OF_NODES_PER_ELEMENT-1
        
          DO K=1+J,NUMBER_OF_NODES_PER_ELEMENT
            THERE_IS_IN_CONNECTIVITY_LIST = 0
            DO N=1,CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))
!    print *,"elem",I,"connectivity number",N

              IF (CONNECTIVITY_LIST(ELEMENT_LIST(I,J),N) .EQ. ELEMENT_LIST(I,K)) THEN
                THERE_IS_IN_CONNECTIVITY_LIST = 1
              ENDIF

            ENDDO

            IF (THERE_IS_IN_CONNECTIVITY_LIST .EQ. 0) THEN

              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,J),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,J))) = ELEMENT_LIST(I,K)
              CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) = CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K)) + 1
              CONNECTIVITY_LIST(ELEMENT_LIST(I,K),CONNECTIVITY_NUMBER(ELEMENT_LIST(I,K))) = ELEMENT_LIST(I,J)

            ENDIF

          ENDDO

        ENDDO

      ENDDO

      CLOSE(11)

! SDEESDSEDSEESDSEEDS      load SEED VALUES data at the nods 
!      set input for *seed points* = LIST
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "LIST") THEN

        DO I=1,totalNumberOfNodes 

          IF (STATUS_MASK(I) .NE. "SEED POINT") THEN
            SEED_VALUE(I) = 1000.0_DP
          ENDIF

        ENDDO
        
      ENDIF

!      read input for *seed points* = FILE
      IF (INPUT_TYPE_FOR_SEED_VALUE .EQ. "FILE") THEN

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".estm"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        DO I=1,totalNumberOfNodes 
 
          READ(11,*) STRING,SEED_VALUE(I)

        ENDDO

        CLOSE(11)

      ENDIF

! NDCONDNCODNCODNDCOM      load CONDUCTIVITY TENSOR data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

          DO I=1,totalNumberOfNodes 
            CONDUCTIVITY_TENSOR(I,1)=1.0_DP
            CONDUCTIVITY_TENSOR(I,2)=0.0_DP
            CONDUCTIVITY_TENSOR(I,3)=0.0_DP
            CONDUCTIVITY_TENSOR(I,4)=0.0_DP
            CONDUCTIVITY_TENSOR(I,5)=1.0_DP
            CONDUCTIVITY_TENSOR(I,6)=0.0_DP
            CONDUCTIVITY_TENSOR(I,7)=0.0_DP
            CONDUCTIVITY_TENSOR(I,8)=0.0_DP
            CONDUCTIVITY_TENSOR(I,9)=1.0_DP
          ENDDO

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

        TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
        STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".fiber"

        OPEN (11,FILE=STRING)
        READ (11,*) STRING

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "TENSOR") THEN
          DO I=1,totalNumberOfNodes 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,9)
          ENDDO
        ENDIF

        IF (INPUT_TYPE_FOR_CONDUCTIVITY .EQ. "VECTOR") THEN
          DO I=1,totalNumberOfNodes 
            READ(11,*) STRING,(CONDUCTIVITY_TENSOR(I,J),J=1,3)

            A=[CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),CONDUCTIVITY_TENSOR(I,3)]
            B=[0.0_DP,0.0_DP,1.0_DP]
            CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)

            IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
              CALL CrossProduct(A,B,C,Err,Error,*999)
            ELSE
              B=[0.0_DP,1.0_DP,0.0_DP]
              CALL VECTOR_VECTOR_PRODUCT(A,B,DOT_PRODUCT_VALUE,Err)
              IF (DOT_PRODUCT_VALUE .LT. 0.8_DP) THEN
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ELSE

                B=[1.0_DP,0.0_DP,0.0_DP]
                CALL CrossProduct(A,B,C,Err,Error,*999)
              ENDIF
            ENDIF

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,4) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,5) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,6) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

            B=[CONDUCTIVITY_TENSOR(I,4),CONDUCTIVITY_TENSOR(I,5),CONDUCTIVITY_TENSOR(I,6)]
            CALL CrossProduct(A,B,C,Err,Error,*999)

            IF (ABS(SQRT(C(1)**2+C(2)**2+C(3)**2))>ZERO_TOLERANCE) THEN

              CONDUCTIVITY_TENSOR(I,7) = C(1)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,8) = C(2)/SQRT(C(1)**2+C(2)**2+C(3)**2)
              CONDUCTIVITY_TENSOR(I,9) = C(3)/SQRT(C(1)**2+C(2)**2+C(3)**2)

            ENDIF

!            PRINT*,(CONDUCTIVITY_TENSOR(I,J),J=1,3),DOT_PRODUCT_VALUE

          ENDDO
        ENDIF

        CLOSE(11)

      ENDIF

! FNCONSJDCNFUCNSUCNF      load VELOCITY FUNCTION data for the nods 
!      for the ISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(1)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1
          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,2) = SPEED_FUNCTION_TABLE(J,1)
            SPEED_FUNCTION_TABLE(J,3) = SPEED_FUNCTION_TABLE(J,1)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

!      for the ANISOTROPIC materials
      IF (MATERIAL_BEHAVIOUR .EQ. "ANISOTROPIC") THEN 

!        if it comes with FIXED format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FIXED") THEN 

          DO I=1,totalNumberOfNodes
            DO J=1,3
              SPEED_FUNCTION_TABLE(I,J) = SPEED_FUNCTION_ALONG_EIGEN_VECTOR(J)
            ENDDO
          ENDDO

        ENDIF

!        if it comes with FILE format
        IF (INPUT_TYPE_FOR_SPEED_FUNCTION .EQ. "FILE") THEN 

          TEXT_LENGTH = INDEX(INPUT_FILE_NAME,' ') - 1

          STRING = INPUT_FILE_NAME(1:TEXT_LENGTH)//".cond"

          OPEN (11,FILE=STRING)
          READ (11,*) N

          DO I=1,N 

            READ(11,*) J,SPEED_FUNCTION_TABLE(J,1),SPEED_FUNCTION_TABLE(J,2),&
                          &SPEED_FUNCTION_TABLE(J,3)

          ENDDO

          CLOSE(11)

        ENDIF

      ENDIF

    ENDIF

    ! to set RAW_INDEX and COLUMN_INDEX
    RAW_INDEX(1) = 0
    DO I=1,totalNumberOfNodes
        
      RAW_INDEX(I+1) = RAW_INDEX(I) + CONNECTIVITY_NUMBER(I)
      DO J = 1,CONNECTIVITY_NUMBER(I)
        COLUMN_INDEX(RAW_INDEX(I)+J) = CONNECTIVITY_LIST(I,J)
      ENDDO          

    ENDDO
    
    DO I=1,totalNumberOfNodes
        
      DO J = RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DO K = 1,3
          SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,K) = &
          & (SPEED_FUNCTION_TABLE(I,K)+SPEED_FUNCTION_TABLE(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
        DO K = 1,9
          CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,K) = &
          & (CONDUCTIVITY_TENSOR(I,K)+CONDUCTIVITY_TENSOR(COLUMN_INDEX(J),K))/2.0_DP
        ENDDO
      ENDDO

    ENDDO

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
999 ERRORS("GENERATE_STATUS_MASK",err,error)
    RETURN

  END SUBROUTINE PRE_PROCESS_INFORMATION


  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM(totalNumberOfNodes,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  &SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(totalNumberOfNodes,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,totalNumberOfNodes

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!    MAIN LOOP 
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=[  NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)]

!            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)/SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)
!            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE = [1.0_DP,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO], SHAPE = [3,3])
            CONDUCTION_RATIO= &
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=[1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO], SHAPE = [3,3])


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =[CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)], &
                               & SHAPE = [3,3])

            CALL MatrixTranspose(F,FT,Err,Error,*999)
            CALL MatrixProduct(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MatrixProduct(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

!	    PRINT *,F(1,1),F(1,2),F(1,3),F(2,1),F(2,2),F(2,3),F(3,1),F(3,2),F(3,3)

            CALL MatrixVectorProduct(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+SQRT(ABS(VMV))*&
          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)
!          &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)


            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,totalNumberOfNodes

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO

    RETURN
999 ERRORS("SOLVE_PROBLEM_FMM",err,error)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY(totalNumberOfNodes,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)  
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER, CHANGED_STATUS, MIN_TRIAL_NODE, TRIAL_STATUS
    REAL(DP), DIMENSION(3)   :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(2)   :: CONDUCTION_RATIO
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    !Start Program

    CALL GENERATE_STATUS_MASK(totalNumberOfNodes,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,totalNumberOfNodes

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=[NODE_LIST(COLUMN_INDEX(I),1)-NODE_LIST(COLUMN_INDEX(J),1)&
                            &,NODE_LIST(COLUMN_INDEX(I),2)-NODE_LIST(COLUMN_INDEX(J),2)&
                            &,NODE_LIST(COLUMN_INDEX(I),3)-NODE_LIST(COLUMN_INDEX(J),3)]

            CONDUCTION_RATIO(1)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
            CONDUCTION_RATIO(2)= SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=[1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)], SHAPE = [3,3])


!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =[CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                               &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)],SHAPE = [3,3])

            CALL MatrixTranspose(F,FT,Err,Error,*999)
            CALL MatrixProduct(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MatrixProduct(F,MFT,FMFT,Err,Error,*999)
!            CALL INVERT(FMFT,INV_FMFT,DET,Err,Error,*999)

            CALL MatrixVectorProduct(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,totalNumberOfNodes

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


    RETURN
999 ERRORS("SOLVE_PROBLEM_FMM_CONNECTIVITY",err,error)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_FMM_CONNECTIVITY


  !
  !================================================================================================================================
  !

  SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY(totalNumberOfNodes,NODE_LIST,CONDUCTIVITY_TENSOR_ON_CONNECTIVITY,&
                       &SPEED_FUNCTION_TABLE_ON_CONNECTIVITY,RAW_INDEX,COLUMN_INDEX,TOTAL_NUMBER_OF_CONNECTIVITY,&
                       &SEED_VALUE,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_CONNECTIVITY
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(:,:)
    INTEGER(INTG), ALLOCATABLE :: COLUMN_INDEX(:)
    INTEGER(INTG), ALLOCATABLE :: RAW_INDEX(:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    REAL(DP), ALLOCATABLE :: CONNECTIVITY_WEIGHT(:)
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW
    REAL(DP), DIMENSION(2) :: CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0
    
    !initialize data
    DO I=1,totalNumberOfNodes
      TRACE_NODE(I)=0
    ENDDO

    !Start Program
    ALLOCATE(CONNECTIVITY_WEIGHT(TOTAL_NUMBER_OF_CONNECTIVITY),STAT=ERR)
    DO I=1,totalNumberOfNodes
      DO J=RAW_INDEX(I)+1,RAW_INDEX(I+1)
        DISTANCE_VECTOR=[NODE_LIST(I,1)-NODE_LIST(COLUMN_INDEX(J),1)&
                        &,NODE_LIST(I,2)-NODE_LIST(COLUMN_INDEX(J),2)&
                        &,NODE_LIST(I,3)-NODE_LIST(COLUMN_INDEX(J),3)]

        CONDUCTION_RATIO(1) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,2)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
        CONDUCTION_RATIO(2) = SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,3)/SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)

        CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=[1.0_DP,0.0_DP,0.0_DP,&
                                            &0.0_DP,CONDUCTION_RATIO(1)*CONDUCTION_RATIO(1),0.0_DP,&
                                            &0.0_DP,0.0_DP,CONDUCTION_RATIO(2)*CONDUCTION_RATIO(2)], SHAPE = [3,3])

        F=RESHAPE(SOURCE =[CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,1),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,2),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,3),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,4),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,5),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,6),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,7),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,8),&
                           &CONDUCTIVITY_TENSOR_ON_CONNECTIVITY(J,9)],SHAPE = [3,3])

        CALL MatrixTranspose(F,FT,Err,Error,*999)
        CALL MatrixProduct(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
        CALL MatrixProduct(F,MFT,FMFT,Err,Error,*999)

        CALL MatrixVectorProduct(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
        CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

        CONNECTIVITY_WEIGHT(J)=SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE_ON_CONNECTIVITY(J,1)
      ENDDO
    ENDDO
          
    CALL GENERATE_STATUS_MASK(totalNumberOfNodes,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,totalNumberOfNodes

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

! ::::::: MAIN LOOP :::::::	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
!    LOOP_NUMBER = 0
    PRINT *,"Running GEODESIC Solver ...."

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=RAW_INDEX(MIN_TRIAL_NODE)+1,RAW_INDEX(MIN_TRIAL_NODE+1)
        TIME_NEW=1000

        DO J=RAW_INDEX(COLUMN_INDEX(I))+1,RAW_INDEX(COLUMN_INDEX(I)+1)
          IF (STATUS_MASK(COLUMN_INDEX(J)) == "KNOWN") THEN 

            TIME_ITER=SEED_VALUE(COLUMN_INDEX(J))+CONNECTIVITY_WEIGHT(J)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=COLUMN_INDEX(J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .LT. SEED_VALUE(COLUMN_INDEX(I))) THEN
          SEED_VALUE(COLUMN_INDEX(I)) = TIME_NEW
          TRACE_NODE(COLUMN_INDEX(I)) = TRACE_NODE_NUMBER
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "KNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(COLUMN_INDEX(I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(COLUMN_INDEX(I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000.0_DP
      CHANGED_STATUS = 0

      DO I=1,totalNumberOfNodes

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE = I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0 .AND. STATUS_MASK(I) .EQ. "SEED POINT" .AND. SEED_VALUE(I) .LT. MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO

    RETURN
999 ERRORS("SOLVE_PROBLEM_GEODESIC_CONNECTIVITY",err,error)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC_CONNECTIVITY

  !
  !================================================================================================================================
  !


  SUBROUTINE SOLVE_PROBLEM_GEODESIC(totalNumberOfNodes,NODE_LIST,CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,&
  & SEED_VALUE,CONNECTIVITY_NUMBER,CONNECTIVITY_LIST,STATUS_MASK,TRACE_NODE)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)
    INTEGER(INTG), ALLOCATABLE :: TRACE_NODE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)
    
    !Local Variables
    INTEGER(INTG) :: I,J
    INTEGER(INTG) :: LOOP_NUMBER,CHANGED_STATUS,MIN_TRIAL_NODE,TRIAL_STATUS,TRACE_NODE_NUMBER
    REAL(DP), DIMENSION(3) :: DISTANCE_VECTOR,MV
    REAL(DP), DIMENSION(3,3) :: CONDUCTIVITY_MATRIX,F,FT,MFT,FMFT
    REAL(DP) :: MIN_TRIAL_VALUE,VMV,MINIMUM_DATA,TIME_ITER,TIME_NEW,CONDUCTION_RATIO
    INTEGER(INTG) :: Err
    TYPE(VARYING_STRING) :: Error
    Err = 0

    !Start Program

    CALL GENERATE_STATUS_MASK(totalNumberOfNodes,SEED_VALUE,STATUS_MASK,Err)

    MIN_TRIAL_VALUE = 1000
    DO I=1,totalNumberOfNodes

      IF (STATUS_MASK(I) == "SEED POINT") THEN

        IF (SEED_VALUE(I) .lt. MIN_TRIAL_VALUE) THEN
          MIN_TRIAL_VALUE=SEED_VALUE(I)
          MIN_TRIAL_NODE = I
          TRIAL_STATUS = 1
        ENDIF

      ENDIF

    ENDDO

!   4444444444444444 MAIN LOOP 55555555555555555	P_NODE_NUMBER=CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)	PP_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
    LOOP_NUMBER = 0

    DO WHILE (TRIAL_STATUS .eq. 1) 

      TRIAL_STATUS = 0
      LOOP_NUMBER = LOOP_NUMBER + 1
      PRINT *,"Running in loop number",LOOP_NUMBER

      ! CALL ASSIGN_MIN_TRIAL_TO_KNOWN(STATUS_MASK,MIN_TRIAL_NODE)
      STATUS_MASK(MIN_TRIAL_NODE) = "KNOWN"

      DO I=1,CONNECTIVITY_NUMBER(MIN_TRIAL_NODE)
        TIME_NEW=1000

        DO J=1,CONNECTIVITY_NUMBER(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))
          IF (STATUS_MASK(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)) == "KNOWN") THEN 

            DISTANCE_VECTOR=[NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),1)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),2)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)&
                             &,NODE_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),3)&
                             &-NODE_LIST(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3)]

            CONDUCTION_RATIO=SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2)/&
                            &SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            CONDUCTIVITY_MATRIX=RESHAPE(SOURCE=[1.0_DP,0.0_DP,0.0_DP,&
                                                &0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO,0.0_DP,&
                                                &0.0_DP,0.0_DP,CONDUCTION_RATIO*CONDUCTION_RATIO], SHAPE = [3,3])

!            CALL LOAD_MATRIX((CONDUCTION_TENSOR(MIN_TRIAL_NODE,K),K=1,9),F)
            F=RESHAPE(SOURCE =[CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),2),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),3),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),4),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),5),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),6),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),7),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),8),&
                               &CONDUCTIVITY_TENSOR(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),9)], &
                               & SHAPE = [3,3])

            CALL MatrixTranspose(F,FT,Err,Error,*999)
            CALL MatrixProduct(CONDUCTIVITY_MATRIX,FT,MFT,Err,Error,*999)
            CALL MatrixProduct(F,MFT,FMFT,Err,Error,*999)

            CALL MatrixVectorProduct(FMFT,DISTANCE_VECTOR,MV,Err,Error,*999)
            CALL VECTOR_VECTOR_PRODUCT(DISTANCE_VECTOR,MV,VMV,Err)

            TIME_ITER=SEED_VALUE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J))+&
                      &SQRT(ABS(VMV))*SPEED_FUNCTION_TABLE(CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J),1)

            IF (TIME_ITER .lt. TIME_NEW) THEN
              TIME_NEW = TIME_ITER 
              TRACE_NODE_NUMBER=CONNECTIVITY_LIST(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I),J)
            ENDIF

          ENDIF
        ENDDO

        IF (TIME_NEW .lt. SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))) THEN
          SEED_VALUE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = TIME_NEW
          TRACE_NODE(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I))=TRACE_NODE_NUMBER
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "KNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "CHANGED"
          ENDIF
          IF (STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) .EQ. "UNKNOWN") THEN
            STATUS_MASK(CONNECTIVITY_LIST(MIN_TRIAL_NODE,I)) = "SEED POINT"
          ENDIF
        ENDIF

      ENDDO

      MINIMUM_DATA = 1000
      CHANGED_STATUS = 0

      DO I=1,totalNumberOfNodes

        IF (STATUS_MASK(I) .EQ. "CHANGED") THEN
          MIN_TRIAL_NODE=I
          CHANGED_STATUS = 1
          TRIAL_STATUS = 1
        ENDIF

        IF (CHANGED_STATUS .EQ. 0.AND.STATUS_MASK(I) .EQ. "SEED POINT".AND.SEED_VALUE(I).LT.MINIMUM_DATA) THEN
          MINIMUM_DATA = SEED_VALUE(I)
          MIN_TRIAL_NODE=I
          TRIAL_STATUS = 1
        ENDIF

      ENDDO

    ENDDO


    RETURN
999 ERRORS("SOLVE_PROBLEM_GEODESIC",err,error)
    RETURN

  END SUBROUTINE SOLVE_PROBLEM_GEODESIC

  
  !
  !================================================================================================================================
  !

  !>Calculates status mask for the local nodes.
  SUBROUTINE GENERATE_STATUS_MASK(totalNumberOfNodes,SEED_VALUE,STATUS_MASK,Err)

    !subroutine variables
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    CHARACTER (LEN=10), ALLOCATABLE :: STATUS_MASK(:)

    !Argument variables
    INTEGER(INTG) :: Err !<The error code
    !    TYPE(VARYING_STRING) :: localError !<The error string

    !Local variables
    INTEGER(INTG) :: I
        
!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)
    
    DO I=1,totalNumberOfNodes

      IF (SEED_VALUE(I) .LT. 100.0_DP) THEN
        STATUS_MASK(I) = "SEED POINT"
      ELSE
        STATUS_MASK(I) = "UNKNOWN"
      ENDIF

    ENDDO

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 ERRORSEXITS("GENERATE_STATUS_MASK",err,error)
!    RETURN 1

  END SUBROUTINE GENERATE_STATUS_MASK


  !
  !================================================================================================================================
  !

  !>Calculates minimum and maximum value at array A.
  SUBROUTINE FIND_MINIMAX(A,N,MIN_VALUE,MAX_VALUE,Err)

    !Argument variables
    REAL(DP), ALLOCATABLE :: A(:)

    REAL(DP), INTENT(OUT) :: MIN_VALUE
    REAL(DP), INTENT(OUT) :: MAX_VALUE
    INTEGER(INTG), INTENT(IN)  :: N
    INTEGER(INTG) :: err !<The error code
    !    TYPE(VARYING_STRING) :: localError !<The error string
    !Local variables
    INTEGER(INTG) :: I
        
!    ENTERS("FIND_MINIMAX",Err,Error,*999)
    
    IF(SIZE(A,1).GT.2) THEN
      MIN_VALUE = A(1)
      MAX_VALUE = A(1)
      DO I=2,N
        IF (MIN_VALUE .GT. A(I)) THEN
          MIN_VALUE = A(I)
        ENDIF
        IF (MAX_VALUE .LT. A(I)) THEN
          MAX_VALUE = A(I)
        ENDIF
      ENDDO
!    ELSE
!      CALL FlagError("Invalid matrix sizes.",Err)
    ENDIF

!    EXITS("FIND_MINIMAX")
!    RETURN
!999 ERRORSEXITS("FIND_MINIMAX",err,error)
!    RETURN 1

  END SUBROUTINE FIND_MINIMAX

  !
  !================================================================================================================================
  !


  !>Calculates and returns the VECTOR-VECTOR-prouct of the double precision VECTOR A*B in C.
  SUBROUTINE VECTOR_VECTOR_PRODUCT(A,B,C,Err)

    !Argument variables
    REAL(DP), INTENT(IN) :: A(3) !<The A VECTOR
    REAL(DP), INTENT(IN) :: B(3) !<The B VECTOR
    REAL(DP), INTENT(OUT) :: C !<On exit, the product SCALAR C=A*B
    INTEGER(INTG) :: err !<The error code
    !    TYPE(VARYING_STRING) :: localError !<The error string
    !Local variables
        
!    ENTERS("VECTOR_VECTOR_PRODUCT",Err,Error,*999)
    
    IF(SIZE(A,1)==SIZE(B,1)) THEN
      SELECT CASE(SIZE(A,1))
        CASE(1)
          C=A(1)*B(1)
        CASE(2)
          C=A(1)*B(1)+A(2)*B(2)
        CASE(3)
          C=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
!        CASE DEFAULT
!          CALL FlagError("Invalid matrix size.",Err)
      END SELECT
!    ELSE
!      CALL FlagError("Invalid matrix sizes.",Err)
    ENDIF

!    EXITS("VECTOR_VECTOR_PRODUCT")
!    RETURN
!999 ERRORSEXITS("VECTOR_VECTOR_PRODUCT",err,error)
!    RETURN 1

  END SUBROUTINE VECTOR_VECTOR_PRODUCT

  !
  !================================================================================================================================
  !


  !>to EXPORT output.
  SUBROUTINE POST_PROCESS_DATA(MATERIAL_BEHAVIOUR,OUTPUT_FILE_NAME,OUTPUT_FILE_FORMAT,totalNumberOfNodes,NODE_LIST,&
&CONDUCTIVITY_TENSOR,SPEED_FUNCTION_TABLE,SEED_VALUE,CONNECTIVITY_NUMBER,OUTPUT_FILE_FIELD_TITLE,&
&CONNECTIVITY_LIST,ELEMENT_LIST,TOTAL_NUMBER_OF_ELEMENTS,NUMBER_OF_NODES_PER_ELEMENT,Err)

    !subroutine variables
    REAL(DP), ALLOCATABLE :: NODE_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SPEED_FUNCTION_TABLE(:,:)
    REAL(DP), ALLOCATABLE :: CONDUCTIVITY_TENSOR(:,:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_LIST(:,:)
    INTEGER(INTG), ALLOCATABLE :: ELEMENT_LIST(:,:)
    REAL(DP), ALLOCATABLE :: SEED_VALUE(:)
    INTEGER(INTG), ALLOCATABLE :: CONNECTIVITY_NUMBER(:)

    INTEGER(INTG), INTENT(IN) :: TOTAL_NUMBER_OF_ELEMENTS
    INTEGER(INTG), INTENT(IN) :: NUMBER_OF_NODES_PER_ELEMENT
    INTEGER(INTG), INTENT(IN) :: totalNumberOfNodes
    CHARACTER (LEN=300) :: OUTPUT_FILE_NAME
    CHARACTER (LEN=300) :: OUTPUT_FILE_FIELD_TITLE
    CHARACTER (LEN=10)  :: OUTPUT_FILE_FORMAT
    CHARACTER (LEN=12)  :: MATERIAL_BEHAVIOUR
    INTEGER(INTG) :: Err
    
    !Local variables
    CHARACTER (LEN=300) :: STRING
    INTEGER(INTG) :: TEXT_LENGTH
    INTEGER(INTG) :: I, J

!    ENTERS("GENERATE_STATUS_MASK",Err,Error,*999)

!    IF (OUTPUT_FILE_FORMAT .NE. "TABC") THEN
!      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,">> OUTPUT FILE FORMAT ERROR: identify a correct output format",Err)
!    ENDIF
    

! in the case the OUTPUT is in TABC format (2)
    IF (OUTPUT_FILE_FORMAT .EQ. "TABC") THEN 

!      EXPORT NODE TABLE list
      TEXT_LENGTH = INDEX(OUTPUT_FILE_NAME,' ') - 1
      STRING = OUTPUT_FILE_NAME(1:TEXT_LENGTH)//".tabc"
!  INQUIRE(FILE=STRING, EXIST=ex)
!      PRINT *, STRING
      OPEN (12,FILE=STRING)
!      OPEN (12,FILE=OUTPUT_FILE_NAME)

      WRITE(12,*)"VARIABLES=""NODE"",""X"",""Y"",""Z"",""U"",""V"",""W"",""Speed function along fibers"", &
 &                ""Speed function in transverse direction"",""Time"""
      WRITE(12,*)"zone i=",totalNumberOfNodes," , DATAPACKING=POINT"

      DO I=1,totalNumberOfNodes
        WRITE(12,*) I,NODE_LIST(I,1),NODE_LIST(I,2),NODE_LIST(I,3),CONDUCTIVITY_TENSOR(I,1),CONDUCTIVITY_TENSOR(I,2),&
            &CONDUCTIVITY_TENSOR(I,3),SPEED_FUNCTION_TABLE(I,1),&
            &SPEED_FUNCTION_TABLE(I,2),SPEED_FUNCTION_TABLE(I,3),&
            &SEED_VALUE(I)
      ENDDO
!      EXPORT NODE CONNECTIVITY list
      WRITE(12,*) "Connectivity"
      DO I=1,totalNumberOfNodes
        WRITE(12,*) I,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
!        PRINT *,STRING,CONNECTIVITY_NUMBER(I),(CONNECTIVITY_LIST(I,J),J=1,CONNECTIVITY_NUMBER(I))
      ENDDO
      CLOSE(12)

    ENDIF

! in the case the OUTPUT is in VTK TETrahedral format (2)
    IF (OUTPUT_FILE_FORMAT .EQ. "VTKTET") THEN 

!      rename to VTK and OPEN the file
      TEXT_LENGTH = INDEX(OUTPUT_FILE_NAME,' ') - 1
      STRING = OUTPUT_FILE_NAME(1:TEXT_LENGTH)//".vtk"
      OPEN (12,FILE=STRING)

!      export HEADER list
      WRITE(12,'(A)')"# vtk DataFile Version 3.0"
      WRITE(12,'(A)')"vtk output"
      WRITE(12,'(A)')"ASCII"
      WRITE(12,'(A)')"DATASET UNSTRUCTURED_GRID"
      WRITE(12,'(A,I8,A6)')"POINTS",totalNumberOfNodes,"float"

!      export NODAL POSITION list
      DO I=1,totalNumberOfNodes
        WRITE(12,*) (NODE_LIST(I,J),J=1,3)
      ENDDO

!      export ELEMENT list
      WRITE(12,'(A,I8,I8)')"CELLS",TOTAL_NUMBER_OF_ELEMENTS,TOTAL_NUMBER_OF_ELEMENTS*(NUMBER_OF_NODES_PER_ELEMENT+1)
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        WRITE(12,*) NUMBER_OF_NODES_PER_ELEMENT,((ELEMENT_LIST(I,J)-1),J=1,NUMBER_OF_NODES_PER_ELEMENT)
      ENDDO

!      export ELEMENT TYPE list
      WRITE(12,'(A,I8)')"CELL_TYPES",TOTAL_NUMBER_OF_ELEMENTS
      DO I=1,TOTAL_NUMBER_OF_ELEMENTS
        WRITE(12,'(A)') "10"
      ENDDO

!      export CELL and POINT DATA list
      WRITE(12,'(A,I8)')"CELL_DATA",TOTAL_NUMBER_OF_ELEMENTS
      WRITE(12,'(A,I8)')"POINT_DATA",totalNumberOfNodes

!      export FIELD information
      WRITE(12,'(A,A)')"FIELD number"," 1"
      WRITE(12,'(A,I3,I8,A6)')OUTPUT_FILE_FIELD_TITLE,1,totalNumberOfNodes,"float"
      DO I=1,totalNumberOfNodes
        WRITE(12,'(F15.10)') SEED_VALUE(I)
      ENDDO

!      export VECTORS information
      WRITE(12,'(A,A,A6)') "VECTORS ","fiber_vector","float"
      DO I=1,totalNumberOfNodes
        WRITE(12,'(3F8.5)') (CONDUCTIVITY_TENSOR(I,J),J=1,3)
      ENDDO

      CLOSE(12)

    ENDIF

!   FILE="cmgui"
!   METHOD="FORTRAN"

!   EXPORT_FIELD=.TRUE.
!   IF(EXPORT_FIELD) THEN
!     WRITE(*,*)'Now export fields...'
!     CALL FLUID_MECHANICS_IO_WRITE_CMGUI(REGION,FILE,Err)
!     CALL FIELD_IO_NODES_EXPORT(REGION%FIELDS, FILE, METHOD, Err)  
!     CALL FIELD_IO_ELEMENTS_EXPORT(REGION%FIELDS, FILE, METHOD, Err)
!     WRITE(*,*)'All fields exported...'
!   ENDIF

!    EXITS("GENERATE_STATUS_MASK")
!    RETURN
!999 ERRORSEXITS("GENERATE_STATUS_MASK",err,error)
!    RETURN 1

  END SUBROUTINE POST_PROCESS_DATA

 
END MODULE HamiltonJacobiRoutines
