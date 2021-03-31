!> \file
!> \author Chris Bradley
!> \brief This module handles all Laplace equations routines.
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

!>This module handles all Laplace equations routines.
MODULE LaplaceEquationsRoutines

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

  PUBLIC Laplace_BoundaryConditionsAnalyticCalculate

  PUBLIC Laplace_EquationsSetSetup
  
  PUBLIC Laplace_EquationsSetSolutionMethodSet

  PUBLIC Laplace_EquationsSetSpecificationSet
  
  PUBLIC Laplace_FiniteElementCalculate

  PUBLIC Laplace_ProblemSetup

  PUBLIC Laplace_ProblemSpecificationSet


CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Laplace_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentInterpolationType,componentIdx,derivativeIdx,dimensionIdx, &
      & globalDerivativeIndex,localDOF,nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfVariables,variableIdx,variableType
    REAL(DP) :: value,x(3)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(VARYING_STRING) :: localError    
    
    ENTERS("Laplace_BoundaryConditionsAnalyticCalculate",err,error,*999)

    NULLIFY(geometricParameters)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
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
        CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,componentInterpolationType,err,error,*999)
        IF(componentInterpolationType/=FIELD_NODE_BASED_INTERPOLATION) &
          & CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)         
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
          DO dimensionIdx=1,numberOfDimensions
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOF,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOF)
          ENDDO !dimensionIdx
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          !Loop over the derivatives
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
              !u=x^2+2.x.y-y^2
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=x(1)*x(1)-2.0_DP*x(1)*x(2)-x(2)*x(2)
                CASE(GLOBAL_DERIV_S1)
                  value=2.0_DP*x(1)+2.0_DP*x(2)
                CASE(GLOBAL_DERIV_S2)
                  value=2.0_DP*x(1)-2.0_DP*x(2)
                CASE(GLOBAL_DERIV_S1_S2)
                  value=2.0_DP
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=0.0_DP !!TODO
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
              !u=cos(x).cosh(y)
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=COS(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S1)
                  value=-SIN(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S2)
                  value=COS(x(1))*SINH(x(2))
                CASE(GLOBAL_DERIV_S1_S2)
                  value=-SIN(x(1))*SINH(x(2))
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=0.0_DP !!TODO
                CASE(GLOBAL_DERIV_S1)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  !CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
              !u=x^2+y^2-2.z^2
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=x(1)*x(1)+x(2)*x(2)-2.0_DP*x(3)*x(3)
                CASE(GLOBAL_DERIV_S1)
                  value=2.0_DP*x(1)
                CASE(GLOBAL_DERIV_S2)
                  value=2.0_DP*x(2)
                CASE(GLOBAL_DERIV_S1_S2)
                  value=0.0_DP
                CASE(GLOBAL_DERIV_S3)
                  value=-4.0_DP*x(3)
                CASE(GLOBAL_DERIV_S1_S3)
                  value=0.0_DP
                CASE(GLOBAL_DERIV_S2_S3)
                  value=0.0_DP
                CASE(GLOBAL_DERIV_S1_S2_S3)
                  value=0.0_DP
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=0.0_DP !!TODO
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
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
              !u=cos(x).cosh(y).z
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=COS(x(1))*COSH(x(2))*x(3)
                CASE(GLOBAL_DERIV_S1)
                  value=-SIN(x(1))*COSH(x(2))*x(3)
                CASE(GLOBAL_DERIV_S2)
                  value=COS(x(1))*SINH(x(2))*x(3)
                CASE(GLOBAL_DERIV_S1_S2)
                  value=-SIN(x(1))*SINH(x(2))*x(3)
                CASE(GLOBAL_DERIV_S3)
                  value=COS(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S1_S3)
                  value=-SIN(x(1))*COSH(x(2))
                CASE(GLOBAL_DERIV_S2_S3)
                  value=COS(x(1))*SINH(x(2))
                CASE(GLOBAL_DERIV_S1_S2_S3)
                  value=-SIN(x(1))*SINH(x(2))
                CASE DEFAULT
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  value=0.0_DP !!TODO
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
                  localError="The global derivative index of "// &
                    & TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The analytic function type of "// &
                & TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE, &
              & 1,derivativeIdx,nodeIdx,componentIdx,VALUE,err,error,*999)
            IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
              CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
              IF(boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
                CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOF,BOUNDARY_CONDITION_FIXED, &
                  & VALUE,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    
    EXITS("Laplace_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Laplace_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Laplace_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Laplace equation finite element equations set.
  SUBROUTINE Laplace_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN)  :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,colsVariableType, &
      & dependentVariableType,esSpecification(3),gaussPointIdx,numberOfColumnComponents,numberOfColumnElementParameters, &
      & numberOfDimensions,numberOfDependentXi,numberOfGauss,numberOfGeometricXi,numberOfRowComponents, &
      & numberOfRowElementParameters,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowsVariableType,rowXiIdx, &
      & scalingType,xiIdx
    REAL(DP) :: columndPhidXi(3),conductivityMaterial(3,3),conductivity(3,3),conductivityTemp(3,3),dNudXi(3,3),dXidNu(3,3), &
      & dXdNu(3,3),dNudX(3,3),gaussWeight,jacobian,jacobianGaussWeight,sum,kValue(3),rowdPhidXi(3)
    LOGICAL :: update,updateMatrix,updateRHS
    TYPE(BasisType), POINTER :: columnBasis,dependentBasis,geometricBasis,rowBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,independentField,materialsField,fibreField
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & independentInterpParameters,materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,independentInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(QuadratureSchemeType), POINTER :: columnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError    

    ENTERS("Laplace_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    !Get the fields and check the equations
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    NULLIFY(fibreField)
    NULLIFY(independentField)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
      !Do nothing
    CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
    CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Laplace equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(equationsInterpolation)
    CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
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
      CALL Basis_NumberOfXiGet(geometricBasis,numberOfGeometricXi,err,error,*999)

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
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColumnComponents,err,error,*999)

      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPoint,err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)

      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      NULLIFY(materialsInterpParameters)
      NULLIFY(materialsInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpParameters,err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & fibreInterpPoint,err,error,*999)
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpPoint,err,error,*999)
      ENDIF

      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
      ENDIF

      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)

      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF
      IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters,err,error,*999)
      ENDIF

      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss

        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfGeometricXi,geometricInterpPointMetrics,err,error,*999)

        IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
            & err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)

          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & materialsInterpPoint%values(1:NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)        
        ENDIF
        IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
            & err,error,*999)
          kValue(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
        ENDIF

        !Calculate Jacobian and Gauss weight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight

        !Loop over field components
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
          NULLIFY(rowQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            DO rowXiIdx=1,numberOfDependentXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,rowdPhidXi(rowXiIdx),err,error,*999)
            ENDDO !rowXiIdx
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColumnComponents
                NULLIFY(columnDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
                NULLIFY(columnDomainTopology)
                CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                NULLIFY(columnDomainElements)
                CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                NULLIFY(columnBasis)
                CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis,err,error,*999)
                NULLIFY(columnQuadratureScheme)
                CALL Basis_QuadratureSchemeGet(columnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,columnQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(columnBasis,numberOfColumnElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColumnElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  DO columnXiIdx=1,numberOfDependentXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,columndPhidXi(columnXiIdx), &
                      & err,error,*999)
                  ENDDO !columnXiIdx
                  sum=0.0_DP
                  SELECT CASE(esSpecification(3))
                  CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
                    DO rowXiIdx=1,numberOfDependentXi
                      DO columnXiIdx=1,numberOfDependentXi
                        sum=sum+rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                          & geometricInterpPointMetrics%gu(rowXiIdx,columnXiIdx)
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
                    DO rowXiIdx=1,numberOfDependentXi
                      DO columnXiIdx=1,numberOfDependentXi
                        DO xiIdx=1,numberOfGeometricXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
                    IF(rowComponentIdx==columnComponentIdx) THEN 
                      DO rowXiIdx=1,numberOfDependentXi
                        DO columnXiIdx=1,numberOfDependentXi
                          sum=sum+kValue(rowComponentIdx)*rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,columnXiIdx)
                        ENDDO !columnXiIdx
                      ENDDO !rowXiIdx
                    ENDIF
                  CASE DEFAULT
                    localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                      & " is not valid for a Laplace equation type of a classical field equations set class."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT

                  linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & linearMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                    & sum*jacobianGaussWeight

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
        !Loop over element rows
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
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColumnComponents
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
      ENDIF ! scale factors

    ENDIF !update

    EXITS("Laplace_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Laplace_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Laplace_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up Laplace equations.
  SUBROUTINE Laplace_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables %
    INTEGER(INTG) :: componentIdx,componentNumber,esSpecification(3),geometricMeshComponent,geometricScalingType, &
      & numberOfDependentComponents,numberOfDimensions,numberOfIndependentComponents,numberOfMaterialsComponents,solutionMethod
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetIndependentType), POINTER :: equationsIndependent
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: componentLabel,localError
    
    ENTERS("Laplace_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "// &
        & TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
        & " is not valid for a Laplace type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
    
    SELECT CASE(equationsSetSetup%setupType)
    CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
      !
      ! Initial setup
      !
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL Laplace_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE)
        CALL EquationsSet_LabelSet(equationsSet,"Standard Laplace equations set",err,error,*999)       
      CASE(EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE)
        CALL EquationsSet_LabelSet(equationsSet,"Generalised Laplace equations set",err,error,*999)       
      CASE(EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
        CALL EquationsSet_LabelSet(equationsSet,"Moving mesh Laplace equations set",err,error,*999)               
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
          & " is not valid for a Laplace type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !
      ! Geometric setup
      !
      ! Do nothing
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !
      ! Dependent setup
      !
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
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
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition, &
            & err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,2,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
            numberOfDependentComponents=numberOfDimensions
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          ELSE
            numberOfDependentComponents=1
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          ENDIF          
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          DO componentIdx=1,numberOfDependentComponents
            componentLabel="Phi"//TRIM(NumberToVString(componentIdx,"*",err,error))
            CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,componentLabel,err,error,*999)
            componentLabel="del Phi"//TRIM(NumberToVString(componentIdx,"*",err,error))//"/del n"
            CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,componentLabel,err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
            localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
            numberOfDependentComponents=numberOfDimensions
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
          ELSE
            numberOfDependentComponents=1
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
          ENDIF
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          SELECT CASE(equationsSet%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                & numberOfDependentComponents,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDependentComponents,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
            localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) &
          & CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation"
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !
      ! Independent setup
      !
      NULLIFY(equationsIndependent)
      IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
      ENDIF      
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          IF(equationsSet%INDEPENDENT%independentFieldAutoCreated) THEN
            !Create the auto created independent field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsIndependent%independentField,err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            !define new created field to be independent
            CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,geometricField,err,error,*999)
            !set number of variables to 1 (1 for U)
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            !calculate number of components with one component for each dimension
            numberOfIndependentComponents=numberOfDimensions
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfIndependentComponents,err,error,*999)
            DO componentIdx=1,numberOfIndependentComponents
              CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx, & 
                & geometricMeshComponent,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField, & 
                & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999) 
            ENDDO !componentIdx
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !calculate number of components with one component for each dimension and one for pressure
            numberOfIndependentComponents=numberOfDimensions
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
              & numberOfIndependentComponents,err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(esSpecification(3)==EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE) THEN
          IF(equationsIndependent%independentFieldAutoCreated) &
            & CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
          CALL Field_ParameterSetEnsureCreated(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, & 
            & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          !Default k values to 1.0
          CALL Field_NumberOfComponentsGet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfIndependentComponents,err,error,*999)
          DO componentIdx=1,numberOfIndependentComponents
            CALL Field_ComponentValuesInitialise(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Laplace equations set"
        CALL FlagError(localError,err,error,*999)
      END SELECT      
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !
      ! Materials setup
      !
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      NULLIFY(equationsMaterials)
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
        CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      ENDIF
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%region,equationsMaterials%materialsField, &
              & err,error,*999)
            CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Conductivity",err,error,*999)
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)
            !Set the number of materials components
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,err,error,*999)        
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            !Default the materials components to the first geometric component
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfMaterialsComponents                
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !components_idx            
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,numberOfMaterialsComponents                
                CALL Field_ComponentInterpolationSetAndLock(equationsMaterials%materialsField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              !Default the scaling to the geometric field scaling
              CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
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
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
              & err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !do nothing
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
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) THEN
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
            !Default conductivity values to 1.0 on the diagonal, 0.0 elsewhere
            CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
              & numberOfDimensions,err,error,*999)
            DO componentIdx=1,numberOfDimensions
              componentNumber=TENSOR_TO_VOIGT(componentIdx,componentIdx,numberOfDimensions)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentNumber,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDIF
        ENDIF
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !
      ! Source setup
      !
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Do nothing
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !
      ! Analytic setup
      !
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(equationsSetSetup%analyticFunctionType)
        CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1)
          !Check that we are in 2D
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
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_1
        CASE(EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2)
          !Check that we are in 2D
          IF(numberOfDimensions/=2) THEN
            localError="The number of geometric dimensions of "// &
              & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid. The analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " requires that there be 2 geometric dimensions."
            CALL FlagError(localError,err,error,*999)
          ENDIF
          !Create analytic field if required
          !Set analytic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_LAPLACE_EQUATION_TWO_DIM_2
        CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1)
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
          !Set analytic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_1
        CASE(EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2)
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
          !Set analytic function type
          equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_LAPLACE_EQUATION_THREE_DIM_2
        CASE DEFAULT
          localError="The specified analytic function type of "// &
            & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
            & " is invalid for a standard Laplace equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
        IF(equationsSet%analytic%analyticFieldAutoCreated) &
          & CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !
      ! Equations setup
      !
      CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE) &
        & CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the equations creation
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
            & err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          SELECT CASE(equations%sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
              & err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
              & err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
              & err,error,*999)
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
          localError="The solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Laplace equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Laplace_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Laplace_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Laplace_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Laplace equation type of an classical field equations set class.
  SUBROUTINE Laplace_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Laplace_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
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
        & " is not valid for a Laplace equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Laplace_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Laplace_EquationsSetSolutionMethodSet",err,error)
    EXITS("Laplace_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Laplace_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Laplace type of a classical field equations set.
  SUBROUTINE Laplace_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Laplace_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    IF(SIZE(specification,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a Laplace type equations set.", &
      & err,error,*999)
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_LAPLACE_SUBTYPE, &
      & EQUATIONS_SET_MOVING_MESH_LAPLACE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a Laplace type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_LAPLACE_EQUATION_TYPE,subtype]

    EXITS("Laplace_EquationsSetSpecificationSet")
    RETURN
999 ERRORSEXITS("Laplace_EquationsSetSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Laplace_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Laplace problem.
  SUBROUTINE Laplace_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
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
    
    ENTERS("Laplace_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(pSpecification(3),"*",err,error))// &
        & " is not valid for a Laplace type of a classical field problem."
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
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      !
      ! Control loop setup
      !
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
          & " is invalid for a standard Laplace problem."
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
          & " is invalid for a standard Laplace equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !
      ! Solver equations setup
      !
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
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
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
          & " is invalid for a standard Laplace problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Laplace problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
        
    EXITS("Laplace_ProblemSetup")
    RETURN
999 ERRORSEXITS("Laplace_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Laplace_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Laplace equation type.
  SUBROUTINE Laplace_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Laplace_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_STANDARD_LAPLACE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Laplace type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_LAPLACE_EQUATION_TYPE,PROBLEM_STANDARD_LAPLACE_SUBTYPE]
 
    EXITS("Laplace_ProblemSpecificationSet")
    RETURN
999 ERRORS("Laplace_ProblemSpecificationSet",err,error)
    EXITS("Laplace_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Laplace_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
  
END MODULE LaplaceEquationsRoutines
