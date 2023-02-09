!> \file
!> \author Chris Bradley
!> \brief This module handles all Helmholtz equations routines.
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

!>This module handles all Helmholtz equations routines.
MODULE HelmholtzEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines
  USE CoordinateSystemAccessRoutines
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

  PUBLIC Helmholtz_BoundaryConditionsAnalyticCalculate
  
  PUBLIC Helmholtz_EquationsSetSolutionMethodSet
  
  PUBLIC Helmholtz_EquationsSetSetup

  PUBLIC Helmholtz_EquationsSetSpecificationSet
  
  PUBLIC Helmholtz_FiniteElementCalculate

  PUBLIC Helmholtz_ProblemSetup

  PUBLIC Helmholtz_ProblemSpecificationSet

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Helmholtz_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,derivativeIdx,dimensionIdx,globalDerivativeIndex,localDOFIdx, &
      & nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodes,numberOfNodeDerivatives,numberOfVariables,variableIdx, &
      & variableType
    REAL(DP) :: analyticValue,x(3),k,mu
    REAL(DP), POINTER :: geometricParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable,materialsVariable
    TYPE(VARYING_STRING) :: localError    
    
    ENTERS("Helmholtz_BoundaryConditionsAnalyticCalculate",err,error,*999)

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
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(materialsVariable)
    CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
    NULLIFY(materialsParameters)
    CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        CALL FieldVariable_ConstantDOFGet(materialsVariable,1,localDOFIdx,err,error,*999)
        k=materialsParameters(localDOFIdx)
        mu=k/SQRT(2.0_DP)
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
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1)
              !u=cos(k*x/sqrt(2))*sin(k*y/sqrt(2))
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(globalDerivativeIndex)
                CASE(NO_GLOBAL_DERIV)
                  analyticValue=COS(mu*x(1))*SIN(mu*x(2))
                CASE(GLOBAL_DERIV_S1)
                  analyticValue=-mu*SIN(mu*x(1))*SIN(mu*x(2))
                CASE(GLOBAL_DERIV_S2)
                  analyticValue=mu*COS(mu*x(1))*COS(mu*x(2))
                CASE(GLOBAL_DERIV_S1_S2)
                  analyticValue=-mu*mu*SIN(mu*x(1))*COS(mu*x(2))
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
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,1,derivativeIdx, &
              & nodeIdx,componentIdx,analyticValue,err,error,*999)
            IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
              IF(boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
                CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDOFIdx, &
                  & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
              ENDIF
            ENDIF
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    CALL FieldVariable_ParameterSetDataRestore(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    
    EXITS("Helmholtz_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("Helmholtz_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("Helmholtz_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Helmholtz_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Helmholtz equation finite element equations set.
  SUBROUTINE Helmholtz_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,colsVariableType,columnXiIdx, &
      & esSpecification(3),gaussPointIdx,numberOfColComponents,numberOfColElementParameters,numberOfColumnXi, &
      & numberOfDependentXi,numberOfDimensions,numberOfGauss,numberOfGeometricXi,numberOfMaterialsComponents, &
      & numberOfRowComponents,numberOfRowElementParameters,numberOfXi,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx,rowsVariableType,rowXiIdx,scalingType,xiIdx,variableType
    REAL(DP) :: columndPhidXi(3),columnPhi,conductivityTensor(3,3),gaussWeight,jacobian,jacobianGaussWeight,k,rowdPhidXi(3), &
      & rowPhi,sum
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
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: fibreInterpPoint,geometricInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,materialsVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: colQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme, &
      & rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("Helmholtz_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Helmholtz equation type of a classical field equations set class."
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
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(linearMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,linearMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
    CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)

    update=(updateMatrix.OR.updateRHS)

    IF(update) THEN

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
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
      
      NULLIFY(geometricInterpParameters)
      CALL EquationsInterpolation_GeometricParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpParameters,err,error,*999)
      NULLIFY(geometricInterpPoint)
      CALL EquationsInterpolation_GeometricPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPoint,err,error,*999)
      NULLIFY(geometricInterpPointMetrics)
      CALL EquationsInterpolation_GeometricPointMetricsGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & geometricInterpPointMetrics,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpParameters,err,error,*999)
          
      NULLIFY(fibreField)
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE) THEN
        CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
        IF(ASSOCIATED(fibreField)) THEN
          CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
            & err,error,*999)
          CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
        ENDIF
      ENDIF

      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(materialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(materialsVariable,numberOfMaterialsComponents,err,error,*999)
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
    
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfGeometricXi,geometricInterpPointMetrics,err,error,*999)        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        
        k = materialsInterpPoint%values(1,NO_PART_DERIV) ! wave number
        IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)
          CALL CoordinateSystem_MaterialTransformVoigtTensor2([COORDINATE_CONTRAVARIANT_INDEX_TYPE, &
            & COORDINATE_COVARIANT_INDEX_TYPE],geometricInterpPointMetrics,fibreInterpPoint, &
            & materialsInterpPoint%values(2:1+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivityTensor, &
            & err,error,*999)
        ENDIF
        
         !Calculate jacobianGaussWeight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
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
          CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowPhi,err,error,*999)
            DO rowXiIdx=1,numberOfDependentXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowXiIdx),gaussPointIdx,rowdPhidXi(rowXiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(updateMatrix) THEN
              columnElementDOFIdx=0
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
                CALL Basis_QuadratureSchemeGet(colBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colBasis,numberOfColElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(colQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & gaussPointIdx,columnPhi,err,error,*999)
                  DO columnXiIdx=1,numberOfDependentXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(columnXiIdx),gaussPointIdx,columndPhidXi(columnXiIdx), &
                      & err,error,*999)
                  ENDDO !columnXiIdx
                  sum=0.0_DP
                  sum=sum-k*k*rowPhi*columnPhi
                  SELECT CASE(esSpecification(3))
                  CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE)
                    DO rowXiIdx=1,numberOfDependentXi
                      DO columnXiIdx=1,numberOfDependentXi
                        sum=sum+rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                          & geometricInterpPointMetrics%gu(rowXiIdx,columnXiIdx)
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  CASE(EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
                    DO rowXiIdx=1,numberOfDependentXi
                      DO columnXiIdx=1,numberOfDependentXi
                        DO xiIdx=1,numberOfGeometricXi
                          sum=sum+conductivityTensor(rowXiIdx,columnXiIdx)*rowdPhidXi(xiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx                    
                  CASE DEFAULT
                    localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                      & " is not valid for a Helmholtz equation type of a classical field equations set class."
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
            IF(updateMatrix) THEN
              !Loop over element columns
              columnElementDOFIdx=0
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
            ENDIF
            IF(rhsVector%updateVector) rhsVector%elementVector%vector(rowElementDOFIdx)= &
              & rhsVector%elementVector%vector(rowElementDOFIdx)* &
              & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling

    ENDIF !update
       
    EXITS("Helmholtz_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Helmholtz_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Helmholtz equation type of an classical field equations set class.
  SUBROUTINE Helmholtz_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Helmholtz_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
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
        & " is not valid for a Helmholtz equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Helmholtz_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Helmholtz_EquationsSetSolutionMethodSet",err,error)
    EXITS("Helmholtz_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Helmholtz equation type of a classical field equations set class.
  SUBROUTINE Helmholtz_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Helmholtz_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    subtype=specification(3)
    
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a Helmholtz type of a classical field equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,subtype]

    EXITS("Helmholtz_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Helmholtz_EquationsSetSpecificationSet",err,error)
    EXITS("Helmholtz_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Helmholtz equation.
  SUBROUTINE Helmholtz_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricComponentNumber,geometricMeshComponent,geometricScalingType, &
      & numberOfmaterialsComponents,numberOfmaterialsVariables,numberOfDimensions,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Helmholtz_EquationsSetSetup",err,error,*999)
   
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
        & " is not valid for a Helmholtz type of a classical field equations set."
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
        CALL Helmholtz_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)               
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Helmholtz equation."
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
      CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
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
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,"Phi",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1,"del Phi/del n", &
            & err,error,*999)
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
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
          & " is invalid for a standard Helmholtz equation"
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
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      numberOfMaterialsVariables=1
      numberOfMaterialsComponents=1 ! wave number k
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE) &
        & numberOfMaterialsComponents=numberOfMaterialsComponents+NUMBER_OF_VOIGT(numberOfDimensions) !conductivity tensor
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
          CALL Field_NumberOfVariablesSet(equationsMaterials%materialsField,numberOfMaterialsVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          !Default to constant interpolation of material field
          DO componentIdx=1,numberOfMaterialsComponents
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
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
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & 1,1.0_DP,err,error,*999) !k
          IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE) THEN
            DO componentIdx=2,numberOfDimensions+1
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999) !diagonal sigma
            ENDDO !componentIdx
            DO componentIdx=numberOfDimensions+2,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999) !off-diagonal sigma
            ENDDO !componentIdx
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Helmholtz equation."
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
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(equationsSetSetup%analyticFunctionType)
        CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1)
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
          equationsAnalytic%analyticFunctionType=EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1
        CASE DEFAULT
          localError="The specified analytic function type of "// &
            & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
            & " is invalid for a standard Helmholtz equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
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
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
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
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Helmholtz equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("Helmholtz_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Helmholtz_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets up the Helmholtz problem.
  SUBROUTINE Helmholtz_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a Helmholtz equation on.
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
    
    ENTERS("Helmholtz_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE, &
      & PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid for a Helmholtz equation type."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    SELECT CASE(problemSetup%setupType)
    CASE(PROBLEM_SETUP_INITIAL_TYPE)
      !-----------------------------------------------------------------
      ! I n i t i a l   S e t u p
      !-----------------------------------------------------------------
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Do nothing????
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Do nothing???
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      !-----------------------------------------------------------------
      ! C o n t r o l   L o o p s
      !-----------------------------------------------------------------
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
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      !-----------------------------------------------------------------
      ! S o l v e r s
      !-----------------------------------------------------------------
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
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! S o l v e r   E q u a t i o n s
      !-----------------------------------------------------------------
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
          & " is invalid for a standard Helmholtz equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Helmholtz equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Helmholtz_ProblemSetup")
    RETURN
999 ERRORSEXITS("Helmholtz_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Helmholtz equation type.
  SUBROUTINE Helmholtz_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Helmholtz_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE, &
      & PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Helmholtz equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_HELMHOLTZ_EQUATION_TYPE,problemSubtype]

    EXITS("Helmholtz_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Helmholtz_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_ProblemSpecificationSet

  !
  !================================================================================================================================
  !
 
END MODULE HelmholtzEquationsRoutines
