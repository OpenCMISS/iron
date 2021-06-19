!> \file
!> \author Chris Bradley
!> \brief This module contains all interface conditions operators routines.
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
!> Contributor(s): Chris Bradley, Xiani (Nancy) Yan, Thiranja Prasad Babarenda Gamage
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

!>This module contains all interface conditions routines. 
MODULE InterfaceOperatorsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE Constants
  USE DecompositionAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE InterfaceAccessRoutines
  USE InterfaceConditionAccessRoutines
  USE InterfaceEquationsRoutines
  USE InterfaceEquationsAccessRoutines
  USE InterfaceMappingRoutines
  USE InterfaceMappingAccessRoutines
  USE InterfaceMatricesRoutines
  USE InterfaceMatricesAccessRoutines
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE MeshAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

  !Module types

  !Module variables

  !Interfaces

  PUBLIC FieldContinuity_FiniteElementCalculate
  
  PUBLIC FrictionlessContact_FiniteElementCalculate

  PUBLIC SolidFluidOperator_FiniteElementCalculate

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for field continuity operator 
  SUBROUTINE FieldContinuity_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,connectedLineFace,coupledElementNumber, &
      & coupledMeshIdx,coupledDependentScalingType,coupledDependentVariableType,dataPointIdx,coupledDecompositionFaceNumber, &
      & coupledDecompositionLineNumber,derivative,derivativeIdx,gaussPointIdx,interfaceConditionMethod, &
      & interfaceDependentScalingType,interfaceDerivative,interfaceDerivativeIdx,interfaceElementNumber,interfaceMatrixIdx, &
      & interfaceNode,interfaceNodeIdx,lagrangeVariableType,integrationType,localElementNode,localElementNodeIdx, &
      & localElementNumber,localFaceNodeIdx,localLineNodeIdx,localNodeIdx,matrixElementIdx,numberOfCoupledElements, &
      & numberOfCoupledMeshes,numberOfCoupledDerivatives,numberOfCoupledElementParameters,numberOfCoupledGeometricComponents, &
      & numberOfCoupledNodes,numberOfCoupledVariableComponents,numberOfCoupledXi,numberOfElementDataPoints,numberOfGauss, &
      & numberOfInterfaceGeometricXi,numberOfInterfaceNodeDerivatives,numberOfInterfaceDependentElementParameters, &
      & numberOfInterfaceDependentNodes,numberOfInterfaceDependentXi,numberOfInterfaceMatrices,numberOfInterfaceMeshXi, &
      & numberOfInterfaceVariableComponents,numberOfLagrangeComponents,numberOfMatrixCoupledElements,numberOfFaceNodeDerivatives, &
      & numberOfLineNodeDerivatives,numberOfNodesInFace,numberOfNodesInLine,numberOfRowElementParameters,numberOfXi, &
      & rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx  
    REAL(DP) :: columnBasisFunction,gaussWeight,jacobian,jacobianGaussWeight,matrixCoefficient,rowBasisFunction,xi(3)
    LOGICAL :: continueSearch,found,updateMatrix,updateMatrices
    TYPE(BasisType), POINTER :: coupledBasis,coupledDependentBasis,coupledFaceBasis,coupledLineBasis,interfaceDependentBasis, &
      & interfaceGeometricBasis,interfacePenaltyBasis,interfaceConnectivityBasis,rowBasis
    TYPE(DecompositionType), POINTER :: coupledDependentDecomposition,interfaceDependentDecomposition, &
      & interfaceGeometricDecomposition,lagrangeDecomposition,penaltyDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionElementDataPointsType), POINTER :: elementDataPoints
    TYPE(DecompositionElementsType), POINTER :: coupledDependentDecompositionElements,interfaceDependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: coupledDependentDecompositionTopology,interfaceDependentDecompositionTopology, &
      & lagrangeDecompositionTopology
    TYPE(DomainType), POINTER :: coupledDependentDomain,interfaceDependentDomain,interfaceGeometricDomain,penaltyDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: coupledDependentDomainElements,interfaceDependentDomainElements, &
      & interfaceGeometricDomainElements,penaltyDomainElements,rowDomainElements
    TYPE(DomainFaceType), POINTER :: coupledDomainFace
    TYPE(DomainFacesType), POINTER :: coupledDependentDomainFaces
    TYPE(DomainLineType), POINTER :: coupledDomainLine
    TYPE(DomainLinesType), POINTER :: coupledDependentDomainLines
    TYPE(DomainTopologyType), POINTER :: coupledDependentDomainTopology,interfaceDependentDomainTopology, &
      & interfaceGeometricDomainTopology,penaltyDomainTopology,rowDomainTopology
    TYPE(FieldType), POINTER :: coupledDependentField,coupledGeometricField,interfaceDependentField,interfaceGeometricField, &
      & interfacePenaltyField,lagrangeField
    TYPE(FieldInterpolationParametersType), POINTER :: coupledDependentInterpParameters,interfaceGeometricInterpParameters, &
      & lagrangeInterpParameters,penaltyInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: interfaceGeometricInterpPoint,penaltyInterpolatedPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interfaceGeometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: coupledGeometricVariable,coupledDependentVariable,lagrangeVariable
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceCoupledElementsType), POINTER :: coupledElements
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations 
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: interfaceInterpolation,coupledInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: coupledDependentInterpSet,interfaceGeometricInterpSet, &
      & lagrangeInterpSet,penaltyInterpSet
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity     
    TYPE(InterfacePointConnectivityType), POINTER :: pointConnectivity
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh     
    TYPE(QuadratureSchemeType), POINTER :: interfaceQuadratureScheme,rowQuadratureScheme
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FieldContinuity_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)
 
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
    updateMatrices=.FALSE.
    DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
      updateMatrices=(updateMatrices.OR.updateMatrix)
    ENDDO !interfaceMatrixIdx
!!TODO: Consider interface RHS

    IF(updateMatrices) THEN
    
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
      NULLIFY(INTERFACE)
      CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
      NULLIFY(interfaceEquationsInterpolation)
      CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
      
      CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        CALL InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*999)
        SELECT CASE(integrationType)
        CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
          NULLIFY(meshConnectivity)
          CALL Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*999)
          NULLIFY(interfaceConnectivityBasis)
          CALL InterfaceMeshConnectivity_BasisGet(meshConnectivity,interfaceConnectivityBasis,err,error,*999)
          !Pointers to interface variables (columns of interface element matrix)
          NULLIFY(interfaceInterpolation)
          CALL InterfaceEquationsInterpolation_InterfaceInterpGet(interfaceEquationsInterpolation,interfaceInterpolation, &
            & err,error,*999)
          NULLIFY(interfaceGeometricField)
          CALL InterfaceDomainInterpolation_GeometricFieldGet(interfaceInterpolation,interfaceGeometricField,err,error,*999)
          NULLIFY(interfaceDependentField)
          CALL InterfaceDomainInterpolation_DependentFieldGet(interfaceInterpolation,interfaceDependentField,err,error,*999)
          NULLIFY(interfaceGeometricDecomposition)
          CALL Field_DecompositionGet(interfaceGeometricField,interfaceGeometricDecomposition,err,error,*999)
          NULLIFY(interfaceGeometricDomain)
          CALL Decomposition_DomainGet(interfaceGeometricDecomposition,0,interfaceGeometricDomain,err,error,*999)
          NULLIFY(interfaceGeometricDomainTopology)
          CALL Domain_DomainTopologyGet(interfaceGeometricDomain,interfaceGeometricDomainTopology,err,error,*999)
          NULLIFY(interfaceGeometricDomainElements)
          CALL DomainTopology_DomainElementsGet(interfaceGeometricDomainTopology,interfaceGeometricDomainElements,err,error,*999)
          NULLIFY(interfaceGeometricBasis)
          CALL DomainElements_ElementBasisGet(interfaceGeometricDomainElements,elementNumber,interfaceGeometricBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(interfaceGeometricBasis,numberOfInterfaceGeometricXi,err,error,*999)
          NULLIFY(interfaceDependentDecomposition)
          CALL Field_DecompositionGet(interfaceDependentField,interfaceDependentDecomposition,err,error,*999)
          NULLIFY(interfaceDependentDomain)
          CALL Decomposition_DomainGet(interfaceDependentDecomposition,0,interfaceDependentDomain,err,error,*999)
          NULLIFY(interfaceDependentDomainTopology)
          CALL Domain_DomainTopologyGet(interfaceDependentDomain,interfaceDependentDomainTopology,err,error,*999)
          NULLIFY(interfaceDependentDomainElements)
          CALL DomainTopology_DomainElementsGet(interfaceDependentDomainTopology,interfaceDependentDomainElements,err,error,*999)
          NULLIFY(interfaceDependentBasis)
          CALL DomainElements_ElementBasisGet(interfaceDependentDomainElements,elementNumber,interfaceDependentBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(interfaceDependentBasis,numberOfInterfaceDependentXi,err,error,*999)
          CALL Basis_NumberOfLocalNodesGet(interfaceDependentBasis,numberOfInterfaceDependentNodes,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(interfaceDependentBasis,numberOfInterfaceDependentElementParameters, &
            & err,error,*999)
          SELECT CASE(interfaceConditionMethod)
          CASE(INTERFACE_CONDITION_PENALTY_METHOD)
            CALL InterfaceDomainInterpolation_PenaltyFieldGet(interfaceInterpolation,interfacePenaltyField,err,error,*999)
            NULLIFY(penaltyDecomposition)
            CALL Field_DecompositionGet(interfacePenaltyField,penaltyDecomposition,err,error,*999)
            NULLIFY(penaltyDomain)
            CALL Decomposition_DomainGet(penaltyDecomposition,0,penaltyDomain,err,error,*999)
            NULLIFY(penaltyDomainTopology)
            CALL Domain_DomainTopologyGet(penaltyDomain,penaltyDomainTopology,err,error,*999)
            NULLIFY(penaltyDomainElements)
            CALL DomainTopology_DomainElementsGet(penaltyDomainTopology,penaltyDomainElements,err,error,*999)
            NULLIFY(interfacePenaltyBasis)
            CALL DomainElements_ElementBasisGet(penaltyDomainElements,elementNumber,interfacePenaltyBasis,err,error,*999)
            NULLIFY(penaltyInterpSet)
            CALL InterfaceDomainInterpolation_PenaltyInterpSetGet(interfaceInterpolation,1,penaltyInterpSet,err,error,*999)
            NULLIFY(penaltyInterpParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(penaltyInterpSet,FIELD_U_VARIABLE_TYPE, &
              & penaltyInterpParameters,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,penaltyInterpParameters, &
              & err,error,*999)
          END SELECT
          !Integrate using the interface quadrature scheme
          NULLIFY(interfaceQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(interfaceGeometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,interfaceQuadratureScheme, &
            & err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(interfaceQuadratureScheme,numberOfGauss,err,error,*999)
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          CALL FieldVariable_VariableTypeGet(lagrangeVariable,lagrangeVariableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(lagrangeVariable,numberOfLagrangeComponents,err,error,*999)
          NULLIFY(lagrangeInterpSet)
          CALL InterfaceDomainInterpolation_DependentInterpSetGet(interfaceInterpolation,1,lagrangeInterpSet, err,error,*999)
          NULLIFY(lagrangeInterpParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(lagrangeInterpSet,lagrangeVariableType, &
            & lagrangeInterpParameters,err,error,*999)
          !Get element interpolation parameters from the first geo interp set (to get Jacobian for interface surface integral)
          NULLIFY(interfaceGeometricInterpSet)
          CALL InterfaceDomainInterpolation_GeometricInterpSetGet(interfaceInterpolation,1,interfaceGeometricInterpSet, &
            & err,error,*999)
          NULLIFY(interfaceGeometricInterpParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpParameters,err,error,*999)
          NULLIFY(interfaceGeometricInterpPoint)
          CALL InterfaceInterpolationSet_InterpolatedPointGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpPoint,err,error,*999)
          NULLIFY(interfaceGeometricInterpPointMetrics)
          CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpPointMetrics,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,interfaceGeometricInterpParameters, &
            & err,error,*999)
          CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
          DO coupledMeshIdx=1,numberOfInterfaceMatrices
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
            IF(updateMatrix) THEN
              !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
              NULLIFY(elementConnectivity)
              CALL InterfaceMeshConnectivity_CoupledElementGet(meshConnectivity,elementNumber,coupledMeshIdx,elementConnectivity, &
                & err,error,*999)
              numberOfXi=SIZE(elementConnectivity%xi,1)
              NULLIFY(coupledInterpolation)
              CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
                & coupledInterpolation,err,error,*999)
              NULLIFY(coupledDependentField)
              CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
              NULLIFY(coupledDependentDecomposition)
              CALL Field_DecompositionGet(coupledDependentField,coupledDependentDecomposition,err,error,*999)
              NULLIFY(coupledDependentDecompositionTopology)
              CALL Decomposition_DecompositionTopologyGet(coupledDependentDecomposition,coupledDependentDecompositionTopology, &
                & err,error,*999)
              NULLIFY(coupledDependentDecompositionElements)
              CALL DecompositionTopology_DecompositionElementsGet(coupledDependentDecompositionTopology, &
                & coupledDependentDecompositionElements,err,error,*999)
              CALL InterfaceMeshConnectivity_CoupledElementNumberGet(meshConnectivity,elementNumber,coupledMeshIdx, &
                & coupledElementNumber,err,error,*999)
              NULLIFY(coupledDependentVariable)
              CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,coupledDependentVariable,err,error,*999)
              CALL FieldVariable_VariableTypeGet(coupledDependentVariable,coupledDependentVariableType,err,error,*999)
              NULLIFY(coupledDependentInterpSet)
              CALL InterfaceDomainInterpolation_DependentInterpSetGet(coupledInterpolation,1,coupledDependentInterpSet, &
                & err,error,*999)
              NULLIFY(coupledDependentInterpParameters)
              CALL InterfaceInterpolationSet_InterpolationParametersGet(coupledDependentInterpSet,coupledDependentVariableType, &
                & coupledDependentInterpParameters,err,error,*999)
              
              !Loop over gauss points
              DO gaussPointIdx=1,numberOfGauss
                CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & interfaceGeometricInterpPoint,err,error,*999)
                CALL Field_InterpolatedPointMetricsCalculate(numberOfInterfaceGeometricXi,interfaceGeometricInterpPointMetrics, &
                  & err,error,*999)
                CALL FieldInterpolatedPointMetrics_JacobianGet(interfaceGeometricInterpPointMetrics,jacobian,err,error,*999)
                CALL BasisQuadratureScheme_GaussWeightGet(interfaceQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
                jacobianGaussWeight=jacobian*gaussWeight
                IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD.AND.coupledMeshIdx==numberOfInterfaceMatrices) THEN
                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                    & penaltyInterpolatedPoint,err,error,*999)
                  !Loop over the Lagrange variable matrix rows (and columns?)
                  rowElementDOFIdx=0
                  DO rowComponentIdx=1,numberOfLagrangeComponents
                    NULLIFY(rowDomain)
                    CALL FieldVariable_ComponentDomainGet(lagrangeVariable,rowComponentIdx,rowDomain,err,error,*999)
                    NULLIFY(rowDomainTopology)
                    CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
                    NULLIFY(rowDomainElements)
                    CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
                    NULLIFY(rowBasis)
                    CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis,err,error,*999)
                    NULLIFY(rowQuadratureScheme)
                    CALL Basis_QuadratureSchemeGet(rowBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowQuadratureScheme,err,error,*999)
                    CALL Basis_NumberOfElementParametersGet(rowBasis,numberOfRowElementParameters,err,error,*999)
                    DO rowElementParameterIdx=1,numberOfRowElementParameters
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                        & gaussPointIdx,rowBasisFunction,err,error,*999)
                      rowElementDOFIdx=rowElementDOFIdx+1
                      interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,rowElementDOFIdx)= &
                        & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,rowElementDOFIdx)- &
                        & (1.0_DP/penaltyInterpolatedPoint%values(1,1))*rowBasisFunction**2.0_DP*jacobianGaussWeight
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ELSE
                  !\todo defaults to first mesh component, generalise
                  CALL InterfaceOperators_InterfaceToCoupledMeshGaussTransform(elementConnectivity,interfaceConnectivityBasis, &
                    & gaussPointIdx,xi,err,error,*999)
                  !XI=interfaceCondition%interface%pointsConnectivity%pointsConnectivity(gaussPointIdx,coupledMeshIdx)%xi
                  !Loop over number of Lagrange variable components as not all components in the dependent field variable may
                  !be coupled
                  !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                  !component numbers. Generalise ordering
                  DO rowComponentIdx=1,numberOfLagrangeComponents
                    NULLIFY(coupledDependentDomain)
                    CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                      & err,error,*999)
                    NULLIFY(coupledDependentDomainTopology)
                    CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                    NULLIFY(coupledDependentDomainElements)
                    CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                      & err,error,*999)
                    NULLIFY(coupledBasis)
                    CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,coupledElementNumber,coupledBasis, &
                      & err,error,*999)
                    CALL Basis_NumberOfXiGet(coupledBasis,numberOfCoupledXi,err,error,*999)
                    CALL Basis_NumberOfLocalNodesGet(coupledBasis,numberOfCoupledNodes,err,error,*999)
                    CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                    
                    SELECT CASE(numberOfInterfaceDependentXi)
                    CASE(1) !1D interface (line)
                      
                      NULLIFY(coupledDependentDomainLines)
                      CALL DomainTopology_DomainLinesGet(coupledDependentDomainTopology,coupledDependentDomainLines,err,error,*999)
                      CALL InterfaceElementConnectivity_CoupledLineFaceNumberGet(elementConnectivity,connectedLineFace, &
                        & err,error,*999)
                      CALL DecompositionElements_ElementLineNumberGet(coupledDependentDecompositionElements,connectedLineFace, &
                        & coupledElementNumber,coupledDecompositionLineNumber,err,error,*999)
                      NULLIFY(coupledDomainLine)
                      CALL DomainLines_LineGet(coupledDependentDomainLines,coupledDecompositionLineNumber,coupledDomainLine, &
                        & err,error,*999)
                      NULLIFY(coupledLineBasis)
                      CALL DomainLines_LineBasisGet(coupledDependentDomainLines,coupledDecompositionLineNumber,coupledLineBasis, &
                        & err,error,*999)
                      CALL Basis_LineNumberOfNodesGet(coupledBasis,connectedLineFace,numberOfNodesInLine,err,error,*999)
                      DO localLineNodeIdx=1,numberOfNodesInLine
                        CALL Basis_LineNodeNumberGet(coupledBasis,localLineNodeIdx,connectedLineFace,localElementNode, &
                          & err,error,*999)
                        CALL Basis_NodeNumberOfDerivativesGet(coupledLineBasis,localLineNodeIdx,numberOfLineNodeDerivatives, &
                          & err,error,*999)
                        DO derivativeIdx=1,numberOfLineNodeDerivatives
                          CALL DomainLine_DerivativeGlobalIndexGet(coupledDomainLine,derivativeIdx,localLineNodeIdx,derivative, &
                            & err,error,*999)
                          CALL Basis_ElementParameterGet(coupledBasis,derivative,localElementNode,rowElementParameterIdx, &
                            & err,error,*999)
                          CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:2),rowBasisFunction, &
                            & err,error,*999)
                          rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                          DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                            CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx, &
                              & numberOfInterfaceNodeDerivatives,err,error,*999)
                            DO interfaceDerivativeIdx=1,numberOfInterfaceNodeDerivatives
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                                & columnElementParameterIdx,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme, &
                                & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                              columnElementDOFIdx=columnElementParameterIdx+ &
                                & numberOfInterfaceDependentElementParameters*(rowComponentIdx-1)
                              interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                                & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                                & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                            ENDDO !interfaceDerivativeIdx
                          ENDDO !interfaceNodeIdx
                        ENDDO !derivativeIdx
                      ENDDO !localLineNodeIdx
                      
                    CASE(2) !2D interface (face)

                      NULLIFY(coupledDependentDomainFaces)
                      CALL DomainTopology_DomainFacesGet(coupledDependentDomainTopology,coupledDependentDomainFaces,err,error,*999)
                      
                      SELECT CASE(numberOfCoupledXi)                        
                      CASE(2) !Coupled Mesh has 2 xi directions
                        DO localElementNodeIdx=1,numberOfCoupledNodes
                          CALL Basis_NodeNumberOfDerivativesGet(coupledBasis,localElementNodeIdx,numberOfCoupledDerivatives, &
                            & err,error,*999)
                          DO derivativeIdx=1,numberOfCoupledDerivatives
                            CALL Basis_ElementParameterGet(coupledBasis,derivativeIdx,localElementNodeIdx, &
                              & rowElementParameterIdx,err,error,*999)
                            CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:numberOfCoupledXi), &
                              & rowBasisFunction,err,error,*999)
                            rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                            DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                              CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx, &
                                & numberOfInterfaceNodeDerivatives,err,error,*999)
                              DO interfaceDerivativeIdx=1,numberOfInterfaceNodeDerivatives
                                !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                                CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                                  & columnElementParameterIdx,err,error,*999)
                                CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme, &
                                  & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                                columnElementDOFIdx=columnElementParameterIdx+ &
                                  & numberOfInterfaceDependentElementParameters*(rowComponentIdx-1)
                                interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                                  & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                                  & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                              ENDDO !interfaceDerivativeIdx
                            ENDDO !interfaceNodeIdx
                          ENDDO !derivativeIdx
                        ENDDO !localElementNodeIdx
                        
                      CASE(3) !Coupled Mesh has 3 xi directions

                        CALL InterfaceElementConnectivity_CoupledLineFaceNumberGet(elementConnectivity,connectedLineFace, &
                          & err,error,*999)
                        CALL DecompositionElements_ElementFaceNumberGet(coupledDependentDecompositionElements,connectedLineFace, &
                          & coupledElementNumber,coupledDecompositionFaceNumber,err,error,*999)
                        NULLIFY(coupledDomainFace)
                        CALL DomainFaces_FaceGet(coupledDependentDomainFaces,coupledDecompositionFaceNumber,coupledDomainFace, &
                          & err,error,*999)
                        NULLIFY(coupledFaceBasis)
                        CALL DomainFace_BasisGet(coupledDomainFace,coupledFaceBasis,err,error,*999)
                        CALL Basis_FaceNumberOfNodesGet(coupledBasis,connectedLineFace,numberOfNodesInFace,err,error,*999)
                        DO localFaceNodeIdx=1,numberOfNodesInFace
                          CALL Basis_FaceNodeNumberGet(coupledBasis,localFaceNodeIdx,connectedLineFace,localElementNode, &
                            & err,error,*999)
                          CALL Basis_NodeNumberOfDerivativesGet(coupledFaceBasis,localFaceNodeIdx,numberOfFaceNodeDerivatives, &
                            & err,error,*999)
                          DO derivativeIdx=1,numberOfFaceNodeDerivatives
                            CALL Basis_FaceNodeDerivativeNumberGet(coupledBasis,derivativeIdx,localFaceNodeIdx, &
                              & connectedLineFace,derivative,err,error,*999)
                            CALL Basis_ElementParameterGet(coupledBasis,derivative,localElementNode,rowElementParameterIdx, &
                              & err,error,*999)
                            CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:numberOfCoupledXi), &
                              & rowBasisFunction,err,error,*999)
                            rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                            DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                              CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx, &
                                & numberOfInterfaceNodeDerivatives,err,error,*999)
                              DO interfaceDerivativeIdx=1,numberOfInterfaceNodeDerivatives
                                !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                                CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                                  & columnElementParameterIdx,err,error,*999)
                                CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme, &
                                  & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                                columnElementDOFIdx=columnElementParameterIdx+ &
                                  & numberOfInterfaceDependentElementParameters*(rowComponentIdx-1)
                                interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                                  & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                                  & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                              ENDDO !interfaceDerivativeIdx
                            ENDDO !interfaceNodeIdx
                          ENDDO !derivativeIdx
                        ENDDO !FaceNodeIdx
                        
                      END SELECT !coupledBasis%numberOfXi
                      
                    END SELECT !interfaceDependentBasis%numberOfXi

                  ENDDO !rowComponentIdx
                ENDIF
              ENDDO !gaussPointIdx

              !Scale factor adjustment
              !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix
              !contribution to the residual for non-linear problems
              !\todo update looping of variables/components for non-zero matrix elements as done above
              CALL Field_ScalingTypeGet(interfaceDependentField,interfaceDependentScalingType,err,error,*999)
              CALL Field_ScalingTypeGet(coupledDependentField,coupledDependentScalingType,err,error,*999)
              IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD.AND.coupledMeshIdx==numberOfInterfaceMatrices) THEN
                !Scale factor adjustment for the Lagrange Variable (columns)
                IF(interfaceDependentScalingType/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpParameters, &
                    & err,error,*999)
                  rowElementDOFIdx=0
                  !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                  !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                  !component numbers. Generalise ordering
                  DO rowComponentIdx=1,numberOfLagrangeComponents
                    !Loop over element Lagrange variable rows
                    DO rowElementParameterIdx=1,numberOfInterfaceDependentElementParameters
                      rowElementDOFIdx=rowElementDOFIdx+1
                      interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,rowElementDOFIdx)= &
                        & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,rowElementDOFIdx)* &
                        & lagrangeInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)**2
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ENDIF
              ELSE
                !Scale factor adjustment for the Lagrange Variable (columns)
                IF(interfaceDependentScalingType/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpParameters,err,error,*999)
                  rowElementDOFIdx=0
                  !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                  !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                  !component numbers. Generalise ordering
                  DO rowComponentIdx=1,numberOfLagrangeComponents
                    NULLIFY(coupledDependentDomain)
                    CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                      & err,error,*999)
                    NULLIFY(coupledDependentDomainTopology)
                    CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                    NULLIFY(coupledDependentDomainElements)
                    CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                      & err,error,*999)
                    NULLIFY(coupledBasis)
                    CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,coupledElementNumber,coupledBasis, &
                      & err,error,*999)
                    !Loop over element rows
                    CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                    DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                      rowElementDOFIdx=rowElementDOFIdx+1
                      columnElementDOFIdx=0
                      !Loop over element columns
                      DO columnComponentIdx=1,numberOfLagrangeComponents
                        DO columnElementParameterIdx=1,numberOfinterfaceDependentElementParameters
                          columnElementDOFIdx=columnElementDOFIdx+1
                          interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                            & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                            & lagrangeInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                        ENDDO !columnElementParameterIdx
                      ENDDO !columnComponentIdx
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ENDIF
                !Scale factor adjustment for the row dependent variable
                IF(coupledDependentScalingType/=FIELD_NO_SCALING) THEN
                  CALL Field_InterpolationParametersScaleFactorsElementGet(coupledElementNumber,coupledDependentInterpParameters, &
                    & err,error,*999)
                  CALL FieldVariable_NumberOfComponentsGet(coupledDependentVariable,numberOfCoupledVariableComponents, &
                    & err,error,*999)
                  !Loop over element rows
                  rowElementDOFIdx=0                 
                  DO rowComponentIdx=1,numberOfCoupledVariableComponents
                    NULLIFY(coupledDependentDomain)
                    CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                      & err,error,*999)
                    NULLIFY(coupledDependentDomainTopology)
                    CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                    NULLIFY(coupledDependentDomainElements)
                    CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                      & err,error,*999)
                    NULLIFY(coupledBasis)
                    CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,coupledElementNumber,coupledBasis, &
                      & err,error,*999)
                    CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                    !Loop over element rows
                    DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                      rowElementDOFIdx=rowElementDOFIdx+1
                      columnElementDOFIdx=0
                      !Loop over element columns
                      DO columnComponentIdx=1,numberOfLagrangeComponents
                        DO columnElementParameterIdx=1,numberOfInterfaceDependentElementParameters
                          columnElementDOFIdx=columnElementDOFIdx+1
                          interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                            & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                            & coupledDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
                        ENDDO !columnElementParameterIdx
                      ENDDO !columnComponentIdx
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ENDIF
              ENDIF
            ENDIF
          ENDDO ! coupledMeshIdx
          
        CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)

          interfaceElementNumber = elementNumber !todo simplify
          NULLIFY(pointsConnectivity)
          CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)  
          NULLIFY(interfaceMesh)
          CALL InterfacePointsConnectivity_InterfaceMeshGet(pointsConnectivity,interfaceMesh,err,error,*999)
          CALL Mesh_NumberOfDimensionsGet(interfaceMesh,numberOfInterfaceMeshXi,err,error,*999)
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          NULLIFY(lagrangeField)
          CALL FieldVariable_FieldGet(lagrangeVariable,lagrangeField,err,error,*999)
          NULLIFY(lagrangeDecomposition)
          CALL Field_DecompositionGet(lagrangeField,lagrangeDecomposition,err,error,*999)
          NULLIFY(lagrangeDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(lagrangeDecomposition,lagrangeDecompositionTopology,err,error,*999)
          NULLIFY(dataPoints)
          CALL DecompositionTopology_DecompositionDataPointsGet(lagrangeDecompositionTopology,dataPoints,err,error,*999)
          NULLIFY(elementDataPoints)
          CALL DecompositionDataPoints_ElementDataPointsGet(dataPoints,interfaceElementNumber,elementDataPoints,err,error,*999)
          CALL DecompositionElementDataPoints_NumberOfDataPointsGet(elementDataPoints,numberOfElementDataPoints,err,error,*999)
          !Calculate row basis function, update interface matrices with row basis function, and update scale factors
          CALL Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*999)
          DO coupledMeshIdx=1,numberOfCoupledMeshes
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
            IF(updateMatrix) THEN
              NULLIFY(coupledElements)
              CALL InterfacePointsConnectivity_CoupledElementsGet(pointsConnectivity,interfaceElementNumber,coupledMeshIdx, &
                & coupledElements,err,error,*999)
              CALL InterfaceCoupledElements_NumberOfCoupledElementsGet(coupledElements,numberOfMatrixCoupledElements, &
                & err,error,*999)
              NULLIFY(coupledMesh)
              CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
              CALL Mesh_NumberOfDimensionsGet(coupledMesh,numberOfCoupledXi,err,error,*999)              
              NULLIFY(coupledInterpolation)
              CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
                & coupledInterpolation,err,error,*999)
              NULLIFY(coupledGeometricField)
              CALL InterfaceDomainInterpolation_GeometricFieldGet(coupledInterpolation,coupledGeometricField,err,error,*999)
              NULLIFY(coupledGeometricVariable)
              CALL Field_VariableGet(coupledGeometricField,FIELD_U_VARIABLE_TYPE,coupledGeometricVariable,err,error,*999)
              NULLIFY(coupledDependentField)
              CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
              NULLIFY(coupledDependentVariable)
              CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,coupledDependentVariable,err,error,*999)
              CALL FieldVariable_VariableTypeGet(coupledDependentVariable,coupledDependentVariableType,err,error,*999)
              NULLIFY(coupledDependentInterpSet)
              CALL InterfaceDomainInterpolation_DependentInterpSetGet(coupledInterpolation,1,coupledDependentInterpSet, &
                & err,error,*999)
              NULLIFY(coupledDependentInterpParameters)
              CALL InterfaceInterpolationSet_InterpolationParametersGet(coupledDependentInterpSet,coupledDependentVariableType, &
                & coupledDependentInterpParameters,err,error,*999)
              CALL FieldVariable_NumberOfComponentsGet(coupledGeometricVariable,numberOfCoupledGeometricComponents,err,error,*999)
              DO dataPointIdx=1,numberOfElementDataPoints
                NULLIFY(pointConnectivity)
                CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,dataPointIdx,coupledMeshIdx, &
                  & pointConnectivity,err,error,*999)
                CALL InterfacePointConnectivity_CoupledElementNumberGet(pointConnectivity,localElementNumber,err,error,*999)
                !Calculate the element index (non-conforming element) for this interface matrix
                CALL InterfaceCoupledElements_NumberOfCoupledElementsGet(coupledElements,numberOfCoupledElements,err,error,*999)
                matrixElementIdx=1
                found=.FALSE.
                continueSearch=.TRUE.
                DO WHILE(continueSearch)
                  CALL InterfaceCoupledElements_CoupledElementNumberGet(coupledElements,matrixElementIdx,coupledElementNumber, &
                    & err,error,*999)
                  IF(localElementNumber==coupledElementNumber) THEN
                    found=.TRUE.
                    continueSearch=.FALSE.
                  ELSE
                    IF(matrixElementIdx==numberOfCoupledElements) THEN
                      continueSearch=.FALSE.
                    ELSE
                      matrixElementIdx=matrixElementIdx+1
                    ENDIF
                  ENDIF
                ENDDO
                CALL InterfacePointConnectivity_XiGet(pointConnectivity,numberOfXi,xi,err,error,*999)
                DO rowComponentIdx=1,numberOfCoupledGeometricComponents
                  NULLIFY(coupledDependentDomain)
                  CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                    & err,error,*999)
                  NULLIFY(coupledDependentDomainTopology)
                  CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                  NULLIFY(coupledDependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                    & err,error,*999)
                  NULLIFY(coupledDependentBasis)
                  CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,localElementNumber,coupledDependentBasis, &
                    & err,error,*999)
                  CALL Basis_NumberOfElementParametersGet(coupledDependentBasis,numberOfCoupledElementParameters, &
                    & err,error,*999)
                  DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                    CALL Basis_EvaluateXi(coupledDependentBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:numberOfXi), &
                      & rowBasisFunction,err,error,*999)
                    rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                    columnElementDOFIdx=dataPointIdx+numberOfElementDataPoints*(rowComponentIdx-1)
                    !Update interface element matrix with contact point contribution
                    interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)=rowBasisFunction
                  ENDDO !rowElementParameterIdx
                ENDDO !rowComponentIdx
              ENDDO !dataPointIdx
              
              !scale factor update
              CALL Field_ScalingTypeGet(coupledDependentField,coupledDependentScalingType,err,error,*999)
              IF(coupledDependentScalingType/=FIELD_NO_SCALING) THEN
                DO dataPointIdx=1,numberOfElementDataPoints
                  NULLIFY(pointConnectivity)
                  CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,dataPointIdx,coupledMeshIdx, &
                    & pointConnectivity,err,error,*999)
                  CALL InterfacePointConnectivity_CoupledElementNumberGet(pointConnectivity,localElementNumber,err,error,*999)
                  CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,coupledDependentInterpParameters, &
                    & err,error,*999)
                  
                  !Calculate the element index (non-conforming element) for this interface matrix
                  CALL InterfaceCoupledElements_NumberOfCoupledElementsGet(coupledElements,numberOfCoupledElements,err,error,*999)
                  matrixElementIdx=1
                  found=.FALSE.
                  continueSearch=.TRUE.
                  DO WHILE(continueSearch)
                    CALL InterfaceCoupledElements_CoupledElementNumberGet(coupledElements,matrixElementIdx,coupledElementNumber, &
                      & err,error,*999)
                    IF(localElementNumber==coupledElementNumber) THEN
                      found=.TRUE.
                      continueSearch=.FALSE.
                    ELSE
                      IF(matrixElementIdx==numberOfCoupledElements) THEN
                        continueSearch=.FALSE.
                      ELSE
                        matrixElementIdx=matrixElementIdx+1
                      ENDIF
                    ENDIF
                  ENDDO
                  DO rowComponentIdx=1,numberOfCoupledGeometricComponents
                    NULLIFY(coupledDependentDomain)
                    CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                      & err,error,*999)
                    NULLIFY(coupledDependentDomainTopology)
                    CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                    NULLIFY(coupledDependentDomainElements)
                    CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                      & err,error,*999)
                    NULLIFY(coupledDependentBasis)
                    CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,localElementNumber,coupledDependentBasis, &
                      & err,error,*999)
                    CALL Basis_NumberOfElementParametersGet(coupledDependentBasis,numberOfCoupledElementParameters, &
                      & err,error,*999)
                    DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                      rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                      columnElementDOFIdx=dataPointIdx+numberOfElementDataPoints*(rowComponentIdx-1)
                      interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                        & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                        & coupledDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
                    ENDDO !rowElementParameterIdx
                  ENDDO !rowComponentIdx
                ENDDO !dataPointIdx
              ENDIF !.NOT. FIELD_NO_SCALING
              
            ENDIF !updateMatrix
          ENDDO !coupledMeshIdx
          
        CASE DEFAULT
          localError="Interface condition integration type "//TRIM(NumberToVString(integrationType,"*",err,error))// &
            & " is not valid."
          CALL FlagError(localError,err,error,*999)
        END SELECT !interfaceCondition%integrationType
        
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Interface condition method "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
          & " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
    ENDIF !update

    EXITS("FieldContinuity_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FieldContinuity_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FieldContinuity_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries for the given element number for frictionless contact operator
  SUBROUTINE FrictionlessContact_FiniteElementCalculate(interfaceCondition,interfaceElementNumber,err,error,*)

    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: interfaceElementNumber !<The interface element number to calcualte the interface element matrix for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnElementDOFIdx,componentIdx,coupledElementNumber,coupledMeshIdx,coupledDependentVariableType, &
      & dataPointIdx,elementLineFaceNumber,globalDataPointNumber,interfaceConditionMethod,interfaceMatrixIdx,integrationType, &
      & interpolationType,localDOF,localElementNumber,localFaceLineNumber,matrixElementIdx,numberOfCoupledElements, &
      & numberOfCoupledElementParameters,numberOfCoupledGeometricComponents,numberOfCoupledMeshes,numberOfCoupledXi, &
      & numberOfElementDataPoints,numberOfInterfaceMatrices,numberOfInterfaceXi,numberOfMatrixCoupledElements, &
      & numberOfPenaltyComponents,numberOfReducedXi,numberOfXi,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,xiIdx
    REAL(DP) :: contactStiffness,normalPoint(3),positionPoint(3),reducedXi(3),rowBasisFunction,tangentsPoint(3,3),xi(3)
    REAL(DP), ALLOCATABLE :: gaps(:),gapsComponents(:,:),normals(:,:)
    LOGICAL :: continueSearch,firstAssembly,found,reverseNormal,updateMatrix,updateMatrices
    LOGICAL, ALLOCATABLE :: orthogonallyProjected(:)       
    TYPE(BasisType), POINTER :: coupledDependentBasis
    TYPE(DecompositionType), POINTER :: coupledDependentDecomposition,lagrangeDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionElementDataPointsType), POINTER :: elementDataPoints
    TYPE(DecompositionElementType), POINTER :: coupledDependentDecompositionElement
    TYPE(DecompositionElementsType), POINTER :: coupledDependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: coupledDependentDecompositionTopology,lagrangeDecompositionTopology
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData 
    TYPE(DomainType), POINTER :: coupledDependentDomain
    TYPE(DomainElementsType), POINTER :: coupledDependentDomainElements
    TYPE(DomainTopologyType), POINTER :: coupledDependentDomainTopology
    TYPE(FieldType), POINTER :: coupledDependentField,coupledGeometricField,lagrangeField,penaltyField
    TYPE(FieldInterpolationParametersType), POINTER :: coupledDependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: coupledDependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: coupledDependentInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: coupledGeometricVariable,coupledDependentVariable,lagrangeVariable,penaltyVariable
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceCoupledElementsType), POINTER :: coupledElements
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations 
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: coupledInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: coupledDependentInterpSet
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix,penaltyMatrix
    TYPE(InterfacePointConnectivityType), POINTER :: pointConnectivity 
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity 
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh     
    TYPE(VARYING_STRING) :: localError

    ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
    updateMatrices=.FALSE.
    DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
      updateMatrices=(updateMatrices.OR.updateMatrix)
    ENDDO !interfaceMatrixIdx

    IF(updateMatrices) THEN
    
      NULLIFY(INTERFACE)
      CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
      CALL Interface_NumberOfCoupledMeshesGet(INTERFACE,numberOfCoupledMeshes,err,error,*999)
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
      NULLIFY(interfaceEquationsInterpolation)
      CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
    
      CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
        CALL InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*999)
        SELECT CASE(integrationType)
        CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
          CALL FlagError("Mesh connectivity is not implemented for frictionless contact.",err,error,*999)
        CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
          NULLIFY(pointsConnectivity)
          CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)  
          NULLIFY(interfaceMesh)
          CALL InterfacePointsConnectivity_InterfaceMeshGet(pointsConnectivity,interfaceMesh,err,error,*999)
          CALL Mesh_NumberOfDimensionsGet(interfaceMesh,numberOfInterfaceXi,err,error,*999)
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          NULLIFY(lagrangeField)
          CALL FieldVariable_FieldGet(lagrangeVariable,lagrangeField,err,error,*999)
          NULLIFY(lagrangeDecomposition)
          CALL Field_DecompositionGet(lagrangeField,lagrangeDecomposition,err,error,*999)
          NULLIFY(lagrangeDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(lagrangeDecomposition,lagrangeDecompositionTopology,err,error,*999)
          NULLIFY(dataPoints)
          CALL DecompositionTopology_DecompositionDataPointsGet(lagrangeDecompositionTopology,dataPoints,err,error,*999)
          NULLIFY(elementDataPoints)
          CALL DecompositionDataPoints_ElementDataPointsGet(dataPoints,interfaceElementNumber,elementDataPoints,err,error,*999)
          CALL DecompositionElementDataPoints_NumberOfDataPointsGet(elementDataPoints,numberOfElementDataPoints,err,error,*999)
      
          !###################################################################################################################
          
          !Test if datapoints were orthogonally projected.  
          !\todo: Allow the user to choose to only include orthogonally projected points or not (check is commented when
          !populating element matrix below).  
          ALLOCATE(orthogonallyProjected(numberOfElementDataPoints),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate orthogonal projected logicals.",err,error,*999)
          orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
          DO coupledMeshIdx=1,numberOfCoupledMeshes
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            NULLIFY(coupledInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & coupledInterpolation,err,error,*999)
            NULLIFY(coupledGeometricField)
            CALL InterfaceDomainInterpolation_GeometricFieldGet(coupledInterpolation,coupledGeometricField,err,error,*999)
            NULLIFY(coupledGeometricVariable)
            CALL Field_VariableGet(coupledGeometricField,FIELD_U_VARIABLE_TYPE,coupledGeometricVariable,err,error,*999)
            NULLIFY(coupledDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
            !mesh component number is the same for all geometric components in elasticity problems
            DO dataPointIdx=1,numberOfElementDataPoints
              CALL DecompositionElementDataPoints_GlobalNumberGet(elementDataPoints,dataPointIdx,globalDataPointNumber, &
                & err,error,*999)
              NULLIFY(pointConnectivity)
              CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,globalDataPointNumber,coupledMeshIdx, &
                & pointConnectivity,err,error,*999)
               CALL InterfacePointConnectivity_ReducedXiGet(pointConnectivity,numberOfReducedXi,reducedXi,err,error,*999)
              DO xiIdx=1,numberOfReducedXi
                IF(ABS(reducedXi(xiIdx))< ZERO_TOLERANCE) orthogonallyProjected(dataPointIdx)=.FALSE.
              ENDDO !xiIdx
            ENDDO !dataPointIdx
          ENDDO !coupledMeshIdx
          
          !Allocate memory for local allocatable variables
          ALLOCATE(gaps(numberOfElementDataPoints),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate gaps.",err,error,*999)
          gaps=0.0_DP !Initialise gap functions
          ALLOCATE(gapsComponents(3,numberOfElementDataPoints),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate component gaps.",err,error,*999)
          gapsComponents=0.0_DP !Initialise gap functions
          ALLOCATE(normals(3,numberOfElementDataPoints),STAT=err)
          IF(err/=0) CALL FlagError("Could not allocate normals.",err,error,*999)
          normals=0.0_DP !Initialise gap functions
          
          !Calculate Gap for each data point 
          !\todo: This is only required if only penetration is to penalized (ie seperation of meshes allowed.)
          ! If a no seperation condition is also required then calculation of the gap is not required.
          ! Need to allow user to choose which type of problem to solve.
          DO coupledMeshIdx=1,numberOfCoupledMeshes
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            NULLIFY(coupledInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & coupledInterpolation,err,error,*999)
            NULLIFY(coupledGeometricField)
            CALL InterfaceDomainInterpolation_GeometricFieldGet(coupledInterpolation,coupledGeometricField,err,error,*999)
            NULLIFY(coupledGeometricVariable)
            CALL Field_VariableGet(coupledGeometricField,FIELD_U_VARIABLE_TYPE,coupledGeometricVariable,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(coupledGeometricVariable,numberOfCoupledGeometricComponents,err,error,*999)
            NULLIFY(coupledDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
            NULLIFY(coupledDependentVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,coupledDependentVariable,err,error,*999)
            CALL FieldVariable_VariableTypeGet(coupledDependentVariable,coupledDependentVariableType,err,error,*999)
            NULLIFY(coupledDependentDecomposition)
            CALL Field_DecompositionGet(coupledDependentField,coupledDependentDecomposition,err,error,*999)
            NULLIFY(coupledDependentDomain)
            CALL Decomposition_DomainGet(coupledDependentDecomposition,0,coupledDependentDomain,err,error,*999)
            NULLIFY(coupledDependentDecompositionTopology)
            CALL Decomposition_DecompositionTopologyGet(coupledDependentDecomposition,coupledDependentDecompositionTopology, &
              & err,error,*999)
            NULLIFY(coupledDependentDomainTopology)
            CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
            NULLIFY(coupledDependentDecompositionElements)
            CALL DecompositionTopology_DecompositionElementsGet(coupledDependentDecompositionTopology, &
              & coupledDependentDecompositionElements,err,error,*999)
            NULLIFY(coupledDependentDecompositionElement)
            CALL DecompositionElements_ElementGet(coupledDependentDecompositionElements,localElementNumber, &
              & coupledDependentDecompositionElement,err,error,*999)
            NULLIFY(coupledDependentDomainElements)
            CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements,err,error,*999)          
            NULLIFY(coupledDependentInterpSet)
            CALL InterfaceDomainInterpolation_DependentInterpSetGet(coupledInterpolation,1,coupledDependentInterpSet, &
              & err,error,*999)
            NULLIFY(coupledDependentInterpParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(coupledDependentInterpSet,coupledDependentVariableType, &
              & coupledDependentInterpParameters,err,error,*999)
            NULLIFY(coupledDependentInterpPoint)
            CALL InterfaceInterpolationSet_InterpolatedPointGet(coupledDependentInterpSet,coupledDependentVariableType, &
              & coupledDependentInterpPoint,err,error,*999)
            NULLIFY(coupledDependentInterpPointMetrics)
            CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(coupledDependentInterpSet,coupledDependentVariableType, &
              & coupledDependentInterpPointMetrics,err,error,*999)
            DO dataPointIdx=1,numberOfElementDataPoints
              CALL DecompositionElementDataPoints_GlobalNumberGet(elementDataPoints,dataPointIdx,globalDataPointNumber, &
                & err,error,*999)
              !Only interpolate if orthogonally projected
              !\todo: Allow the user to choose to only include orthogonally projected points or not (currenlty commented out).  
              !IF(orthogonallyProjected(dataPointIdx)) THEN
              NULLIFY(pointConnectivity)
              CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,globalDataPointNumber,coupledMeshIdx, &
                & pointConnectivity,err,error,*999)
              CALL InterfacePointConnectivity_CoupledElementNumberGet(pointConnectivity,localElementNumber,err,error,*999)
              CALL InterfacePointConnectivity_CoupledLineFaceNumberGet(pointConnectivity,elementLineFaceNumber,err,error,*999)
              CALL DecompositionElement_FaceNumberGet(coupledDependentDecompositionElement,elementLineFaceNumber, &
                & localFaceLineNumber,err,error,*999)
              SELECT CASE(numberOfInterfaceXi) !Use face/line interpolation parameters for normal calculation
              CASE(1)
                CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                  & coupledDependentInterpParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              CASE(2)
                SELECT CASE(elementLineFaceNumber)
                CASE(1,3,5)
                  reverseNormal=.FALSE.
                CASE(2,4,6)
                  reverseNormal=.TRUE.
                END SELECT
                CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                  & coupledDependentInterpParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
              END SELECT
              !Determine the gap. 
              ! \todo: Note that Field_InterpolateXi(FIRST_PART_DERIV by default calculates NO_PART_DERIV too
              ! and is used because the FIRST_PART_DERIV is also need for the surface normal calculation. However the
              ! normal is only calculated for one of the coupled bodies so unnecessary computation. Need to generalize
              ! Field_InterpolateXi to allow the user to specify which PART_DERIV to calculate.
              CALL InterfacePointConnectivity_ReducedXiGet(pointConnectivity,numberOfReducedXi,reducedXi,err,error,*999)
              CALL Field_InterpolateXi(FIRST_PART_DERIV,reducedXi(1:numberOfReducedXi),coupledDependentInterpPoint, &
                & err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
              !Calculate 3 components gap function for each contact point
              gapsComponents(1:numberOfCoupledGeometricComponents,dataPointIdx)= &
                & gapsComponents(1:numberOfCoupledGeometricComponents,dataPointIdx)+ &
                & coupledDependentInterpPoint%values(1:numberOfCoupledGeometricComponents,NO_PART_DERIV)
              !Calculate surface normal (use 2nd coupled mesh surface normal)
              !\todo: Allow the user to choose which surface normal to calculate or alternatively allow for a weighted
              !       average of the two.              
              IF(coupledMeshIdx==2) THEN
                CALL Field_InterpolatedPointMetricsCalculate(numberOfCoupledGeometricComponents, &
                  & coupledDependentInterpPointMetrics,err,error,*999)
                CALL Field_PositionNormalTangentsCalculateIntPtMetric(coupledDependentInterpPointMetrics, &
                  & reverseNormal,positionPoint,normalPoint,tangentsPoint,err,error,*999)
                normals(1:numberOfCoupledGeometricComponents,dataPointIdx)=normalPoint(1:numberOfCoupledGeometricComponents)
              ENDIF !coupledMeshIdx==1
              !ENDIF !orthogonallyProjected(dataPointIdx)
            ENDDO !dataPointIdx
          ENDDO !coupledMeshIdx
          
          !###################################################################################################################
          
          !Calcualte 1 component gap
          DO dataPointIdx=1,numberOfElementDataPoints
            gaps(dataPointIdx)=DOT_PRODUCT(gapsComponents(1:numberOfCoupledGeometricComponents,dataPointIdx), &
              & normals(1:numberOfCoupledGeometricComponents,dataPointIdx))
          ENDDO !dataPointIdx
          
          !###################################################################################################################
          
          !Calculate basis function, update interface matrices with basis function, and update scale factors
          DO coupledMeshIdx=1,numberOfCoupledMeshes
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
            IF(updateMatrix) THEN
              NULLIFY(coupledElements)
              CALL InterfacePointsConnectivity_CoupledElementsGet(pointsConnectivity,interfaceElementNumber,coupledMeshIdx, &
                & coupledElements,err,error,*999)
              CALL InterfaceCoupledElements_NumberOfCoupledElementsGet(coupledElements,numberOfMatrixCoupledElements, &
                & err,error,*999)
              NULLIFY(coupledInterpolation)
              CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
                & coupledInterpolation,err,error,*999)
              NULLIFY(coupledGeometricField)
              CALL InterfaceDomainInterpolation_GeometricFieldGet(coupledInterpolation,coupledGeometricField,err,error,*999)
              NULLIFY(coupledGeometricVariable)
              CALL Field_VariableGet(coupledGeometricField,FIELD_U_VARIABLE_TYPE,coupledGeometricVariable,err,error,*999)
              CALL FieldVariable_NumberOfComponentsGet(coupledGeometricVariable,numberOfCoupledGeometricComponents,err,error,*999)
              NULLIFY(coupledDependentField)
              CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
              NULLIFY(coupledDependentVariable)
              CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,coupledDependentVariable,err,error,*999)
              CALL FieldVariable_VariableTypeGet(coupledDependentVariable,coupledDependentVariableType,err,error,*999)
              NULLIFY(coupledMesh)
              CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
              CALL Mesh_NumberOfDimensionsGet(coupledMesh,numberOfCoupledXi,err,error,*999)
              NULLIFY(coupledDependentInterpSet)
              CALL InterfaceDomainInterpolation_DependentInterpSetGet(coupledInterpolation,1,coupledDependentInterpSet, &
                & err,error,*999)
              NULLIFY(coupledDependentInterpParameters)
              CALL InterfaceInterpolationSet_InterpolationParametersGet(coupledDependentInterpSet,coupledDependentVariableType, &
              & coupledDependentInterpParameters,err,error,*999)

              DO dataPointIdx=1,numberOfElementDataPoints
                CALL DecompositionElementDataPoints_GlobalNumberGet(elementDataPoints,dataPointIdx,globalDataPointNumber, &
                  & err,error,*999)
                !\todo: Allow the user to choose gap tolerance or default to zero tolerance (currently commented out).  
                !IF(gaps(dataPointIdx)>1.0E-10) THEN !Only add contact point contribution if the gap is a penetration
                NULLIFY(pointConnectivity)
                CALL InterfacePointsConnectivity_CoupledPointGet(pointsConnectivity,globalDataPointNumber,coupledMeshIdx, &
                  & pointConnectivity,err,error,*999)
                CALL InterfacePointConnectivity_CoupledElementNumberGet(pointConnectivity,localElementNumber,err,error,*999)
                !Calculate the element index (non-conforming element) for this interface matrix
                CALL InterfaceCoupledElements_NumberOfCoupledElementsGet(coupledElements,numberOfCoupledElements,err,error,*999)
                matrixElementIdx=1
                found=.FALSE.
                continueSearch=.TRUE.
                DO WHILE(continueSearch)
                  CALL InterfaceCoupledElements_CoupledElementNumberGet(coupledElements,matrixElementIdx,coupledElementNumber, &
                    & err,error,*999)
                  IF(localElementNumber==coupledElementNumber) THEN
                    found=.TRUE.
                    continueSearch=.FALSE.
                  ELSE
                    IF(matrixElementIdx==numberOfCoupledElements) THEN
                      continueSearch=.FALSE.
                    ELSE
                      matrixElementIdx=matrixElementIdx+1
                    ENDIF
                  ENDIF
                ENDDO
                CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,coupledDependentInterpParameters, &
                  & err,error,*999)
                CALL InterfacePointConnectivity_XiGet(pointConnectivity,numberOfXi,xi(1:numberOfCoupledXi),err,error,*999)
                DO rowComponentIdx=1,numberOfCoupledGeometricComponents
                  NULLIFY(coupledDependentDomain)
                  CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                    & err,error,*999)
                  NULLIFY(coupledDependentDomainTopology)
                  CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                  NULLIFY(coupledDependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                    & err,error,*999)
                  NULLIFY(coupledDependentBasis)
                  CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,localElementNumber,coupledDependentBasis, &
                    & err,error,*999)
                  CALL Basis_NumberOfElementParametersGet(coupledDependentBasis,numberOfCoupledElementParameters, &
                    & err,error,*999)
                  DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                    CALL Basis_EvaluateXi(coupledDependentBasis,rowElementParameterIdx,NO_PART_DERIV, &
                      & xi(1:numberOfCoupledXi),rowBasisFunction,err,error,*999)
                    rowBasisFunction=rowBasisFunction*normals(rowComponentIdx,dataPointIdx)
                    rowElementDOFIdx=numberOfCoupledElementParameters*numberOfMatrixCoupledElements*(rowComponentIdx-1)+ &
                      & numberOfCoupledElementParameters*(matrixElementIdx-1)+rowElementParameterIdx
                    columnElementDOFIdx=dataPointIdx
                    !Update interface element matrix with contact point contribution
                    !\todo: Seperate multiplication of scale factors if required.  
                    interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)=rowBasisFunction* &
                      & coupledDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
                  ENDDO !rowElementParameterIdx
                ENDDO !rowComponentIdx
                !ENDIF !gaps(dataPointIdx)>ZERO_TOLERANCE
              ENDDO !dataPointIdx
              CALL InterfaceMatrix_FirstAssemblyGet(interfaceMatrix,firstAssembly,err,error,*999)
              IF(firstAssembly) CALL InterfaceMatrix_FirstAssemblySet(interfaceMatrix,.FALSE.,err,error,*999)
            ENDIF !updateMatrix
          ENDDO !coupledMeshIdx
        
          !###################################################################################################################
          
          !Deallocate memory
          IF(ALLOCATED(orthogonallyProjected)) DEALLOCATE(orthogonallyProjected)
          IF(ALLOCATED(gapsComponents)) DEALLOCATE(gapsComponents)
          IF(ALLOCATED(gaps)) DEALLOCATE(gaps)
          
          !###################################################################################################################
        
          !Calculate penalty matrix if required
          IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD) THEN
            NULLIFY(penaltyField)
            CALL InterfaceCondition_PenaltyFieldGet(interfaceCondition,penaltyField,err,error,*999)
            NULLIFY(penaltyMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,numberOfInterfaceMatrices,penaltyMatrix,err,error,*999)
            CALL InterfaceMatrix_FirstAssemblyGet(penaltyMatrix,firstAssembly,err,error,*999)
            CALL InterfaceMatrix_UpdateMatrixGet(penaltyMatrix,updateMatrix,err,error,*999)
            IF(firstAssembly.AND.updateMatrix) THEN
              NULLIFY(penaltyVariable)
              CALL Field_VariableGet(penaltyField,FIELD_U_VARIABLE_TYPE,penaltyVariable,err,error,*999)
              CALL FieldVariable_NumberOfComponentsGet(penaltyVariable,numberOfPenaltyComponents,err,error,*999)
              DO componentIdx=1,numberOfPenaltyComponents
                CALL FieldVariable_ComponentInterpolationGet(penaltyVariable,componentIdx,interpolationType,err,error,*999)
                SELECT CASE(interpolationType)
                CASE(FIELD_CONSTANT_INTERPOLATION)
                  CALL FieldVariable_ParameterSetGetConstant(penaltyVariable,FIELD_VALUES_SET_TYPE,componentIdx, &
                    & contactStiffness,err,error,*999)
                CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                  CALL FieldVariable_LocalElementDOFGet(penaltyVariable,interfaceElementNumber,componentIdx,localDOF, &
                    & err,error,*999)
                  CALL FieldVariable_ParameterSetGetLocalDOF(penaltyVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                    & contactStiffness,err,error,*999)
                  DO dataPointIdx=1,numberOfElementDataPoints
                    penaltyMatrix%elementMatrix%matrix(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                  ENDDO !dataPointIdx
                CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                  DO dataPointIdx=1,numberOfElementDataPoints
                    CALL FieldVariable_LocalDataPointDOFGet(penaltyVariable,dataPointIdx,componentIdx,localDOF,err,error,*999)
                    CALL FieldVariable_ParameterSetGetLocalDOF(penaltyVariable,FIELD_VALUES_SET_TYPE,localDOF, &
                      & contactStiffness,err,error,*999)
                    penaltyMatrix%elementMatrix%matrix(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                  ENDDO !dataPointIdx
                CASE DEFAULT
                  localError="The interpolation type for component number "// &
                    & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                    & " of variable type "//TRIM(NumberToVString(FIELD_U_VARIABLE_TYPE,"*",err,error))// &
                    & " of field number "//TRIM(NumberToVString(penaltyField%userNumber,"*",err,error))//" is "// &
                    & TRIM(NumberToVString(interpolationType,"*", err,error))//" which is invalid for penalty field."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDDO !componentIdx
            ENDIF
          ENDIF
        CASE DEFAULT
          localError="Interface condition integration type "//TRIM(NumberToVString(integrationType,"*",err,error))// &
            & " is not valid."
          CALL FlagError(localError,err,error,*999)
        END SELECT !interfaceCondition%integrationType
      CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Interface condition method "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
          & " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceConditionMethod
    ENDIF !update matrices
    
    EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FrictionlessContact_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices for the given element number for solid fluid operator.
  !>Note: First interface matrix must be solid equations set's interface matrix
  SUBROUTINE SolidFluidOperator_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)
  
    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,connectedLineFace,coupledElementNumber, &
      & coupledMeshIdx,coupledDependentScalingType,coupledDependentVariableType,decompositionFaceNumber,decompositionLineNumber, &
      & derivative,derivativeIdx,gaussPointIdx,integrationType,interfaceConditionMethod,interfaceDependentScalingType, &
      & interfaceDerivative,interfaceDerivativeIdx,interfaceMatrixIdx,interfaceNode,interfaceNodeIdx,lagrangeVariableType, &
      & localElementNode,localElementNodeIdx,localFaceNodeIdx,localLineNodeIdx,localNodeIdx, &
      & numberOfCoupledElementParameters,numberOfCoupledNodes,numberOfCoupledNodeDerivatives,numberOfCoupledXi, &
      & numberOfFaceNodes,numberOfFaceNodeDerivatives,numberOfGauss,numberOfInterfaceGeometricXi, &
      & numberOfInterfaceDependentElementParameters,numberOfInterfaceDependentNodes,numberOfInterfaceNodeDerivatives, &
      & numberOfInterfaceDependentXi,numberOfInterfaceMatrices,numberOfCoupledDependentComponents,numberOfLagrangeComponents, &
      & numberOfLineNodeDerivatives,numberOfNodeDerivatives,numberOfNodesInLine,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx
    REAL(DP) :: columnBasisFunction,jacobian,jacobianGaussWeight,gaussWeight,rowBasisFunction,xi(3)
    LOGICAL :: updateMatrix,updateMatrices
    TYPE(BasisType), POINTER :: faceBasis,interfaceDependentBasis,coupledBasis,interfaceGeometricBasis, &
      & interfaceConnectivityBasis,lineBasis
    TYPE(QuadratureSchemeType), POINTER :: interfaceQuadratureScheme
    TYPE(DecompositionType), POINTER :: coupledDependentDecomposition,interfaceDependentDecomposition, &
      & interfaceGeometricDecomposition
    TYPE(DecompositionElementType), POINTER :: coupledDependentDecompositionElement,interfaceDependentDecompositionElement
    TYPE(DecompositionElementsType), POINTER :: coupledDependentDecompositionElements,interfaceDependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: coupledDependentDecompositionTopology,interfaceDependentDecompositionTopology
    TYPE(DomainType), POINTER :: coupledDependentDomain,interfaceDependentDomain,interfaceGeometricDomain
    TYPE(DomainElementsType), POINTER :: coupledDependentDomainElements,interfaceDependentDomainElements, &
      & interfaceGeometricDomainElements
    TYPE(DomainLineType), POINTER :: coupledDomainLine
    TYPE(DomainLinesType), POINTER :: coupledDependentDomainLines
    TYPE(DomainFaceType), POINTER :: coupledDomainFace
    TYPE(DomainFacesType), POINTER :: coupledDependentDomainFaces
    TYPE(DomainTopologyType), POINTER :: coupledDependentDomainTopology,interfaceDependentDomainTopology, &
      & interfaceGeometricDomainTopology
    TYPE(FieldType), POINTER :: coupledDependentField,interfaceDependentField,interfaceGeometricField
    TYPE(FieldInterpolationParametersType), POINTER :: coupledDependentInterpParameters,interfaceGeometricInterpParameters, &
      & lagrangeInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: interfaceGeometricInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: interfaceGeometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: coupledDependentVariable,lagrangeVariable
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: interfaceInterpolation,coupledInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: coupledDependentInterpSet,interfaceGeometricInterpSet, &
      & lagrangeInterpSet
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity
    TYPE(VARYING_STRING) :: localError

    ENTERS("SolidFluidOperator_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    
    CALL InterfaceMatrices_NumberOfInterfaceMatricesGet(interfaceMatrices,numberOfInterfaceMatrices,err,error,*999)
    updateMatrices=.FALSE.
    DO interfaceMatrixIdx=1,numberOfInterfaceMatrices
      NULLIFY(interfaceMatrix)
      CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrixIdx,interfaceMatrix,err,error,*999)
      CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
      updateMatrices=(updateMatrices.OR.updateMatrix)
    ENDDO !interfaceMatrixIdx
!!TODO: Consider interface RHS

    IF(updateMatrices) THEN
    
      NULLIFY(INTERFACE)
      CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
      NULLIFY(interfaceMapping)
      CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
      !Select Interface method
      CALL InterfaceCondition_MethodGet(interfaceCondition,interfaceConditionMethod,err,error,*999)
      SELECT CASE(interfaceConditionMethod)
      CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
        !Select Integration type
        CALL InterfaceCondition_IntegrationTypeGet(interfaceCondition,integrationType,err,error,*999)
        SELECT CASE(integrationType)
        CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
          NULLIFY(meshConnectivity)
          CALL Interface_MeshConnectivityGet(INTERFACE,meshConnectivity,err,error,*999)
          NULLIFY(interfaceConnectivityBasis)
          CALL InterfaceMeshConnectivity_BasisGet(meshConnectivity,interfaceConnectivityBasis,err,error,*999)
          !Pointers to interface variables (columns of interface element matrix)
          NULLIFY(interfaceEquationsInterpolation)
          CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
          NULLIFY(interfaceInterpolation)
          CALL InterfaceEquationsInterpolation_InterfaceInterpGet(interfaceEquationsInterpolation,interfaceInterpolation, &
            & err,error,*999)
          NULLIFY(interfaceGeometricField)
          CALL InterfaceDomainInterpolation_GeometricFieldGet(interfaceInterpolation,interfaceGeometricField,err,error,*999)
          NULLIFY(interfaceGeometricDecomposition)
          CALL Field_DecompositionGet(interfaceGeometricField,interfaceGeometricDecomposition,err,error,*999)
          NULLIFY(interfaceGeometricDomain)
          CALL Decomposition_DomainGet(interfaceGeometricDecomposition,0,interfaceGeometricDomain,err,error,*999)
          NULLIFY(interfaceGeometricDomainTopology)
          CALL Domain_DomainTopologyGet(interfaceGeometricDomain,interfaceGeometricDomainTopology,err,error,*999)
          NULLIFY(interfaceGeometricDomainElements)
          CALL DomainTopology_DomainElementsGet(interfaceGeometricDomainTopology,interfaceGeometricDomainElements,err,error,*999)
          NULLIFY(interfaceGeometricBasis)
          CALL DomainElements_ElementBasisGet(interfaceGeometricDomainElements,elementNumber,interfaceGeometricBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(interfaceGeometricBasis,numberOfInterfaceGeometricXi,err,error,*999)
          NULLIFY(interfaceDependentField)
          CALL InterfaceDomainInterpolation_DependentFieldGet(interfaceInterpolation,interfaceDependentField,err,error,*999)
          NULLIFY(interfaceDependentDecomposition)
          CALL Field_DecompositionGet(interfaceDependentField,interfaceDependentDecomposition,err,error,*999)
          NULLIFY(interfaceDependentDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(interfaceDependentDecomposition,interfaceDependentDecompositionTopology, &
            & err,error,*999)
          NULLIFY(interfaceDependentDecompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(interfaceDependentDecompositionTopology, &
            & interfaceDependentDecompositionElements,err,error,*999)
          NULLIFY(interfaceDependentDomain)
          CALL Decomposition_DomainGet(interfaceDependentDecomposition,0,interfaceDependentDomain,err,error,*999)
          NULLIFY(interfaceDependentDomainTopology)
          CALL Domain_DomainTopologyGet(interfaceDependentDomain,interfaceDependentDomainTopology,err,error,*999)
          NULLIFY(interfaceDependentDomainElements)
          CALL DomainTopology_DomainElementsGet(interfaceDependentDomainTopology,interfaceDependentDomainElements,err,error,*999)
          NULLIFY(interfaceDependentBasis)
          CALL DomainElements_ElementBasisGet(interfaceDependentDomainElements,elementNumber,interfaceDependentBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(interfaceDependentBasis,numberOfInterfaceDependentXi,err,error,*999)
          CALL Basis_NumberOfLocalNodesGet(interfaceDependentBasis,numberOfInterfaceDependentNodes,err,error,*999)
          !Integrate using the interface quadrature scheme
          NULLIFY(interfaceQuadratureScheme)
          CALL Basis_QuadratureSchemeGet(interfaceGeometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,interfaceQuadratureScheme, &
            & err,error,*999)
          CALL BasisQuadratureScheme_NumberOfGaussGet(interfaceQuadratureScheme,numberOfGauss,err,error,*999)
          NULLIFY(lagrangeVariable)
          CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
          CALL FieldVariable_VariableTypeGet(lagrangeVariable,lagrangeVariableType,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(lagrangeVariable,numberOfLagrangeComponents,err,error,*999)
          NULLIFY(lagrangeInterpSet)
          CALL InterfaceDomainInterpolation_DependentInterpSetGet(interfaceInterpolation,1,lagrangeInterpSet,err,error,*999)
          NULLIFY(lagrangeInterpParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(lagrangeInterpSet,lagrangeVariableType, &
            & lagrangeInterpParameters,err,error,*999)
          !Get element interpolation param from the first geometric interp set (to get Jacobian for interface surface integral)
          NULLIFY(interfaceGeometricInterpSet)
          CALL InterfaceDomainInterpolation_GeometricInterpSetGet(interfaceInterpolation,1,interfaceGeometricInterpSet, &
            & err,error,*999)
          NULLIFY(interfaceGeometricInterpParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpParameters,err,error,*999)
          NULLIFY(interfaceGeometricInterpPoint)
          CALL InterfaceInterpolationSet_InterpolatedPointGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpPoint,err,error,*999)
          NULLIFY(interfaceGeometricInterpPointMetrics)
          CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(interfaceGeometricInterpSet,FIELD_U_VARIABLE_TYPE, &
            & interfaceGeometricInterpPointMetrics,err,error,*999)
          !Get element interpolation param from the first geometric interp set (to get Jacobian for interface surface integral)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,interfaceGeometricInterpParameters, &
            & err,error,*999)
          !Loop over interface matrices (1st solid, 2nd fluid)
          DO coupledMeshIdx=1,numberOfInterfaceMatrices
            NULLIFY(interfaceMatrix)
            CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
            CALL InterfaceMatrix_UpdateMatrixGet(interfaceMatrix,updateMatrix,err,error,*999)
            IF(updateMatrix) THEN
              !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
              NULLIFY(coupledInterpolation)
              CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
                & coupledInterpolation,err,error,*999)
              NULLIFY(coupledDependentField)
              CALL InterfaceDomainInterpolation_DependentFieldGet(coupledInterpolation,coupledDependentField,err,error,*999)
              NULLIFY(coupledDependentDecomposition)
              CALL Field_DecompositionGet(coupledDependentField,coupledDependentDecomposition,err,error,*999)
              NULLIFY(elementConnectivity)
              CALL InterfaceMeshConnectivity_CoupledElementGet(meshConnectivity,elementNumber,coupledMeshIdx, &
                & elementConnectivity,err,error,*999)
              CALL InterfaceElementConnectivity_CoupledElementNumberGet(elementConnectivity,coupledElementNumber,err,error,*999)
              CALL InterfaceElementConnectivity_CoupledLineFaceNumberGet(elementConnectivity,connectedLineFace,err,error,*999)
              CALL DecompositionElements_ElementGet(interfaceDependentDecompositionElements,coupledElementNumber, &
                & interfaceDependentDecompositionElement,err,error,*999)
              NULLIFY(coupledDependentVariable)
              CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,coupledDependentVariable,err,error,*999)
              CALL FieldVariable_VariableTypeGet(coupledDependentVariable,coupledDependentVariableType,err,error,*999)
              CALL FieldVariable_NumberOfComponentsGet(coupledDependentVariable,numberOfCoupledDependentComponents,err,error,*999)
              NULLIFY(coupledDependentInterpSet)
              CALL InterfaceDomainInterpolation_DependentInterpSetGet(coupledInterpolation,1,coupledDependentInterpSet, &
                & err,error,*999)
              NULLIFY(coupledDependentInterpParameters)
              CALL InterfaceInterpolationSet_InterpolationParametersGet(coupledDependentInterpSet,coupledDependentVariableType, &
                & coupledDependentInterpParameters,err,error,*999)
 
              !Loop over gauss points
              DO gaussPointIdx=1,numberOfGauss
                !Interpolate the geometric field at given gauss point, includes first partial derivatives
                CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & interfaceGeometricInterpPoint,err,error,*999)
                !Calculate the interpolated point metrics and the associated interpolated point
                CALL Field_InterpolatedPointMetricsCalculate(numberOfInterfaceGeometricXi,interfaceGeometricInterpPointMetrics, &
                  & err,error,*999)
                CALL FieldInterpolatedPointMetrics_JacobianGet(interfaceGeometricInterpPointMetrics,jacobian,err,error,*999)
                CALL BasisQuadratureScheme_GaussWeightGet(interfaceQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
                jacobianGaussWeight=jacobian*gaussWeight
                
                IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD.AND. &
                  & coupledMeshIdx==numberOfInterfaceMatrices) CALL FlagError("Not implemented.",err,error,*999) 
                !\todo defaults to first mesh component, generalise
                !TODO Originally xi=...
                CALL InterfaceOperators_InterfaceToCoupledMeshGaussTransform(elementConnectivity,interfaceConnectivityBasis, &
                  & gaussPointIdx,xi(1:SIZE(elementConnectivity%xi,1)),err,error,*999)
                !Loop over number of Lagrange variable components as not all components in the dependent field variable may be
                !coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                !component numbers. Generalise ordering
                DO rowComponentIdx=1,numberOfLagrangeComponents                
                  NULLIFY(coupledDependentDomain)
                  CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                    & err,error,*999)
                  NULLIFY(coupledDependentDomainTopology)
                  CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                  NULLIFY(coupledDependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                    & err,error,*999)
                  NULLIFY(coupledDependentDomainLines)
                  CALL DomainTopology_DomainLinesGet(coupledDependentDomainTopology,coupledDependentDomainLines,err,error,*999)
                  NULLIFY(coupledBasis)
                  CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,elementNumber,coupledBasis,err,error,*999)
                  CALL Basis_NumberOfXiGet(coupledBasis,numberOfCoupledXi,err,error,*999)
                  CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                  CALL Basis_NumberOfLocalNodesGet(coupledBasis,numberOfCoupledNodes,err,error,*999)
                  
                  SELECT CASE(numberOfInterfaceDependentXi)
                  
                  CASE(1) !1D interface (line)
                    NULLIFY(coupledDependentDomainLines)
                    CALL DomainTopology_DomainLinesGet(coupledDependentDomainTopology,coupledDependentDomainLines,err,error,*999)
                    CALL DecompositionElement_LineNumberGet(coupledDependentDecompositionElement,connectedLineFace, &
                      & decompositionLineNumber,err,error,*999)
                    NULLIFY(coupledDomainLine)
                    CALL DomainLines_LineGet(coupledDependentDomainLines,decompositionLineNumber,coupledDomainLine,err,error,*999)
                    NULLIFY(lineBasis)
                    CALL DomainLines_LineBasisGet(coupledDependentDomainLines,decompositionLineNumber,lineBasis,err,error,*999)
                    CALL Basis_LineNumberOfNodesGet(coupledBasis,connectedLineFace,numberOfNodesInLine,err,error,*999)
                    DO localLineNodeIdx=1,numberOfNodesInLine
                      CALL Basis_LineNodeNumberGet(coupledBasis,localLineNodeIdx,connectedLineFace,localElementNode, &
                        & err,error,*999)
                      CALL Basis_NodeNumberOfDerivativesGet(lineBasis,localLineNodeIdx,numberOfNodeDerivatives,err,error,*999)
                      DO derivativeIdx=1,numberOfNodeDerivatives
                        CALL DomainLine_DerivativeGlobalIndexGet(coupledDomainLine,derivativeIdx,localLineNodeIdx, &
                          & derivative,err,error,*999)
                        CALL Basis_ElementParameterGet(coupledBasis,derivative,localElementNode,rowElementParameterIdx, &
                          & err,error,*999)
                        !Evaluates the appropriate partial derivative index at position xi for the row basis (solid,fluid)
                        CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:2),rowBasisFunction, &
                          & err,error,*999)
                        rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                        DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                          CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx,numberOfNodeDerivatives, &
                            & err,error,*999)
                          DO interfaceDerivativeIdx=1,numberOfNodeDerivatives
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                              & columnElementParameterIdx,err,error,*999)
                            !Evaluates the appropriate partial derivative index at position xi for the row basis (lambda)
                            CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme,columnElementParameterIdx, &
                              & NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                            columnElementDOFIdx=columnElementParameterIdx+numberOfInterfaceDependentElementParameters* &
                              & (rowComponentIdx-1)
                            interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                              & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                              & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                          ENDDO !interfaceDerivativeIdx
                        ENDDO !interfaceNodeIdx
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx
                    
                  CASE(2) !2D interface (face)
                                       
                    SELECT CASE(numberOfCoupledXi)
                      
                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNodeIdx=1,numberOfCoupledNodes
                        CALL Basis_NodeNumberOfDerivativesGet(coupledBasis,localElementNodeIdx,numberOfCoupledNodeDerivatives, &
                          & err,error,*999)
                        DO derivativeIdx=1,numberOfCoupledNodeDerivatives
                          CALL Basis_ElementParameterGet(coupledBasis,derivativeIdx,localElementNodeIdx, &
                            & rowElementParameterIdx,err,error,*999)
                          CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV, &
                            & xi(1:numberOfCoupledXi),rowBasisFunction,err,error,*999)
                          rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                          DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                            CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx, &
                              & numberOfInterfaceNodeDerivatives,err,error,*999)
                            DO interfaceDerivative=1,numberOfInterfaceNodeDerivatives
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                                & columnElementParameterIdx,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme, &
                                & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                              columnElementDOFIdx=columnElementParameterIdx+numberOfInterfaceDependentElementParameters* &
                                & (rowComponentIdx-1)
                              interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                                & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                                & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                            ENDDO !interfaceDerivativeIdx
                          ENDDO !interfaceNodeIdx
                        ENDDO !derivativeIdx
                      ENDDO !localElementNodeIdx
                    
                    CASE(3) !Coupled Mesh has 3 xi directions
                      NULLIFY(coupledDependentDomainFaces)
                      CALL DomainTopology_DomainFacesGet(coupledDependentDomainTopology,coupledDependentDomainFaces,err,error,*999)
                      CALL DecompositionElement_FaceNumberGet(coupledDependentDecompositionElement,connectedLineFace, &
                        & decompositionFaceNumber,err,error,*999)
                      NULLIFY(coupledDomainFace)
                      CALL DomainFaces_FaceGet(coupledDependentDomainFaces,decompositionFaceNumber,coupledDomainFace, &
                        & err,error,*999)
                      NULLIFY(faceBasis)
                      CALL DomainFace_BasisGet(coupledDomainFace,faceBasis,err,error,*999)
                      CALL Basis_FaceNumberOfNodesGet(coupledBasis,connectedLineFace,numberOfFaceNodes,err,error,*999)
                      DO localFaceNodeIdx=1,numberOfFaceNodes
                        CALL Basis_FaceNodeNumberGet(coupledBasis,localFaceNodeIdx,connectedLineFace,localElementNode, &
                          & err,error,*999)
                        CALL Basis_NodeNumberOfDerivativesGet(faceBasis,localFaceNodeIdx,numberOfFaceNodeDerivatives,err,error,*999)
                        DO derivativeIdx=1,numberOfFaceNodeDerivatives
                          CALL Basis_FaceNodeDerivativeNumberGet(coupledBasis,derivativeIdx,localFaceNodeIdx, &
                            & connectedLineFace,derivative,err,error,*999)
                          CALL Basis_ElementParameterGet(coupledBasis,derivative,localElementNode,rowElementParameterIdx, &
                            & err,error,*999)
                          CALL Basis_EvaluateXi(coupledBasis,rowElementParameterIdx,NO_PART_DERIV,xi(1:numberOfCoupledXi), &
                            & rowBasisFunction,err,error,*999)
                          rowElementDOFIdx=rowElementParameterIdx+numberOfCoupledElementParameters*(rowComponentIdx-1)
                          DO interfaceNodeIdx=1,numberOfInterfaceDependentNodes
                            CALL Basis_NodeNumberOfDerivativesGet(interfaceDependentBasis,interfaceNodeIdx, &
                              & numberOfInterfaceNodeDerivatives,err,error,*999)
                            DO interfaceDerivativeIdx=1,numberOfInterfaceNodeDerivatives
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              CALL Basis_ElementParameterGet(interfaceDependentBasis,interfaceDerivativeIdx,interfaceNodeIdx, &
                                & columnElementParameterIdx,err,error,*999)
                              CALL BasisQuadratureScheme_GaussBasisFunctionGet(interfaceQuadratureScheme, &
                                & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnBasisFunction,err,error,*999)
                              columnElementDOFIdx=columnElementParameterIdx+numberOfInterfaceDependentElementParameters* &
                                & (rowComponentIdx-1)
                              interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                                & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                                & rowBasisFunction*columnBasisFunction*jacobianGaussWeight
                            ENDDO !interfaceDerivativeIdx
                          ENDDO !interfaceNodeIdx
                        ENDDO !derivativeIdx
                      ENDDO !faceNodeIdx                      
                    END SELECT !numberOfCoupledXi
                  END SELECT !numberOfInterfaceDependentXi
                
                ENDDO !rowComponentIdx
              ENDDO !gaussPointIdx
              
               IF(interfaceConditionMethod==INTERFACE_CONDITION_PENALTY_METHOD.AND.coupledMeshIdx==numberOfInterfaceMatrices) &
                 & CALL FlagError("Not implemented.",err,error,*999)
               
               !Scale factor adjustment
               !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix
               !contribution to the residual for non-linear problems
               !\todo update looping of variables/components for non-zero matrix elements as done above 
               !Scale factor adjustment for the Lagrange Variable (columns)
               CALL Field_ScalingTypeGet(interfaceDependentField,interfaceDependentScalingType,err,error,*999)
               CALL Field_ScalingTypeGet(coupledDependentField,coupledDependentScalingType,err,error,*999)
               IF(interfaceDependentScalingType/=FIELD_NO_SCALING) THEN
                 CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpParameters,err,error,*999)
                 rowElementDOFIdx=0
                 !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                 !\todo Currently Lagrange field variable component numbers must match each coupled dependent field
                 !variable component numbers. Generalise ordering
                 DO rowComponentIdx=1,numberOfLagrangeComponents
                   NULLIFY(coupledDependentDomain)
                   CALL FieldVariable_ComponentDomainGet(lagrangeVariable,rowComponentIdx,coupledDependentDomain,err,error,*999)
                   NULLIFY(coupledDependentDomainTopology)
                   CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                   NULLIFY(coupledDependentDomainElements)
                   CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                     & err,error,*999)
                   NULLIFY(coupledBasis)
                   CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,coupledElementNumber,coupledBasis, &
                     & err,error,*999)
                   CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                   !Loop over element rows
                   DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                     rowElementDOFIdx=rowElementDOFIdx+1
                     columnElementDOFIdx=0
                     !Loop over element columns
                     DO columnComponentIdx=1,numberOfLagrangeComponents
                       DO columnElementParameterIdx=1,numberOfInterfaceDependentElementParameters
                         columnElementDOFIdx=columnElementDOFIdx+1
                         interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                           & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                           & lagrangeInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                       ENDDO !columnElementParameterIdx
                     ENDDO !columnComponentIdx
                   ENDDO !rowElementParameterIdx
                 ENDDO !rowComponentIdx
                 !Scale factor adjustment for the row dependent variable
                 IF(coupledDependentScalingType/=FIELD_NO_SCALING) THEN
                   CALL Field_InterpolationParametersScaleFactorsElementGet(coupledElementNumber, &
                     & coupledDependentInterpParameters,err,error,*999)
                   rowElementDOFIdx=0
                   DO rowComponentIdx=1,numberOfCoupledDependentComponents
                     NULLIFY(coupledDependentDomain)
                     CALL FieldVariable_ComponentDomainGet(coupledDependentVariable,rowComponentIdx,coupledDependentDomain, &
                       & err,error,*999)
                     NULLIFY(coupledDependentDomainTopology)
                     CALL Domain_DomainTopologyGet(coupledDependentDomain,coupledDependentDomainTopology,err,error,*999)
                     NULLIFY(coupledDependentDomainElements)
                     CALL DomainTopology_DomainElementsGet(coupledDependentDomainTopology,coupledDependentDomainElements, &
                       & err,error,*999)
                     NULLIFY(coupledBasis)
                     CALL DomainElements_ElementBasisGet(coupledDependentDomainElements,coupledElementNumber,coupledBasis, &
                       & err,error,*999)
                     CALL Basis_NumberOfElementParametersGet(coupledBasis,numberOfCoupledElementParameters,err,error,*999)
                     !Loop over element rows
                     DO rowElementParameterIdx=1,numberOfCoupledElementParameters
                       rowElementDOFIdx=rowElementDOFIdx+1
                       columnElementDOFIdx=0
                       !Loop over element columns
                       DO columnComponentIdx=1,numberOfLagrangeComponents
                         DO columnElementParameterIdx=1,numberOfInterfaceDependentElementParameters
                           columnElementDOFIdx=columnElementDOFIdx+1
                           interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                             & interfaceMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                             & coupledDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
                         ENDDO !columnElementParameterIdx
                       ENDDO !columnComponentIdx
                     ENDDO !rowElementParameterIdx
                   ENDDO !rowComponentIdx
                 ENDIF
               ENDIF
             ENDIF
           ENDDO ! coupledMeshIdx
           
         CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
           CALL FlagError("Not implemented.",err,error,*999)
         CASE DEFAULT
           localError="Interface condition integration type "//TRIM(NumberToVString(integrationType,"*",err,error))//  &
             & " is not valid."
           CALL FlagError(localError,err,error,*999)
         END SELECT !interfaceCondition%integrationType
         
       CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
         CALL FlagError("Not implemented.",err,error,*999)
       CASE DEFAULT
         localError="Interface condition method "//TRIM(NumberToVString(interfaceConditionMethod,"*",err,error))// &
           & " is not valid."
         CALL FlagError(localError,err,error,*999)
       END SELECT
       
    ENDIF !update matrices

    EXITS("SolidFluidOperator_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("SolidFluidOperator_FiniteElementCalculate",err,error)
    RETURN 1
  
  END SUBROUTINE SolidFluidOperator_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  SUBROUTINE InterfaceOperators_InterfaceToCoupledMeshGaussTransform(elementConnectivity,interfaceConnectivityBasis, &
    & gaussPointIdx,transformedGaussPoint,err,error,*)
    
    !Argument variables
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity !<A pointer to the interface element connectivity for the transform
    TYPE(BasisType), POINTER :: interfaceConnectivityBasis !<The basis for the interface connectivity
   INTEGER(INTG), INTENT(IN) :: gaussPointIdx !<The gauss point index which needs to be transformed
    REAL(DP), INTENT(OUT) :: transformedGaussPoint(:) !<On exit the transformed Gauss point
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfElementParameters,numberOfXi,rowElementParameterIdx
    REAL(DP) :: gaussBasisFunction
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
#ifdef WITH_PRECHECKS
    TYPE(VARYING_STRING) :: localError
#endif

    ENTERS("InterfaceOperators_InterfaceToCoupledMeshGaussTransform",err,error,*999)
    
#ifdef WITH_PRECHECKS
    IF(.NOT.ASSOCIATED(elementConnectivity)) CALL FlagError("Element connectivity is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(interfaceConnectivityBasis)) CALL FlagError("Interface connectivity basis is not associated.", &
      & err,error,*999)
    IF(.NOT.ALLOCATED(elementConnectivity%xi)) CALL FlagError("Element connectivity xi is not allocated.",err,error,*999)
    IF(SIZE(transformedGaussPoint,1)<SIZE(elementConnectivity%xi,1)) THEN
      localError="The size of the specified transformed Gauss point array of "// &
        & TRIM(NumberToVString(SIZE(transformedGaussPoint,1),"*",err,error))// &
        & " is too small. The size of the array should be >= "// &
        & TRIM(NumberToVString(SIZE(elementConnectivity%xi,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
#endif
    
    transformedGaussPoint=0.0_DP
    numberOfXi=SIZE(elementConnectivity%xi,1)
    CALL Basis_NumberOfElementParametersGet(interfaceConnectivityBasis,numberOfElementParameters,err,error,*999)
    NULLIFY(quadratureScheme)
    CALL Basis_QuadratureSchemeGet(interfaceConnectivityBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
    DO rowElementParameterIdx=1,numberOfElementParameters
      CALL BasisQuadratureScheme_GaussBasisFunctionGet(quadratureScheme,rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx, &
        & gaussBasisFunction,err,error,*999)
      transformedGaussPoint(1:numberOfXi)=transformedGaussPoint(1:numberOfXi)+ &
        & gaussBasisFunction*elementConnectivity%xi(1:numberOfXi,1,rowElementParameterIdx)
    ENDDO !rowElementParameterIdx
     
    EXITS("InterfaceOperators_InterfaceToCoupledMeshGaussTransform")
    RETURN
999 ERRORS("InterfaceOperators_InterfaceToCoupledMeshGaussTransform",err,error)
    EXITS("InterfaceOperators_InterfaceToCoupledMeshGaussTransform")
    RETURN
    
  END SUBROUTINE InterfaceOperators_InterfaceToCoupledMeshGaussTransform

  !
  !================================================================================================================================
  !

END MODULE InterfaceOperatorsRoutines
