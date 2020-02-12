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
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations !<A pointer to the interface equations
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointIdx, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLineFace,decompositionLineNumber,localLineNodeIdx,decompositionFaceNumber,localFaceNodeIdx
    INTEGER(INTG) :: numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi,numberOfMatrixCoupledElements
    INTEGER(INTG) :: dataPointIdx,localElementNumber,matrixElementIdx
    INTEGER(INTG) :: matrixCoefficients(2),interfaceelementnumber
    REAL(DP) :: xi(3),jacobianGaussWeight,rowBasisFunction,columnBasisFunction,matrixCoefficient
    TYPE(BasisType), POINTER :: coupledMeshBasis,coupledMeshDependentBasis,faceBasis,interfaceDependentBasis, &
      & interfaceGeometricBasis,interfacePenaltyBasis,interfaceConnectivityBasis,lineBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition,lagrangeDecomposition,penaltyDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData
    TYPE(DecompositionElementsType), POINTER :: dependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: dependentDecompositionTopology,lagrangeDecompositionTopology
    TYPE(DomainType), POINTER :: dependentDomain,geometricDomain,penaltyDomain
    TYPE(DomainElementsType), POINTER :: dependentDomainElements,geometricDomainElements,penaltyDomainElements
    TYPE(DomainFaceType), POINTER :: coupledMeshDomainFace
    TYPE(DomainFacesType), POINTER :: dependentDomainFaces
    TYPE(DomainLineType), POINTER :: coupledMeshDomainLine
    TYPE(DomainLinesType), POINTER :: dependentDomainLines
    TYPE(DomainTopologyType), POINTER :: dependentDomainTopology,geometricDomainTopology,penaltyDomainTopology
    TYPE(FieldType), POINTER :: coupledMeshGeometricField,coupledMeshDependentField,interfaceDependentField, &
      & interfaceGeometricField,interfacePenaltyField,lagrangeField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpolationParameters,geometricInterpolationParameters, &
      & lagrangeInterpolationParameters,penaltyInterpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpolatedPoint,penaltyInterpolatedPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpolatedPointMetrics
    TYPE(FieldVariableType), POINTER :: coupledMeshGeometricVariable,interfaceMatrixVariable,lagrangeVariable
    TYPE(InterfaceType), POINTER :: INTERFACE
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: interfaceInterpolation,variableInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: dependentInterpolationSet,geometricInterpolationSet, &
      & lagrangeInterpolationSet,penaltyInterpolationSet
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix
    TYPE(InterfaceMeshConnectivityType), POINTER :: meshConnectivity     
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh     
    TYPE(QuadratureSchemeType), POINTER :: interfaceQuadratureScheme
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("FieldContinuity_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FLAG_error("Interface condition is not associated.",err,error,*999)
 
    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(interfaceMapping)
    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    NULLIFY(interfaceEquationsInterpolation)
    CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)

    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      SELECT CASE(interfaceCondition%integrationType)
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
        NULLIFY(geometricDecomposition)
        CALL Field_DecompositionGet(interfaceGeometricField,geometricDecomposition,err,error,*999)
        NULLIFY(geometricDomain)
        CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
        NULLIFY(geometricDomainTopology)
        CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
        NULLIFY(geometricDomainElements)
        CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
        NULLIFY(interfaceGeometricBasis)
        CALL DomainElements_BasisGet(geometricDomainElements,elementNumber,interfaceGeometricBasis,err,error,*999)
        NULLIFY(dependentDecomposition)
        CALL Field_DecompositionGet(interfaceDependentField,dependentDecomposition,err,error,*999)
        NULLIFY(dependentDomain)
        CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
        NULLIFY(dependentDomainTopology)
        CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
        NULLIFY(dependentDomainElements)
        CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
        NULLIFY(interfaceDependentBasis)
        CALL DomainElements_BasisGet(dependentDomainElements,elementNumber,interfaceDependentBasis,err,error,*999)
        SELECT CASE(interfaceCondition%method)
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
          CALL DomainElements_BasisGet(penaltyDomainElements,elementNumber,interfacePenaltyBasis,err,error,*999)
          NULLIFY(penaltyInterpolationSet)
          CALL InterfaceDomainInterpolation_PenaltyInterpSetGet(interfaceInterpolation,1,penaltyInterpolationSet, &
            & err,error,*999)
          NULLIFY(penaltyInterpolationParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(penaltyInterpolationSet,FIELD_U_VARIABLE_TYPE, &
            & penaltyInterpolationParameters,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,penaltyInterpolationParameters, &
            & err,error,*999)
        ENDSELECT
        !Integrate using the interface quadrature scheme
        NULLIFY(interfaceQuadratureScheme)
        CALL Basis_QuadratureSchemeGet(interfaceGeometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,interfaceQuadratureScheme, &
          & err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
        lagrangeVariableType=lagrangeVariable%variableType
        NULLIFY(lagrangeInterpolationSet)
        CALL InterfaceDomainInterpolation_DependentInterpSetGet(interfaceInterpolation,1,lagrangeInterpolationSet, &
          & err,error,*999)
        NULLIFY(lagrangeInterpolationParameters)
        CALL InterfaceInterpolationSet_InterpolationParametersGet(lagrangeInterpolationSet,lagrangeVariableType, &
          & lagrangeInterpolationParameters,err,error,*999)
        !Get element interpolation parameters from the first geometric interp set (to get Jacobian for interface surface integral)
        NULLIFY(geometricInterpolationSet)
        CALL InterfaceDomainInterpolation_GeometricInterpSetGet(interfaceInterpolation,1,geometricInterpolationSet, &
          & err,error,*999)
        NULLIFY(geometricInterpolationParameters)
        CALL InterfaceInterpolationSet_InterpolationParametersGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolationParameters,err,error,*999)
        NULLIFY(geometricInterpolatedPoint)
        CALL InterfaceInterpolationSet_InterpolatedPointGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolatedPoint,err,error,*999)
        NULLIFY(geometricInterpolatedPointMetrics)
        CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolatedPointMetrics,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
          & err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP        
        DO coupledMeshIdx=1,interfaceMatrices%numberOfInterfaceMatrices
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          IF(interfaceMatrix%updateMatrix) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            NULLIFY(variableInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & variableInterpolation,err,error,*999)
            NULLIFY(coupledMeshDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
            NULLIFY(dependentDecomposition)
            CALL Field_DecompositionGet(interfaceDependentField,dependentDecomposition,err,error,*999)
            NULLIFY(dependentDecompositionTopology)
            CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
            NULLIFY(dependentDecompositionElements)
            CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
              & err,error,*999)
            elementConnectivity=>meshConnectivity%elementConnectivity(elementNumber,coupledMeshIdx)
            coupledElementNumber=elementConnectivity%coupledElementNumber
            NULLIFY(interfaceMatrixVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,interfaceMatrixVariable,err,error,*999)
            coupledMeshVariableType=interfaceMatrixVariable%variableType
            NULLIFY(dependentInterpolationSet)
            !CALL InterfaceDomainInterpolation_DependentInterpSetGet(variableInterpolation,1,dependentInterpolationSet, &
            !  & err,error,*999)
            NULLIFY(dependentInterpolationParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(dependentInterpolationSet,coupledMeshVariableType, &
              & dependentInterpolationParameters,err,error,*999)

            !Loop over gauss points
            DO gaussPointIdx=1,interfaceQuadratureScheme%numberOfGauss
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & geometricInterpolatedPoint,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(interfaceGeometricBasis%numberOfXi,geometricInterpolatedPointMetrics, &
                & err,error,*999)
              jacobianGaussWeight=geometricInterpolatedPointMetrics%jacobian*interfaceQuadratureScheme%gaussWeights(gaussPointIdx)
              IF(interfaceCondition%method==INTERFACE_CONDITION_PENALTY_METHOD.AND. &
                & coupledMeshIdx==interfaceMatrices%numberOfInterfaceMatrices) THEN
                CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & penaltyInterpolatedPoint,err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                  !Loop over the Lagrange variable matrix rows
                  DO rowParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                    columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(rowParameterIdx,NO_PART_DERIV,gaussPointIdx)
                    rowIdx=rowIdx+1
                    interfaceMatrix%elementMatrix%matrix(rowIdx,rowIdx)= &
                      & interfaceMatrix%elementMatrix%matrix(rowIdx,rowIdx)- &
                      & (1.0_DP/penaltyInterpolatedPoint%values(1,1))*columnBasisFunction**2.0_DP*jacobianGaussWeight
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ELSE
                !\todo defaults to first mesh component, generalise
                 xi(1:SIZE(elementConnectivity%xi,1))=InterfaceOperators_InterfToCoupledMeshGaussTransform( &
                  & elementConnectivity,interfaceConnectivityBasis,gaussPointIdx,err,error)
                !XI=interfaceCondition%interface%pointsConnectivity%pointsConnectivity(gaussPointIdx,coupledMeshIdx)%xi
                !Loop over number of Lagrange variable components as not all components in the dependent field variable may
                !be coupled
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                !component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                  NULLIFY(dependentDomain)
                  CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                  NULLIFY(dependentDomainTopology)
                  CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                  NULLIFY(dependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                  NULLIFY(coupledMeshBasis)
                  CALL DomainElements_BasisGet(dependentDomainElements,coupledElementNumber,coupledMeshBasis,err,error,*999)

                  SELECT CASE(interfaceDependentBasis%numberOfXi)
                  CASE(1) !1D interface (line)
                    NULLIFY(dependentDomainLines)
                    CALL DomainTopology_DomainLinesGet(dependentDomainTopology,dependentDomainLines,err,error,*999)
                    connectedLineFace=elementConnectivity%connectedLineFace
                    decompositionLineNumber=dependentDecompositionElements%elements(coupledElementNumber)% &
                      & elementLines(connectedLineFace)
                    NULLIFY(coupledMeshDomainLine)
                    CALL DomainLines_LineGet(dependentDomainLines,decompositionLineNumber,coupledMeshDomainLine,err,error,*999)
                    NULLIFY(lineBasis)
                    CALL DomainLines_BasisGet(dependentDomainLines,decompositionLineNumber,lineBasis,err,error,*999)
                     DO localLineNodeIdx=1,coupledMeshBasis%numberOfNodesInLocalLine(connectedLineFace)
                      localElementNode=coupledMeshBasis%nodeNumbersInLocalLine(localLineNodeIdx,connectedLineFace)
                      DO derivativeIdx=1,lineBasis%numberOfDerivatives(localLineNodeIdx)
                        derivative=coupledMeshDomainLine%derivativesInLine(1,derivativeIdx,localLineNodeIdx)
                        rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                        rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,XI(1:2),err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                          DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                            columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                              & gaussPointIdx)
                            colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                            !multiplying them here
                            interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                              & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                              & rowBasisFunction*columnBasisFunction*jacobianGaussWeight*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !localLineNodeIdx

                  CASE(2) !2D interface (face)

                    NULLIFY(dependentDomainFaces)
                    CALL DomainTopology_DomainFacesGet(dependentDomainTopology,dependentDomainFaces,err,error,*999)
                    
                    SELECT CASE(coupledMeshBasis%numberOfXi)

                    CASE(2) !Coupled Mesh has 2 xi directions
                      DO localElementNode=1,coupledMeshBasis%numberOfNodes
                        DO derivative=1,coupledMeshBasis%numberOfDerivatives(localElementNode)
                          rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                          rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & xi(1:coupledMeshBasis%numberOfXi),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                            DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                              columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                                & gaussPointIdx)
                              colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                              !multiplying them here
                              interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                                & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                                & rowBasisFunction*columnBasisFunction*jacobianGaussWeight*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivative
                      ENDDO !localElementNode

                    CASE(3) !Coupled Mesh has 3 xi directions
                      connectedLineFace = elementConnectivity%connectedLineFace
                      decompositionFaceNumber=dependentDecompositionElements%elements(coupledElementNumber)% &
                        & elementFaces(connectedLineFace)
                      NULLIFY(coupledMeshDomainFace)
                      CALL DomainFaces_FaceGet(dependentDomainFaces,decompositionFaceNumber,coupledMeshDomainFace,err,error,*999)
                      NULLIFY(faceBasis)
                      CALL DomainFaces_BasisGet(dependentDomainFaces,decompositionFaceNumber,faceBasis,err,error,*999)
                      DO localFaceNodeIdx=1,coupledMeshBasis%numberOfNodesInLocalFace(connectedLineFace)
                        localElementNode=coupledMeshBasis%nodeNumbersInLocalFace(localFaceNodeIdx,connectedLineFace)
                        DO derivativeIdx=1,faceBasis%numberOfDerivatives(localFaceNodeIdx)
                          derivative=coupledMeshBasis% &
                            & derivativeNumbersInLocalFace(derivativeIdx,localFaceNodeIdx,connectedLineFace)
                          rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                          rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                            & xi(1:coupledMeshBasis%numberOfXi),err,error)
                          rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                          DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                            DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                              !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                              colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                              columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                                & gaussPointIdx)
                              colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                              !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                              !multiplying them here
                              interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                                & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                                & rowBasisFunction*columnBasisFunction*jacobianGaussWeight*matrixCoefficient
                            ENDDO !interfaceDerivative
                          ENDDO !interfaceNode
                        ENDDO !derivativeIdx
                      ENDDO !FaceNodeIdx

                    END SELECT !coupledMeshBasis%numberOfXi

                  END SELECT !interfaceDependentBasis%numberOfXi

                ENDDO !rowComponentIdx
              ENDIF
            ENDDO !gaussPointIdx

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix
            !contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%method==INTERFACE_CONDITION_PENALTY_METHOD .AND. &
              & coupledMeshIdx==interfaceEquations%interfaceMatrices%numberOfInterfaceMatrices) THEN
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpolationParameters, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                !component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                  !Loop over element Lagrange variable rows
                  DO rowParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                    rowIdx=rowIdx+1
                    interfaceMatrix%elementMatrix%matrix(rowIdx,rowIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,rowIdx)* &
                      & lagrangeInterpolationParameters%scaleFactors(rowParameterIdx,rowComponentIdx)**2
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ELSE
              !Scale factor adjustment for the Lagrange Variable (columns)
              IF(interfaceDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpolationParameters, &
                  & err,error,*999)
                rowIdx=0
                !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
                !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
                !component numbers. Generalise ordering
                DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                  NULLIFY(dependentDomain)
                  CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                  NULLIFY(dependentDomainTopology)
                  CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                  NULLIFY(dependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                  NULLIFY(coupledMeshBasis)
                  CALL DomainElements_BasisGet(dependentDomainElements,coupledElementNumber,coupledMeshBasis,err,error,*999)
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%numberOfElementParameters
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%numberOfComponents
                      DO colParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                        colIdx=colIdx+1
                        interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)* &
                          & lagrangeInterpolationParameters%scaleFactors(colParameterIdx,colComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(coupledElementNumber,dependentInterpolationParameters, &
                  & err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%numberOfComponents
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%numberOfElementParameters
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%numberOfComponents
                      DO colParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                        colIdx=colIdx+1
                        interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)* &
                        & dependentInterpolationParameters%scaleFactors(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx

      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)

        matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
        matrixCoefficients(2)=-1;
        interfaceElementNumber = elementNumber !todo simplify
        NULLIFY(pointsConnectivity)
        CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)  
        NULLIFY(interfaceMesh)
        CALL InterfacePointsConnectivity_InterfaceMeshGet(pointsConnectivity,interfaceMesh,err,error,*999)
        numberOfInterfaceMeshXi=interfaceMesh%numberOfDimensions
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
        decompositionElementData=>dataPoints%elementDataPoints(interfaceElementNumber)
        
        !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
        DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          IF(interfaceMatrix%updateMatrix) THEN
            numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
              & numberOfCoupledElements
            NULLIFY(coupledMesh)
            CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
            numberOfCoupledMeshXi=coupledMesh%numberOfDimensions

            NULLIFY(variableInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & variableInterpolation,err,error,*999)
            NULLIFY(coupledMeshGeometricField)
            CALL InterfaceDomainInterpolation_GeometricFieldGet(variableInterpolation,coupledMeshGeometricField,err,error,*999)
            NULLIFY(coupledMeshGeometricVariable)
            CALL Field_VariableGet(coupledMeshGeometricField,FIELD_U_VARIABLE_TYPE,coupledMeshGeometricVariable,err,error,*999)
            NULLIFY(coupledMeshDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
            NULLIFY(interfaceMatrixVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,interfaceMatrixVariable,err,error,*999)
            coupledMeshVariableType=interfaceMatrixVariable%variableType
            NULLIFY(dependentInterpolationSet)
            CALL InterfaceDomainInterpolation_DependentInterpSetGet(variableInterpolation,1,dependentInterpolationSet, &
              & err,error,*999)
            NULLIFY(dependentInterpolationParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(dependentInterpolationSet,coupledMeshVariableType, &
              & dependentInterpolationParameters,err,error,*999)

            numberOfCoupledMeshGeoComp=coupledMeshGeometricVariable%numberOfComponents
            DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
              localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%coupledElementNumber
              
              !Calculate the element index (non-conforming element) for this interface matrix
              matrixElementIdx=1
              DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                matrixElementIdx=matrixElementIdx+1
              ENDDO
              xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)% &
                & xi(1:numberOfCoupledMeshXi)
             DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                NULLIFY(dependentDomain)
                CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                NULLIFY(dependentDomainTopology)
                CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                NULLIFY(dependentDomainElements)
                CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                NULLIFY(coupledMeshDependentBasis)
                CALL DomainElements_BasisGet(dependentDomainElements,localElementNumber,coupledMeshDependentBasis,err,error,*999)
                
                DO rowParameterIdx=1,coupledMeshDependentBasis%numberOfElementParameters
                  rowBasisFunction=Basis_EvaluateXi(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                    & xi(1:numberOfCoupledMeshXi),err,error)*matrixCoefficients(coupledMeshIdx)
                  rowIdx=rowParameterIdx+coupledMeshDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                  colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                  interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=rowBasisFunction !Update interface element matrix with contact point contribution
                ENDDO !rowParameterIdx
              ENDDO !rowComponentIdx
            ENDDO !dataPointIdx
            
            !scale factor update
            IF(coupledMeshDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
              DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                localElementNumber=pointsConnectivity%pointsConnectivity(dataPointIdx,coupledMeshIdx)%coupledElementNumber
                CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,dependentInterpolationParameters, &
                  & err,error,*999)
                
                !Calculate the element index (non-conforming element) for this interface matrix
                matrixElementIdx=1
                DO WHILE ((localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                  & elementNumbers(matrixElementIdx)).AND.(matrixElementIdx/= &
                  & pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)%numberOfCoupledElements))
                  matrixElementIdx=matrixElementIdx+1
                ENDDO
                DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                  NULLIFY(dependentDomain)
                  CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                  NULLIFY(dependentDomainTopology)
                  CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                  NULLIFY(dependentDomainElements)
                  CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                  NULLIFY(coupledMeshDependentBasis)
                  CALL DomainElements_BasisGet(dependentDomainElements,localElementNumber,coupledMeshDependentBasis,err,error,*999)
                  DO rowParameterIdx=1,coupledMeshDependentBasis%numberOfElementParameters
                    rowIdx=rowParameterIdx+coupledMeshDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                    colIdx=dataPointIdx+decompositionElementData%numberOfProjectedData*(rowComponentIdx-1)
                    interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)* &
                      & dependentInterpolationParameters%scaleFactors(rowParameterIdx,rowComponentIdx)
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDDO !dataPointIdx
            ENDIF !.NOT. FIELD_NO_SCALING
            
          ENDIF !updateMatrix
        ENDDO !coupledMeshIdx

      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NumberToVString(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FLAG_error("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is not valid."
      CALL FLAG_error(localError,err,error,*999)
    END SELECT

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
    INTEGER(INTG) :: numberOfCoupledMeshGeoComp,numberOfInterfaceMeshXi,numberOfCoupledMeshXi, &
      & numberOfMatrixCoupledElements,localDof,elementLineFaceNumber
    INTEGER(INTG) :: dataPointIdx,coupledMeshIdx,coupledMeshVariableType,xiIdx,localElementNumber,localFaceLineNumber, &
      & matrixElementIdx,rowComponentIdx,rowParameterIdx,rowIdx,colIdx,componentIdx,globalDataPointNumber
    INTEGER(INTG) :: matrixCoefficients(2)
    REAL(DP) :: rowBasisFunction,contactStiffness
    REAL(DP) :: positionPoint(3),normalPoint(3),tangentsPoint(3,3),xi(3)
    LOGICAL :: reverseNormal
    REAL(DP), ALLOCATABLE :: gaps(:),gapsComponents(:,:),normals(:,:)
    LOGICAL, ALLOCATABLE :: orthogonallyProjected(:)       
    TYPE(BasisType), POINTER :: coupledMeshDependentBasis
    TYPE(DecompositionType), POINTER :: dependentDecomposition,lagrangeDecomposition
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(DecompositionElementsType), POINTER :: dependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: dependentDecompositionTopology,lagrangeDecompositionTopology
    TYPE(DecompositionElementDataPointsType), POINTER :: decompositionElementData 
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainElementsType), POINTER :: dependentDomainElements
    TYPE(DomainTopologyType), POINTER :: dependentDomainTopology
    TYPE(FieldType), POINTER :: coupledMeshDependentField,coupledMeshGeometricField,lagrangeField,penaltyField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpolatedPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: dependentInterpolatedPointMetrics
    TYPE(FieldVariableType), POINTER :: coupledMeshGeometricVariable,interfaceMatrixVariable,lagrangeVariable,penaltyVariable
    TYPE(InterfaceType), POINTER :: interface 
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations 
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: variableInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: dependentInterpolationSet
    TYPE(InterfaceMappingType), POINTER :: interfaceMapping
    TYPE(InterfaceMatricesType), POINTER :: interfaceMatrices
    TYPE(InterfaceMatrixType), POINTER :: interfaceMatrix,penaltyMatrix
    TYPE(InterfacePointsConnectivityType), POINTER :: pointsConnectivity 
    TYPE(MeshType), POINTER :: coupledMesh,interfaceMesh     
    TYPE(VARYING_STRING) :: localError

    ENTERS("FrictionlessContact_FiniteElementCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(interfaceCondition)) CALL FlagError("Interface condition is not associated.",err,error,*999)

    NULLIFY(interfaceEquations)
    CALL InterfaceCondition_InterfaceEquationsGet(interfaceCondition,interfaceEquations,err,error,*999)
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    NULLIFY(interfaceMapping)
    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)
    NULLIFY(interfaceEquationsInterpolation)
    CALL InterfaceEquations_EquationsInterpolationGet(interfaceEquations,interfaceEquationsInterpolation,err,error,*999)
    
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD,INTERFACE_CONDITION_PENALTY_METHOD)
      SELECT CASE(interfaceCondition%integrationType)
      CASE(INTERFACE_CONDITION_GAUSS_INTEGRATION)
        CALL FlagError("Mesh connectivity is not implemented for frictionless contact.",err,error,*999)
      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
        matrixCoefficients(1)=1; !\todo: Change to interface mapping matrix coefficients
        matrixCoefficients(2)=-1;
        NULLIFY(pointsConnectivity)
        CALL Interface_PointsConnectivityGet(INTERFACE,pointsConnectivity,err,error,*999)  
        NULLIFY(interfaceMesh)
        CALL InterfacePointsConnectivity_InterfaceMeshGet(pointsConnectivity,interfaceMesh,err,error,*999)
        numberOfInterfaceMeshXi=interfaceMesh%numberOfDimensions
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
        decompositionElementData=>dataPoints%elementDataPoints(interfaceElementNumber)
      
        !###################################################################################################################
          
        !Test if datapoints were orthogonally projected.  
        !\todo: Allow the user to choose to only include orthogonally projected points or not (check is commented when
        !populating element matrix below).  
        ALLOCATE(orthogonallyProjected(decompositionElementData%numberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate orthogonal projected logicals.",err,error,*999)
        orthogonallyProjected=.TRUE. !Initialise orthogonal projected logicals
        DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          NULLIFY(variableInterpolation)
          CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
            & variableInterpolation,err,error,*999)
          NULLIFY(coupledMeshGeometricField)
          CALL InterfaceDomainInterpolation_GeometricFieldGet(variableInterpolation,coupledMeshGeometricField,err,error,*999)
          NULLIFY(coupledMeshGeometricVariable)
          CALL Field_VariableGet(coupledMeshGeometricField,FIELD_U_VARIABLE_TYPE,coupledMeshGeometricVariable,err,error,*999)
          NULLIFY(coupledMeshDependentField)
          CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
          !mesh component number is the same for all geometric components in elasticity problems
          DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
            globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
            DO xiIdx=1,SIZE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi,1)
              IF(ABS(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%reducedXi(xiIdx)) &
                & < ZERO_TOLERANCE) THEN
                orthogonallyProjected(dataPointIdx)=.FALSE.
              ENDIF
            ENDDO !xiIdx
          ENDDO !dataPointIdx
        ENDDO !coupledMeshIdx
                
        !Allocate memory for local allocatable variables
        ALLOCATE(gaps(decompositionElementData%numberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate gaps.",err,error,*999)
        gaps=0.0_DP !Initialise gap functions
        ALLOCATE(gapsComponents(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate component gaps.",err,error,*999)
        gapsComponents=0.0_DP !Initialise gap functions
        ALLOCATE(normals(3,decompositionElementData%numberOfProjectedData),STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate normals.",err,error,*999)
        normals=0.0_DP !Initialise gap functions
        
        !Calculate Gap for each data point 
        !\todo: This is only required if only penetration is to penalized (ie seperation of meshes allowed.)
        ! If a no seperation condition is also required then calculation of the gap is not required.
        ! Need to allow user to choose which type of problem to solve.
        DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          NULLIFY(variableInterpolation)
          CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
            & variableInterpolation,err,error,*999)
          NULLIFY(coupledMeshGeometricField)
          CALL InterfaceDomainInterpolation_GeometricFieldGet(variableInterpolation,coupledMeshGeometricField,err,error,*999)
          NULLIFY(coupledMeshGeometricVariable)
          CALL Field_VariableGet(coupledMeshGeometricField,FIELD_U_VARIABLE_TYPE,coupledMeshGeometricVariable,err,error,*999)
          numberOfCoupledMeshGeoComp=coupledMeshGeometricVariable%numberOfComponents
          NULLIFY(coupledMeshDependentField)
          CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
          NULLIFY(interfaceMatrixVariable)
          CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,interfaceMatrixVariable,err,error,*999)
          coupledMeshVariableType=interfaceMatrixVariable%variableType
          NULLIFY(dependentDecomposition)
          CALL Field_DecompositionGet(coupledMeshDependentField,dependentDecomposition,err,error,*999)
          NULLIFY(dependentDomain)
          CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
          NULLIFY(dependentDecompositionTopology)
          CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
          NULLIFY(dependentDomainTopology)
          CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
          NULLIFY(dependentDecompositionElements)
          CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
            & err,error,*999)
          NULLIFY(dependentDomainElements)
          CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)          
          NULLIFY(dependentInterpolationSet)
          CALL InterfaceDomainInterpolation_DependentInterpSetGet(variableInterpolation,1,dependentInterpolationSet, &
            & err,error,*999)
          NULLIFY(dependentInterpolationParameters)
          CALL InterfaceInterpolationSet_InterpolationParametersGet(dependentInterpolationSet,coupledMeshVariableType, &
            & dependentInterpolationParameters,err,error,*999)
          NULLIFY(dependentInterpolatedPoint)
          CALL InterfaceInterpolationSet_InterpolatedPointGet(dependentInterpolationSet,coupledMeshVariableType, &
            & dependentInterpolatedPoint,err,error,*999)
          NULLIFY(dependentInterpolatedPointMetrics)
          CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(dependentInterpolationSet,coupledMeshVariableType, &
            & dependentInterpolatedPointMetrics,err,error,*999)
 
          DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
            globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
            !Only interpolate if orthogonally projected
            !\todo: Allow the user to choose to only include orthogonally projected points or not (currenlty commented out).  
            !IF(orthogonallyProjected(dataPointIdx)) THEN
            localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%coupledElementNumber
            elementLineFaceNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
              & elementLineFaceNumber
            localFaceLineNumber=dependentDecompositionElements%elements(localElementNumber)%elementFaces(elementLineFaceNumber)
            SELECT CASE(numberOfInterfaceMeshXi) !Use face/line interpolation parameters for normal calculation
            CASE(1)
              CALL Field_InterpolationParametersLineGet(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                & dependentInterpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            CASE(2)
              SELECT CASE(pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)%elementLineFaceNumber)
              CASE(1,3,5)
                reverseNormal=.FALSE.
              CASE(2,4,6)
                reverseNormal=.TRUE.
              END SELECT
              CALL Field_InterpolationParametersFaceGet(FIELD_VALUES_SET_TYPE,localFaceLineNumber, &
                & dependentInterpolationParameters,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE)
            END SELECT
            ! Determine the gap. 
            ! \todo: Note that Field_InterpolateXi(FIRST_PART_DERIV by default calculates NO_PART_DERIV too
            ! and is used because the FIRST_PART_DERIV is also need for the surface normal calculation. However the
            ! normal is only calculated for one of the coupled bodies so unnecessary computation. Need to generalize
            ! Field_InterpolateXi to allow the user to specify which PART_DERIV to calculate.  
            CALL Field_InterpolateXi(FIRST_PART_DERIV,pointsConnectivity%pointsConnectivity(globalDataPointNumber, &
              & coupledMeshIdx)%reducedXi(:),dependentInterpolatedPoint,err,error,*999,FIELD_GEOMETRIC_COMPONENTS_TYPE) !Interpolate contact data points on each surface
            gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx)=gapsComponents(1:numberOfCoupledMeshGeoComp, &
              & dataPointIdx)+dependentInterpolatedPoint%values(1:numberOfCoupledMeshGeoComp,NO_PART_DERIV)* &
              & matrixCoefficients(coupledMeshIdx) !Calculate 3 components gap function for each contact point
            !Calculate surface normal (use 2nd coupled mesh surface normal)
            !\todo: Allow the user to choose which surface normal to calculate or alternatively allow for a weighted average of the two.  
            IF (coupledMeshIdx==2) THEN
             CALL Field_InterpolatedPointMetricsCalculate(numberOfCoupledMeshGeoComp,dependentInterpolatedPointMetrics, &
                & err,error,*999)
              CALL Field_PositionNormalTangentsCalculateIntPtMetric(dependentInterpolatedPointMetrics, &
                & reverseNormal,positionPoint,normalPoint,tangentsPoint,err,error,*999)
              normals(1:numberOfCoupledMeshGeoComp,dataPointIdx)=normalPoint(1:numberOfCoupledMeshGeoComp)
            ENDIF !coupledMeshIdx==1
            !ENDIF !orthogonallyProjected(dataPointIdx)
          ENDDO !dataPointIdx
        ENDDO !coupledMeshIdx
                
        !###################################################################################################################
                
        !Calcualte 1 component gap
        DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
          gaps(dataPointIdx)=DOT_PRODUCT(gapsComponents(1:numberOfCoupledMeshGeoComp,dataPointIdx), &
            & normals(1:numberOfCoupledMeshGeoComp,dataPointIdx))
        ENDDO !dataPointIdx
        
        !###################################################################################################################
        
        !Calculate PGSMI, update interface matrices with PGSMI, and update scale factors
        DO coupledMeshIdx=1,INTERFACE%numberOfCoupledMeshes
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          IF(interfaceMatrix%updateMatrix) THEN
            numberOfMatrixCoupledElements=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
              & numberOfCoupledElements

            NULLIFY(variableInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & variableInterpolation,err,error,*999)
            NULLIFY(coupledMeshGeometricField)
            CALL InterfaceDomainInterpolation_GeometricFieldGet(variableInterpolation,coupledMeshGeometricField,err,error,*999)
            NULLIFY(coupledMeshGeometricVariable)
            CALL Field_VariableGet(coupledMeshGeometricField,FIELD_U_VARIABLE_TYPE,coupledMeshGeometricVariable,err,error,*999)
            numberOfCoupledMeshGeoComp=coupledMeshGeometricVariable%numberOfComponents
            NULLIFY(coupledMeshDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
            NULLIFY(interfaceMatrixVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,interfaceMatrixVariable,err,error,*999)
            coupledMeshVariableType=interfaceMatrixVariable%variableType
            NULLIFY(coupledMesh)
            CALL Interface_CoupledMeshGet(INTERFACE,coupledMeshIdx,coupledMesh,err,error,*999)
            numberOfCoupledMeshXi=coupledMesh%numberOfDimensions
            NULLIFY(dependentInterpolationSet)
            CALL InterfaceDomainInterpolation_DependentInterpSetGet(variableInterpolation,1,dependentInterpolationSet, &
              & err,error,*999)
            NULLIFY(dependentInterpolationParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(dependentInterpolationSet,coupledMeshVariableType, &
              & dependentInterpolationParameters,err,error,*999)

            DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
              globalDataPointNumber=decompositionElementData%dataIndices(dataPointIdx)%globalNumber
              !\todo: Allow the user to choose gap tolerance or default to zero tolerance (currently commented out).  
              !IF(gaps(dataPointIdx)>1.0E-10) THEN !Only add contact point contribution if the gap is a penetration
              localElementNumber=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                & coupledElementNumber
              !Calculate the element index (non-conforming element) for this interface matrix
              matrixElementIdx=1
              DO WHILE (localElementNumber/=pointsConnectivity%coupledElements(interfaceElementNumber,coupledMeshIdx)% &
                & elementNumbers(matrixElementIdx))
                matrixElementIdx=matrixElementIdx+1
              ENDDO
              CALL Field_InterpolationParametersScaleFactorsElementGet(localElementNumber,dependentInterpolationParameters, &
                & err,error,*999)
              xi(1:numberOfCoupledMeshXi)=pointsConnectivity%pointsConnectivity(globalDataPointNumber,coupledMeshIdx)% &
                & xi(1:numberOfCoupledMeshXi)                  
              DO rowComponentIdx=1,numberOfCoupledMeshGeoComp
                NULLIFY(dependentDomain)
                CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                NULLIFY(dependentDomainTopology)
                CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                NULLIFY(dependentDomainElements)
                CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                NULLIFY(coupledMeshDependentBasis)
                CALL DomainElements_BasisGet(dependentDomainElements,localElementNumber,coupledMeshDependentBasis,err,error,*999)
                DO rowParameterIdx=1,coupledMeshDependentBasis%numberOfElementParameters
                  rowBasisFunction=Basis_EvaluateXi(coupledMeshDependentBasis,rowParameterIdx,NO_PART_DERIV, &
                    & xi(1:numberOfCoupledMeshXi),err,error)*normals(rowComponentIdx,dataPointIdx)* &
                    & matrixCoefficients(coupledMeshIdx)
                  rowIdx=coupledMeshDependentBasis%numberOfElementParameters*numberOfMatrixCoupledElements* &
                    & (rowComponentIdx-1)+coupledMeshDependentBasis%numberOfElementParameters* &
                    & (matrixElementIdx-1)+rowParameterIdx
                  colIdx=dataPointIdx
                  !Update interface element matrix with contact point contribution
                  !\todo: Seperate multiplication of scale factors if required.  
                  interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=rowBasisFunction* &
                    & dependentInterpolationParameters%scaleFactors(rowParameterIdx,rowComponentIdx)
                ENDDO !rowParameterIdx
              ENDDO !rowComponentIdx
              !ENDIF !gaps(dataPointIdx)>ZERO_TOLERANCE
            ENDDO !dataPointIdx            
            IF(interfaceMatrix%firstAssembly) interfaceMatrix%firstAssembly=.FALSE.
          ENDIF !updateMatrix
        ENDDO !coupledMeshIdx
        
        !###################################################################################################################
        
        !Deallocate memory
        IF(ALLOCATED(orthogonallyProjected)) DEALLOCATE(orthogonallyProjected)
        IF(ALLOCATED(gapsComponents)) DEALLOCATE(gapsComponents)
        IF(ALLOCATED(gaps)) DEALLOCATE(gaps)
        
        !###################################################################################################################
        
        !Calculate penalty matrix if required
        IF(interfaceCondition%method==INTERFACE_CONDITION_PENALTY_METHOD) THEN
          NULLIFY(penaltyField)
          CALL InterfaceCondition_PenaltyFieldGet(interfaceCondition,penaltyField,err,error,*999)
          NULLIFY(penaltyMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,interfaceMatrices%numberOfInterfaceMatrices, &
            & penaltyMatrix,err,error,*999)
          IF(penaltyMatrix%firstAssembly.AND.penaltyMatrix%updateMatrix) THEN
            NULLIFY(penaltyVariable)
            CALL Field_VariableGet(penaltyField,FIELD_U_VARIABLE_TYPE,penaltyVariable,err,error,*999)
            DO componentIdx=1,penaltyVariable%numberOfComponents
              SELECT CASE(penaltyVariable%components(componentIdx)%interpolationType)
              CASE(FIELD_CONSTANT_INTERPOLATION)
                CALL Field_ParameterSetGetConstant(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & componentIdx,contactStiffness,err,error,*999)
               CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                CALL FieldVariable_LocalElementDOFGet(penaltyVariable,interfaceElementNumber,componentIdx,localDof,err,error,*999)
                CALL Field_ParameterSetGetLocalDOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & localDof,contactStiffness,err,error,*999)
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  penaltyMatrix%elementMatrix%matrix(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                ENDDO !dataPointIdx
              CASE(FIELD_DATA_POINT_BASED_INTERPOLATION)
                DO dataPointIdx=1,decompositionElementData%numberOfProjectedData
                  CALL FieldVariable_LocalDataPointDOFGet(penaltyVariable,dataPointIdx,componentIdx,localDof,err,error,*999)
                  CALL Field_ParameterSetGetLocalDOF(penaltyField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                    & localDof,contactStiffness,err,error,*999)
                  penaltyMatrix%elementMatrix%matrix(dataPointIdx,dataPointIdx)=-(1.0_DP/contactStiffness)
                ENDDO !dataPointIdx
              CASE DEFAULT
                localError="The interpolation type for component number "// &
                  & TRIM(NumberToVString(componentIdx,"*",err,error))// &
                  & " of variable type "//TRIM(NumberToVString(FIELD_U_VARIABLE_TYPE,"*",err,error))// &
                  & " of field number "//TRIM(NumberToVString(penaltyField%userNumber,"*",err,error))//" is "// &
                  & TRIM(NumberToVString(penaltyField%variableTypeMap(FIELD_U_VARIABLE_TYPE)%PTR%COMPONENTS &
                  & (componentIdx)%interpolationType,"*", err,error))// " which is invalid for penalty field."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDDO !componentIdx
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NumberToVString(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType
    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT !interfaceCondition%method
    
    EXITS("FrictionlessContact_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("FrictionlessContact_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE FrictionlessContact_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices for the given element number for solid fluid operator 
  !Note First interface matrix must be solid equations set's interface matrix
  SUBROUTINE SolidFluidOperator_FiniteElementCalculate(interfaceCondition,elementNumber,err,error,*)
  
    !Argument variables
    TYPE(InterfaceConditionType), POINTER :: interfaceCondition !<A pointer to the interface condition
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: gaussPointIdx, rowComponentIdx, rowIdx, rowParameterIdx, colComponentIdx, colIdx, colParameterIdx
    INTEGER(INTG) :: derivativeIdx,derivative,localElementNode,interfaceNode,interfaceDerivative
    INTEGER(INTG) :: coupledElementNumber,coupledMeshIdx,coupledMeshVariableType,lagrangeVariableType
    INTEGER(INTG) :: connectedLineFace,decompositionLineNumber,localLineNodeIdx,decompositionFaceNumber,localFaceNodeIdx
    REAL(DP) :: xi(3),jacobianGaussWeight,rowBasisFunction,columnBasisFunction,matrixCoefficient
    TYPE(BasisType), POINTER :: faceBasis,interfaceDependentBasis,coupledMeshBasis,interfaceGeometricBasis, &
      & interfaceConnectivityBasis,lineBasis
    TYPE(QuadratureSchemeType), POINTER :: interfaceQuadratureScheme
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DecompositionElementsType), POINTER :: dependentDecompositionElements
    TYPE(DecompositionTopologyType), POINTER :: dependentDecompositionTopology
    TYPE(DomainType), POINTER :: dependentDomain,geometricDomain
    TYPE(DomainElementsType), POINTER :: dependentDomainElements,geometricDomainElements
    TYPE(DomainLineType), POINTER :: coupledMeshDomainLine
    TYPE(DomainLinesType), POINTER :: dependentDomainLines
    TYPE(DomainFaceType), POINTER :: coupledMeshDomainFace
    TYPE(DomainFacesType), POINTER :: dependentDomainFaces
    TYPE(DomainTopologyType), POINTER :: dependentDomainTopology,geometricDomainTopology
    TYPE(FieldType), POINTER :: coupledMeshDependentField,interfaceDependentField,interfaceGeometricField
    TYPE(FieldInterpolationParametersType), POINTER :: dependentInterpolationParameters,geometricInterpolationParameters, &
      & lagrangeInterpolationParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpolatedPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpolatedPointMetrics
    TYPE(FieldVariableType), POINTER :: interfaceMatrixVariable,lagrangeVariable
    TYPE(InterfaceType), POINTER :: interface
    TYPE(InterfaceEquationsType), POINTER :: interfaceEquations
    TYPE(InterfaceEquationsDomainInterpolationType), POINTER :: interfaceInterpolation,variableInterpolation
    TYPE(InterfaceEquationsInterpolationType), POINTER :: interfaceEquationsInterpolation
    TYPE(InterfaceEquationsInterpolationSetType), POINTER :: dependentInterpolationSet,geometricInterpolationSet, &
      & lagrangeInterpolationSet
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
    NULLIFY(INTERFACE)
    CALL InterfaceCondition_InterfaceGet(interfaceCondition,INTERFACE,err,error,*999)
    NULLIFY(interfaceMapping)
    CALL InterfaceEquations_InterfaceMappingGet(interfaceEquations,interfaceMapping,err,error,*999)
    NULLIFY(interfaceMatrices)
    CALL InterfaceEquations_InterfaceMatricesGet(interfaceEquations,interfaceMatrices,err,error,*999)

    !Select Interface method
    SELECT CASE(interfaceCondition%method)
    CASE(INTERFACE_CONDITION_POINT_TO_POINT_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(INTERFACE_CONDITION_LAGRANGE_MULTIPLIERS_METHOD)
     !Select Integration type
      SELECT CASE(interfaceCondition%integrationType)
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
        NULLIFY(geometricDecomposition)
        CALL Field_DecompositionGet(interfaceGeometricField,geometricDecomposition,err,error,*999)
        NULLIFY(geometricDomain)
        CALL Decomposition_DomainGet(geometricDecomposition,0,geometricDomain,err,error,*999)
        NULLIFY(geometricDomainTopology)
        CALL Domain_DomainTopologyGet(geometricDomain,geometricDomainTopology,err,error,*999)
        NULLIFY(geometricDomainElements)
        CALL DomainTopology_DomainElementsGet(geometricDomainTopology,geometricDomainElements,err,error,*999)
        NULLIFY(interfaceGeometricBasis)
        CALL DomainElements_BasisGet(geometricDomainElements,elementNumber,interfaceGeometricBasis,err,error,*999)
        NULLIFY(interfaceDependentField)
        CALL InterfaceDomainInterpolation_DependentFieldGet(interfaceInterpolation,interfaceDependentField,err,error,*999)
        NULLIFY(dependentDecomposition)
        CALL Field_DecompositionGet(interfaceDependentField,dependentDecomposition,err,error,*999)
        NULLIFY(dependentDecompositionTopology)
        CALL Decomposition_DecompositionTopologyGet(dependentDecomposition,dependentDecompositionTopology,err,error,*999)
        NULLIFY(dependentDecompositionElements)
        CALL DecompositionTopology_DecompositionElementsGet(dependentDecompositionTopology,dependentDecompositionElements, &
          & err,error,*999)
        NULLIFY(dependentDomain)
        CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
        NULLIFY(dependentDomainTopology)
        CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
        NULLIFY(dependentDomainElements)
        CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
        NULLIFY(interfaceDependentBasis)
        CALL DomainElements_BasisGet(dependentDomainElements,elementNumber,interfaceDependentBasis,err,error,*999)
        !Integrate using the interface quadrature scheme
        NULLIFY(interfaceQuadratureScheme)
        CALL Basis_QuadratureSchemeGet(interfaceGeometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,interfaceQuadratureScheme, &
          & err,error,*999)
        NULLIFY(lagrangeVariable)
        CALL InterfaceMapping_LagrangeVariableGet(interfaceMapping,lagrangeVariable,err,error,*999)
        lagrangeVariableType=lagrangeVariable%variableType
        NULLIFY(lagrangeInterpolationSet)
        CALL InterfaceDomainInterpolation_DependentInterpSetGet(interfaceInterpolation,1,lagrangeInterpolationSet, &
          & err,error,*999)
        NULLIFY(lagrangeInterpolationParameters)
        CALL InterfaceInterpolationSet_InterpolationParametersGet(lagrangeInterpolationSet,lagrangeVariableType, &
          & lagrangeInterpolationParameters,err,error,*999)
        !Get element interpolation parameters from the first geometric interp set (to get Jacobian for interface surface integral)
        NULLIFY(geometricInterpolationSet)
        CALL InterfaceDomainInterpolation_GeometricInterpSetGet(interfaceInterpolation,1,geometricInterpolationSet, &
          & err,error,*999)
        NULLIFY(geometricInterpolationParameters)
        CALL InterfaceInterpolationSet_InterpolationParametersGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolationParameters,err,error,*999)
        NULLIFY(geometricInterpolatedPoint)
        CALL InterfaceInterpolationSet_InterpolatedPointGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolatedPoint,err,error,*999)
        NULLIFY(geometricInterpolatedPointMetrics)
        CALL InterfaceInterpolationSet_InterpolatedPointMetricsGet(geometricInterpolationSet,FIELD_U_VARIABLE_TYPE, &
          & geometricInterpolatedPointMetrics,err,error,*999)
        !Get element interpolation parameters from the first geometric interp set (to get Jacobian for interface surface integral)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,geometricInterpolationParameters, &
          & err,error,*999)
        !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
        matrixCoefficient=1.0_DP
        !Loop over interface matrices (1st solid, 2nd fluid)
        DO coupledMeshIdx=1,interfaceMatrices%numberOfInterfaceMatrices
          NULLIFY(interfaceMatrix)
          CALL InterfaceMatrices_InterfaceMatrixGet(interfaceMatrices,coupledMeshIdx,interfaceMatrix,err,error,*999)
          IF(interfaceMatrix%updateMatrix) THEN
            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of multiplying them here
            IF(coupledMeshIdx>1) THEN
              matrixCoefficient=-1.0_DP
            ENDIF 
            !Pointers to the coupledMeshIdx'th coupled mesh variables (rows of interface element matrix)
            NULLIFY(variableInterpolation)
            CALL InterfaceEquationsInterpolation_VariableInterpGet(interfaceEquationsInterpolation,coupledMeshIdx, &
              & variableInterpolation,err,error,*999)
            NULLIFY(coupledMeshDependentField)
            CALL InterfaceDomainInterpolation_DependentFieldGet(variableInterpolation,coupledMeshDependentField,err,error,*999)
            NULLIFY(dependentDecomposition)
            CALL Field_DecompositionGet(interfaceDependentField,dependentDecomposition,err,error,*999)
            elementConnectivity=>meshConnectivity%elementConnectivity(elementNumber,coupledMeshIdx)
            coupledElementNumber=elementConnectivity%coupledElementNumber
            coupledMeshVariableType=interfaceMatrixVariable%variableType
            NULLIFY(interfaceMatrixVariable)
            CALL InterfaceMapping_MatrixVariableGet(interfaceMapping,coupledMeshIdx,interfaceMatrixVariable,err,error,*999)
            coupledMeshVariableType=interfaceMatrixVariable%variableType
            NULLIFY(dependentInterpolationSet)
            CALL InterfaceDomainInterpolation_DependentInterpSetGet(variableInterpolation,1,dependentInterpolationSet, &
              & err,error,*999)
            NULLIFY(dependentInterpolationParameters)
            CALL InterfaceInterpolationSet_InterpolationParametersGet(dependentInterpolationSet,coupledMeshVariableType, &
              & dependentInterpolationParameters,err,error,*999)
 
            !Loop over gauss points
            DO gaussPointIdx=1,interfaceQuadratureScheme%numberOfGauss
              !Interpolate the geometric field at given gauss point, includes first partial derivatives
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & geometricInterpolatedPoint,err,error,*999)
              !Calculate the interpolated point metrics and the associated interpolated point
              CALL Field_InterpolatedPointMetricsCalculate(interfaceGeometricBasis%numberOfXi,geometricInterpolatedPointMetrics, &
                & err,error,*999)
              jacobianGaussWeight=geometricInterpolatedPointMetrics%jacobian*interfaceQuadratureScheme%gaussWeights(gaussPointIdx)
              IF(interfaceCondition%method==INTERFACE_CONDITION_PENALTY_METHOD.AND. &
                & coupledMeshIdx==interfaceMatrices%numberOfInterfaceMatrices) CALL FlagError("Not implemented.",err,error,*999) 
              !\todo defaults to first mesh component, generalise
              !TODO Originally xi=...
              xi(1:SIZE(elementConnectivity%xi,1))=InterfaceOperators_InterfToCoupledMeshGaussTransform(elementConnectivity, &
                & interfaceConnectivityBasis,gaussPointIdx,err,error)
              !Loop over number of Lagrange variable components as not all components in the dependent field variable may be
              !coupled
              !\todo Currently Lagrange field variable component numbers must match each coupled dependent field variable
              !component numbers. Generalise ordering
              DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                
                NULLIFY(dependentDomain)
                CALL FieldVariable_ComponentDomainGet(interfaceMatrixVariable,rowComponentIdx,dependentDomain,err,error,*999)
                NULLIFY(dependentDomainTopology)
                CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                NULLIFY(dependentDomainElements)
                CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                NULLIFY(dependentDomainLines)
                CALL DomainTopology_DomainLinesGet(dependentDomainTopology,dependentDomainLines,err,error,*999)
                NULLIFY(interfaceDependentBasis)
                CALL DomainElements_BasisGet(dependentDomainElements,elementNumber,coupledMeshBasis,err,error,*999)
                
                SELECT CASE(interfaceDependentBasis%numberOfXi)
                  
                CASE(1) !1D interface (line)
                  NULLIFY(dependentDomainLines)
                  CALL DomainTopology_DomainLinesGet(dependentDomainTopology,dependentDomainLines,err,error,*999)
                  connectedLineFace=elementConnectivity%connectedLineFace
                  decompositionLineNumber=dependentDecompositionElements%elements(coupledElementNumber)% &
                    & elementLines(connectedLineFace)                     
                  NULLIFY(coupledMeshDomainLine)
                  CALL DomainLines_LineGet(dependentDomainLines,decompositionLineNumber,coupledMeshDomainLine,err,error,*999)
                  NULLIFY(lineBasis)
                  CALL DomainLines_BasisGet(dependentDomainLines,decompositionLineNumber,lineBasis,err,error,*999)
                  DO localLineNodeIdx=1,coupledMeshBasis%numberOfNodesInLocalLine(connectedLineFace)
                    localElementNode=coupledMeshBasis%nodeNumbersInLocalLine(localLineNodeIdx,connectedLineFace)
                    DO derivativeIdx=1,lineBasis%numberOfDerivatives(localLineNodeIdx)                     
                      derivative=coupledMeshDomainLine%derivativesInLine(1,derivativeIdx,localLineNodeIdx)                    
                      rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                      !Evaluates the appropriate partial derivative index at position xi for the row basis (solid,fluid)
                      rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV,xi(1:2),err,error)
                      rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                      DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                        DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                          !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                          colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                          !Evaluates the appropriate partial derivative index at position xi for the row basis (lambda)
                          columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                            & gaussPointIdx)
                          colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                          !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                          !multiplying them here
                          interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                            & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                            & columnBasisFunction*rowBasisFunction*jacobianGaussWeight*matrixCoefficient
                        ENDDO !interfaceDerivative
                      ENDDO !interfaceNode
                    ENDDO !derivativeIdx
                  ENDDO !localLineNodeIdx
                  
                CASE(2) !2D interface (face)
                  
                  NULLIFY(dependentDomainFaces)
                  CALL DomainTopology_DomainFacesGet(dependentDomainTopology,dependentDomainFaces,err,error,*999)
                  
                  SELECT CASE(coupledMeshBasis%numberOfXi)
                    
                  CASE(2) !Coupled Mesh has 2 xi directions
                    DO localElementNode=1,coupledMeshBasis%numberOfNodes
                      DO derivative=1,coupledMeshBasis%numberOfDerivatives(localElementNode)
                        rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                        rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                          & xi(1:coupledMeshBasis%numberOfXi),err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                          DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                            columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                              & gaussPointIdx)
                            colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                            !multiplying them here
                            interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                              & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                              & columnBasisFunction*rowBasisFunction*jacobianGaussWeight*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivative
                    ENDDO !localElementNode
                    
                  CASE(3) !Coupled Mesh has 3 xi directions
                    connectedLineFace = elementConnectivity%connectedLineFace
                    decompositionFaceNumber=dependentDecompositionElements%elements(coupledElementNumber)% &
                      & elementFaces(connectedLineFace)
                    NULLIFY(coupledMeshDomainFace)
                    CALL DomainFaces_FaceGet(dependentDomainFaces,decompositionFaceNumber,coupledMeshDomainFace,err,error,*999)
                    NULLIFY(faceBasis)
                    CALL DomainFaces_BasisGet(dependentDomainFaces,decompositionFaceNumber,faceBasis,err,error,*999)
                    DO localFaceNodeIdx=1,coupledMeshBasis%numberOfNodesInLocalFace(connectedLineFace)
                      localElementNode=coupledMeshBasis%nodeNumbersInLocalFace(localFaceNodeIdx,connectedLineFace)
                      DO derivativeIdx=1,faceBasis%numberOfDerivatives(localFaceNodeIdx)
                        derivative=coupledMeshBasis%derivativeNumbersInLocalFace(derivativeIdx,localFaceNodeIdx, &
                          & connectedLineFace)
                        rowParameterIdx=coupledMeshBasis%elementParameterIndex(derivative,localElementNode)
                        rowBasisFunction=Basis_EvaluateXi(coupledMeshBasis,rowParameterIdx,NO_PART_DERIV, &
                          & XI(1:coupledMeshBasis%numberOfXi),err,error)
                        rowIdx=rowParameterIdx+coupledMeshBasis%numberOfElementParameters*(rowComponentIdx-1)
                        DO interfaceNode=1,interfaceDependentBasis%numberOfNodes
                          DO interfaceDerivative=1,interfaceDependentBasis%numberOfDerivatives(interfaceNode)
                            !\todo requires equal number of nodes between interface mesh and coupled mesh. Generalize
                            colParameterIdx=interfaceDependentBasis%elementParameterIndex(interfaceDerivative,interfaceNode)
                            columnBasisFunction=interfaceQuadratureScheme%gaussBasisFunctions(colParameterIdx,NO_PART_DERIV, &
                              & gaussPointIdx)
                            colIdx=colParameterIdx+interfaceDependentBasis%numberOfElementParameters*(rowComponentIdx-1)
                            !\todo Use matrix coefficients in solver routines when assembling solver matrices instead of
                            !multiplying them here
                            interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)= &
                              & interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)+ &
                              & rowBasisFunction*columnBasisFunction*jacobianGaussWeight*matrixCoefficient
                          ENDDO !interfaceDerivative
                        ENDDO !interfaceNode
                      ENDDO !derivativeIdx
                    ENDDO !faceNodeIdx
                    
                  END SELECT !coupledMeshBasis%numberOfXi

                END SELECT !interfaceDependentBasis%numberOfXi
                
              ENDDO !rowComponentIdx
            ENDDO !gaussPointIdx

            !Scale factor adjustment
            !\todo check if scale factor adjustments are already made elsewhere eg when calculating the interface matrix
            !contribution to the residual for non-linear problems
            !\todo update looping of variables/components for non-zero matrix elements as done above 
            IF(interfaceCondition%method==INTERFACE_CONDITION_PENALTY_METHOD.AND. &
              & coupledMeshIdx==interfaceEquations%interfaceMatrices%numberOfInterfaceMatrices) &
              & CALL FlagError("Not implemented.",err,error,*999)
            !Scale factor adjustment for the Lagrange Variable (columns)
            IF(interfaceDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
              CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,lagrangeInterpolationParameters, &
                & err,error,*999)
              rowIdx=0
              !Use Lagrange variable number of components here since we are only dealing with Lagrange variable scale factors 
              !\todo Currently Lagrange field variable component numbers must match each coupled dependent field
              !variable component numbers. Generalise ordering
              DO rowComponentIdx=1,lagrangeVariable%numberOfComponents
                NULLIFY(dependentDomain)
                CALL FieldVariable_ComponentDomainGet(lagrangeVariable,rowComponentIdx,dependentDomain,err,error,*999)
                NULLIFY(dependentDomainTopology)
                CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
                NULLIFY(dependentDomainElements)
                CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
                NULLIFY(interfaceDependentBasis)
                CALL DomainElements_BasisGet(dependentDomainElements,coupledElementNumber,coupledMeshBasis,err,error,*999)
                !Loop over element rows
                DO rowParameterIdx=1,coupledMeshBasis%numberOfElementParameters
                  rowIdx=rowIdx+1
                  colIdx=0
                  !Loop over element columns
                  DO colComponentIdx=1,lagrangeVariable%numberOfComponents
                    DO colParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                      colIdx=colIdx+1
                      interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)* &
                        & lagrangeInterpolationParameters%scaleFactors(colParameterIdx,colComponentIdx)
                    ENDDO !colParameterIdx
                  ENDDO !colComponentIdx
                ENDDO !rowParameterIdx
              ENDDO !rowComponentIdx
              !Scale factor adjustment for the row dependent variable
              IF(coupledMeshDependentField%scalings%scalingType/=FIELD_NO_SCALING) THEN
                CALL Field_InterpolationParametersScaleFactorsElementGet(coupledElementNumber,dependentInterpolationParameters, &
                  & err,error,*999)
                rowIdx=0
                DO rowComponentIdx=1,interfaceMatrixVariable%numberOfComponents
                  !Loop over element rows
                  DO rowParameterIdx=1,coupledMeshBasis%numberOfElementParameters
                    rowIdx=rowIdx+1
                    colIdx=0
                    !Loop over element columns
                    DO colComponentIdx=1,lagrangeVariable%numberOfComponents
                      DO colParameterIdx=1,interfaceDependentBasis%numberOfElementParameters
                        colIdx=colIdx+1
                        interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)=interfaceMatrix%elementMatrix%matrix(rowIdx,colIdx)* &
                        & dependentInterpolationParameters%scaleFactors(rowParameterIdx,rowComponentIdx)
                      ENDDO !colParameterIdx
                    ENDDO !colComponentIdx
                  ENDDO !rowParameterIdx
                ENDDO !rowComponentIdx
              ENDIF
            ENDIF
          ENDIF
        ENDDO ! coupledMeshIdx
        
      CASE(INTERFACE_CONDITION_DATA_POINTS_INTEGRATION)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Interface condition integration type "//TRIM(NumberToVString(interfaceCondition%integrationType, &
          & "*",err,error))// " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT !interfaceCondition%integrationType

    CASE(INTERFACE_CONDITION_AUGMENTED_LAGRANGE_METHOD)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Interface condition method "//TRIM(NumberToVString(interfaceCondition%method,"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("SolidFluidOperator_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("SolidFluidOperator_FiniteElementCalculate",err,error)
    RETURN 1
  
  END SUBROUTINE SolidFluidOperator_FiniteElementCalculate
  
  !
  !================================================================================================================================
  !

  FUNCTION InterfaceOperators_InterfToCoupledMeshGaussTransform(elementConnectivity,interfaceConnectivityBasis,gaussPointIdx, &
    & err,error)
  
    !Argument variables
    TYPE(InterfaceElementConnectivityType), POINTER :: elementConnectivity !<A pointer to the element connectivity
    TYPE(BasisType), POINTER :: interfaceConnectivityBasis !<A pointer to the interface mesh connectivity basis
    INTEGER(INTG), INTENT(IN) :: gaussPointIdx !< Index to the gauss point which needs to be transformed
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Function variable
    REAL(DP) :: InterfaceOperators_InterfToCoupledMeshGaussTransform(SIZE(elementConnectivity%xi,1))
    !Local Variables
    INTEGER(INTG) :: rowParameterIdx
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme

    ENTERS("InterfaceOperators_InterfToCoupledMeshGaussTransform",err,error,*999)
    
    InterfaceOperators_InterfToCoupledMeshGaussTransform=0.0_DP
    NULLIFY(quadratureScheme)
    CALL Basis_QuadratureSchemeGet(interfaceConnectivityBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
    DO rowParameterIdx = 1,interfaceConnectivityBasis%numberOfElementParameters
      InterfaceOperators_InterfToCoupledMeshGaussTransform(:)= InterfaceOperators_InterfToCoupledMeshGaussTransform(:) + &
        & quadratureScheme%gaussBasisFunctions(rowParameterIdx,NO_PART_DERIV,gaussPointIdx)* &
        & elementConnectivity%xi(:,1,rowParameterIdx)
    ENDDO !rowParameterIdx
     
    EXITS("InterfaceOperators_InterfToCoupledMeshGaussTransform")
    RETURN
999 ERRORS("InterfaceOperators_InterfToCoupledMeshGaussTransform",err,error)
    EXITS("InterfaceOperators_InterfToCoupledMeshGaussTransform")
    RETURN
    
  END FUNCTION InterfaceOperators_InterfToCoupledMeshGaussTransform

END MODULE InterfaceOperatorsRoutines
