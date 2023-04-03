!> \file
!> \author Adam Reeve
!> \brief This module handles all Darcy pressure equations routines to be coupled with finite elasticity.
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
!> Contributor(s): Adam Reeve, Chris Bradley
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

!>This module handles all Darcy pressure equations routines.
MODULE DarcyPressureEquationsRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
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
  USE FiniteElasticityRoutines
  USE FiniteElasticityUtilityRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE ProblemAccessRoutines
  USE SolverRoutines
  USE SolverAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DarcyPressure_FiniteElementResidualEvaluate
  
  PUBLIC DarcyPressure_EquationsSetSetup
  
  PUBLIC DarcyPressure_EquationsSetSolutionMethodSet
  
  PUBLIC DarcyPressure_EquationsSetSpecificationSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element residual vector and RHS for a Darcy pressure equation finite element equations set.
  SUBROUTINE DarcyPressure_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx,esSpecification(3), &
      & fluidColsVariableType,gaussPointIdx,materialsVariableType,numberOfColsElementParameters,numberOfColsXi, &
      & numberOfDimensions,numberOfFluidColsComponents,numberOfFluidComponents,numberOfFluidXi,numberOfGauss, &
      & numberOfRowsComponents,numberOfRowsElementParameters,numberOfRowsXi,numberOfSolidColsComponents,numberOfSolidXi, &
      & numberOfXi,rowComponentIdx,rowElementDOFIdx,rowXiIdx,rowElementParameterIdx,rowsVariableType,scalingType, &
      & solidComponentNumber,solidColsVariableType,solidNumberOfXi,xiIdx   
    REAL(DP) :: colsdPhidXi(3),density,dNudXi(3,3),dNudZ(3,3),dNudZT(3,3),dXidNu(3,3),dZdNu(3,3),dZdX(3,3),fibreVectors(3,3), &
      & gaussWeight,jacobian,jacobianGaussWeight,JZ,JZNu,K(3,3),rowsdPhidXi(3),sigma(3,3),sum,tempMatrix(3,3)
    REAL(DP), POINTER :: elementResidualVector(:)
    LOGICAL :: update,updateResidual,updateRHS
    TYPE(BasisType), POINTER :: dependentBasis,fluidBasis,geometricBasis,rowsBasis,solidBasis
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(DecompositionTopologyType), POINTER :: geometricDecompositionTopology
    TYPE(DomainType), POINTER :: colsDomain,geometricDomain,rowsDomain   
    TYPE(DomainElementsType), POINTER :: fluidDomainElements,geometricDomainElements,rowsDomainElements,solidDomainElements
    TYPE(DomainTopologyType), POINTER :: colsDomainTopology,geometricDomainTopology,rowsDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: fieldVariable,fluidColsVariable,geometricVariable,materialsVariable,rowsVariable, &
      & solidColsVariable
    TYPE(FieldInterpolationParametersType), POINTER :: geometricInterpParameters,dependentInterpParameters, &
      & fibreInterpParameters,fluidDependentInterpParameters,materialsInterpParameters,rowsDependentInterpParameters, &
      & solidDependentInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpPoint,fibreInterpPoint,fluidDependentInterpPoint, &
      & materialsInterpPoint,solidDependentInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics,solidDependentInterpPointMetrics
    TYPE(QuadratureSchemeType), POINTER :: colsQuadratureScheme,fluidQuadratureScheme,geometricQuadratureScheme, &
      & quadratureScheme,rowsQuadratureScheme
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_FiniteElementResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a coupled Darcy fluid pressure type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateResidual,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    update=(updateResidual.OR.updateRHS)

    IF(update) THEN
      
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(fibreField)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      
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

      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(solidColsVariable)
      solidColsVariableType=FIELD_U_VARIABLE_TYPE
      CALL Field_VariableGet(dependentField,solidColsVariableType,solidColsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(solidColsVariable,numberOfSolidColsComponents,err,error,*999)
      NULLIFY(colsDomain)
      CALL FieldVariable_DomainGet(solidColsVariable,1,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(solidDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,solidDomainElements,err,error,*999)
      NULLIFY(solidBasis)
      CALL DomainElements_ElementBasisGet(solidDomainElements,elementNumber,solidBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(solidBasis,numberOfSolidXi,err,error,*999)
      
      NULLIFY(fluidColsVariable)
      fluidColsVariableType=FIELD_V_VARIABLE_TYPE
      CALL Field_VariableGet(dependentField,fluidColsVariableType,fluidColsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(fluidColsVariable,numberOfFluidColsComponents,err,error,*999)
      NULLIFY(colsDomain)
      CALL FieldVariable_DomainGet(solidColsVariable,1,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(fluidDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,fluidDomainElements,err,error,*999)
      NULLIFY(fluidBasis)
      CALL DomainElements_ElementBasisGet(fluidDomainElements,elementNumber,fluidBasis,err,error,*999)
      CALL Basis_NumberOfXiGet(fluidBasis,numberOfFluidXi,err,error,*999)
      NULLIFY(quadratureScheme)
      CALL Basis_QuadratureSchemeGet(fluidBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
            
      NULLIFY(materialsVariable)
      materialsVariableType=FIELD_U1_VARIABLE_TYPE
      CALL Field_VariableGet(dependentField,materialsVariableType,materialsVariable,err,error,*999)
      
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
      
      NULLIFY(rowsDependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,rowsVariableType,rowsDependentInterpParameters, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,rowsDependentInterpParameters,err,error,*999)
      
      NULLIFY(solidDependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,solidColsVariableType, &
        & solidDependentInterpParameters,err,error,*999)
      NULLIFY(solidDependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,solidColsVariableType,solidDependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,solidDependentInterpParameters, &
        & err,error,*999)
    
      NULLIFY(fluidDependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,fluidColsVariableType, &
        & fluidDependentInterpParameters,err,error,*999)
      NULLIFY(fluidDependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,fluidColsVariableType,fluidDependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fluidDependentInterpParameters, &
        & err,error,*999)
    
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,materialsVariableType,materialsInterpParameters, &
        & err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,materialsVariableType,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      IF(ASSOCIATED(fibreField)) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF

      CALL Field_ScalingTypeGet(dependentField,scalingType,err,error,*999)
      IF(scalingType/=FIELD_NO_SCALING) THEN
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,rowsDependentInterpParameters,err,error,*999)
        CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,fluidDependentInterpParameters,err,error,*999)
      ENDIF
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
         
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        IF(ASSOCIATED(fibreField)) &
          & CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,solidDependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfSolidXi,solidDependentInterpPointMetrics,err,error,*999)
        
        !Get the permeability tensor and density
!!TODO: CHECK TENSOR TYPE OF PERMEABILITY TENSOR
        CALL CoordinateSystem_MaterialTransformTensor([TENSOR_CONTRAVARIANT_INDEX,TENSOR_CONTRAVARIANT_INDEX], &
          & solidDependentInterpPointMetrics,fibreInterpPoint, &
          & materialsInterpPoint%values(1:NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),sigma,err,error,*999)        
        density=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+1,NO_PART_DERIV)
        
        !Calculate material fibre coordinate system
        CALL CoordinateSystem_MaterialFibreSystemCalculate(geometricInterpPointMetrics,fibreInterpPoint, &
          & dNudXi(1:numberOfDimensions,1:numberOfXi),dXidNu(1:numberOfXi,1:numberOfDimensions), &
          & fibreVectors(1:numberOfDimensions,1:numberOfDimensions),err,error,*999)
        !Calculate F=dZ/dNu, the deformation gradient tensor at the gauss point
        CALL FiniteElasticity_DeformationGradientTensorCalculate(solidDependentInterpPointMetrics,geometricInterpPointMetrics, &
          & dXidNu(1:numberOfXi,1:numberOfDimensions),dZdX,JZ,dZdNu,JZNu,err,error,*999)

        sigma=density*JZNu*sigma
        
        !Calculate jacobianGaussWeight.      
!!TODO: Think about symmetric problems.
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
        IF(diagnostics1) THEN
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
            & 3,3,dZdNu,WRITE_STRING_MATRIX_NAME_AND_INDICES,'(" dZdNu','(",I1,",:)',' :",3(X,E13.6))', &
            & '(17X,3(X,E13.6))',err,error,*999)
          CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Jznu",Jznu,err,error,*999)
          CALL WriteStringMatrix(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
            & 3,3,sigma,WRITE_STRING_MATRIX_NAME_AND_INDICES,'(" sigma','(",I1,",:)',' :",3(X,E13.6))', &
            & '(17X,3(X,E13.6))',err,error,*999)
        ENDIF
        
        CALL MatrixProduct(sigma,geometricInterpPointMetrics%gu,tempMatrix,err,error,*999)
        
        !Loop over field components
        rowElementDOFIdx=0
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_NumberOfXiGet(rowsBasis,numberOfRowsXi,err,error,*999)
          CALL Basis_QuadratureSchemeGet(rowsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowsQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,error,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,dependentBasis%numberOfElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            DO xiIdx=1,numberOfRowsXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowsdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(updateResidual) THEN
              sum=0.0_DP
              columnElementDOFIdx=0
              !Loop over element columns
              DO columnComponentIdx=1,numberOfFluidComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(fluidColsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(fluidDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,fluidDomainElements,err,error,*999)
                NULLIFY(fluidBasis)
                CALL DomainElements_ElementBasisGet(fluidDomainElements,elementNumber,fluidBasis,err,error,*999)
                CALL Basis_NumberOfXiGet(fluidBasis,numberOfColsXi,err,error,*999)
                CALL Basis_QuadratureSchemeGet(fluidBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,fluidQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(fluidBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  DO xiIdx=1,numberOfColsXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,colsdPhidXi(xiIdx),err,error,*999)
                  ENDDO !xiIdx
                  IF(scalingType/=FIELD_NO_SCALING) THEN
                    DO rowXiIdx=1,numberOfRowsXi
                      DO columnXiIdx=1,numberOfColsXi
                        sum=sum+rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)*tempMatrix(rowXiIdx,columnXiIdx)* &
                          & fluidDependentInterpParameters%parameters(columnElementParameterIdx,NO_PART_DERIV)* &
                          & rowsDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                          & fluidDependentInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  ELSE
                    DO rowXiIdx=1,numberOfRowsXi
                      DO columnXiIdx=1,numberOfColsXi
                        sum=sum+rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)*tempMatrix(rowXiIdx,columnXiIdx)* &
                          & fluidDependentInterpParameters%parameters(columnElementParameterIdx,NO_PART_DERIV)
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
              residualVector%elementResidual%vector(rowElementDOFIdx)=residualVector%elementResidual%vector(rowElementDOFIdx)+ &
                & sum*jacobianGaussWeight
            ENDIF
            IF(updateRHS) THEN
              IF(scalingType/=FIELD_NO_SCALING) THEN
                rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP* &
                  & rowsDependentInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
              ELSE
                rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
              ENDIF
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDDO !GaussPointIdx
            
    ENDIF !update
       
    EXITS("DarcyPressure_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("DarcyPressure_FiniteElementResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy pressure equation type of a fluid mechanics equations set class.
  SUBROUTINE DarcyPressure_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a Darcy pressure equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricMeshComponent,geometricScalingType,numberOfComponents, &
      & numberOfDarcyComponents,numberOfDimensions,numberOfSolidComponents,solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: geometricFIeld
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy pressure type of a fluid mechanics equation set."
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
        CALL DarcyPressure_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Darcy pressure equation."
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
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      numberOfComponents=numberOfDimensions !Solid components
      numberOfDarcyComponents=1 !Only solving for the fluid pressure at the moment
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          !Create the auto created dependent field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSet%dependent%dependentField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsSet%dependent%dependentField,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
          CALL Field_LabelSet(equationsSet%dependent%dependentField,"Dependent Field",err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSet%dependent%dependentField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%dependentField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%dependentField,4,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & FIELD_DP_TYPE,err,error,*999)
          
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
            & numberOfDarcyComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
            & numberOfDarcyComponents,err,error,*999)
          
          !Set labels
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
            & err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE,"V",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,"del V/del n", &
            & err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,"x1",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,2,"x2",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,3,"x3",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & "del x1/del n",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,2, &
            & "del x2/del n",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,3, &
            & "del x3/del n",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,"p",err,error,*999)
          CALL Field_ComponentLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & "del p/del n",err,error,*999)
          
          !Elasticity: Default to the geometric interpolation setup
          DO componentIdx=1,numberOfDimensions
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          !Darcy: Default pressure and mass increase to the first geometric component
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          DO componentIdx=1,numberOfDarcyComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Elasticity: Set the displacement components to node based interpolation
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Darcy: Set the pressure and mass increase components to node based interpolation
            DO componentIdx=1,numberOfDarcyComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField, &
                & FIELD_DELVDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
         
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfDarcyComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,numberOfDarcyComponents, &
            & err,error,*999)

          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Elasticity:
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Darcy:
            DO componentIdx=1,numberOfDarcyComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,componentIdx, &
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
          & " is invalid for a finite elasticity equation"
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
      numberOfComponents=NUMBER_OF_VOIGT(numberOfDimensions)+2 !permeability tensor + density + original porosity
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
        numberOfSolidComponents=6
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
        numberOfSolidComponents=4
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
        numberOfSolidComponents=6
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Darcy pressure type of a fluid mechanics equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field, which is shared with the solid equations
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,3,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)          
          !U variable, elasticity parameters
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfSolidComponents,err,error,*999)
          !V variable, solid density
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
          !U1 variable, fluid parameters
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE,numberOfComponents, &
            & err,error,*999)

          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          DO componentIdx=1,numberOfSolidComponents
            !Default the materials components to the geometric interpolation setup with constant interpolation
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          !Solid density field variable
          CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
            & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_V_VARIABLE_TYPE, &
            & 1,geometricMeshComponent,err,error,*999)
          DO componentIdx=1,numberOfComponents
            !Default the materials components to the geometric interpolation setup with constant interpolation
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx          
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,3,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
            & FIELD_U1_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfSolidComponents,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U1_VARIABLE_TYPE,numberOfComponents,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          !Default to K=I
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,5,0.0_DP,err,error,*999)
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,6,1.0_DP,err,error,*999)
          !Density
          CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U1_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,7,1.0_DP,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a standard Darcy pressure equation."
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
          & " is invalid for a finite elasticity fluid pressure equation."
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
        CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
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
          NULLIFY(vectorMapping)
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,0,err,error,*999)
          CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,2,err,error,*999)
          CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,[FIELD_V_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE], &
            & err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
              & err,error,*999)
            CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
              & err,error,*999)
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
          & " is invalid for a standard Darcy pressure equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a standard Darcy pressure equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("DarcyPressure_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("DarcyPressure_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Darcy pressure equation type of an fluid mechanics equations set class.
  SUBROUTINE DarcyPressure_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
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
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy pressure type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
      
    EXITS("DarcyPressure_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("DarcyPressure_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation specification for a Darcy pressure type of a fluid mechanics equations set.
  SUBROUTINE DarcyPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("DarcyPressure_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    subtype=specification(3)
    
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
      & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a Darcy pressure type of a fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE,subtype]

    EXITS("DarcyPressure_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("DarcyPressure_EquationsSetSpecificationSet",err,error)
    EXITS("DarcyPressure_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Darcy pressure problem post solve.
  SUBROUTINE DarcyPressure_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("DarcyPressure_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      !nothing
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy pressure type of a fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DarcyPressure_PostSolve")
    RETURN
999 ERRORSEXITS("DarcyPressure_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_PostSolve

  !
  !================================================================================================================================
  !

  !>Darcy pressure problem pre solve.
  SUBROUTINE DarcyPressure_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("DarcyPressure_PreSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
      !nothing
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Darcy pressure type of a fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("DarcyPressure_PreSolve")
    RETURN
999 ERRORSEXITS("DarcyPressure_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_PreSolve

  !
  !================================================================================================================================
  !

END MODULE DarcyPressureEquationsRoutines

