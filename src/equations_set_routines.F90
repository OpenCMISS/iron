!!> \file
!> \author Chris Bradley
!> \brief This module handles all equations set routines.
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

!>This module handles all equations set routines.
MODULE EquationsSetRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BioelectricRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE ClassicalFieldRoutines
  USE CmissMPI
  USE ComputationRoutines
  USE ComputationAccessRoutines
  USE Constants
  USE CoordinateSystemAccessRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE ElasticityRoutines
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FittingRoutines
  USE FluidMechanicsRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
  USE MeshAccessRoutines
  USE MonodomainEquationsRoutines
#ifndef NOMPIMOD
  USE MPI
#endif
  USE MultiPhysicsRoutines
  USE ProfilingRoutines
  USE RegionAccessRoutines
  USE Strings
  USE Timer
  USE Types

#include "macros.h"  

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC EquationsSet_AnalyticCreateStart,EquationsSet_AnalyticCreateFinish

  PUBLIC EquationsSet_AnalyticDestroy

  PUBLIC EquationsSet_AnalyticEvaluate

  PUBLIC EquationsSet_Assemble
  
  PUBLIC EquationsSet_Backsubstitute
  
  PUBLIC EquationsSet_BoundaryConditionsAnalytic

  PUBLIC EquationsSet_CreateStart,EquationsSet_CreateFinish

  PUBLIC EquationsSet_Destroy

  PUBLIC EquationsSets_Finalise,EquationsSets_Initialise

  PUBLIC EquationsSet_EquationsCreateFinish,EquationsSet_EquationsCreateStart

  PUBLIC EquationsSet_EquationsDestroy
  
  PUBLIC EquationsSet_MaterialsCreateStart,EquationsSet_MaterialsCreateFinish

  PUBLIC EquationsSet_MaterialsDestroy
  
  PUBLIC EquationsSet_DependentCreateStart,EquationsSet_DependentCreateFinish

  PUBLIC EquationsSet_DependentDestroy

  PUBLIC EquationsSet_DerivedCreateStart,EquationsSet_DerivedCreateFinish

  PUBLIC EquationsSet_DerivedDestroy
  
  PUBLIC EquationSet_IndependentCreateStart,EquationsSet_IndependentCreateFinish

  PUBLIC EquationsSet_IndependentDestroy
  
  PUBLIC EquationsSet_JacobianEvaluate,EquationsSet_ResidualEvaluate

  PUBLIC EquationsSet_OutputTypeSet
  
  PUBLIC EquationsSet_SolutionMethodSet
  
  PUBLIC EquationsSet_SourceCreateStart,EquationsSet_SourceCreateFinish

  PUBLIC EquationsSet_SourceDestroy

  PUBLIC EquationsSet_TensorInterpolateGaussPoint

  PUBLIC EquationsSet_TensorInterpolateXi

  PUBLIC EquationsSet_DerivedVariableCalculate,EquationsSet_DerivedVariableSet

  PUBLIC EquationsSet_LoadIncrementApply
  
  PUBLIC EquationsSet_AnalyticUserParamSet,EquationsSet_AnalyticUserParamGet

  PUBLIC EquationsSet_TimesSet

CONTAINS

  !
  !================================================================================================================================
  !
      
  !>Finish the creation of a analytic solution for equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_AnalyticCreateFinish
  SUBROUTINE EquationsSet_AnalyticCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to create the analytic for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: analyticField

    ENTERS("EquationsSet_AnalyticCreateFinish",err,error,*999)

    CALL EquationsSet_AssertAnalyticNotFinished(equationsSet,err,error,*999)

    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    analyticField=>equationsSet%analytic%analyticField
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_UserNumberGet(analyticField,equationsSetSetupInfo%fieldUserNumber,err,error,*999)
      equationsSetSetupInfo%field=>analyticField
    ENDIF
    !Finish the equations set specific analytic setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish the analytic creation
    equationsSet%analytic%analyticFinished=.TRUE.
      
    EXITS("EquationsSet_AnalyticCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of a analytic solution for a equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_AnalyticCreateStart
  SUBROUTINE EquationsSet_AnalyticCreateStart(equationsSet,analyticFunctionType,analyticFieldUserNumber,analyticField, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of an analytic for.
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The analytic function type to setup \see EquationsSetRoutines_AnalyticFunctionTypes,EquationsSetRoutines
    INTEGER(INTG), INTENT(IN) :: analyticFieldUserNumber !<The user specified analytic field number
    TYPE(FieldType), POINTER :: analyticField !<If associated on entry, a pointer to the user created analytic field which has the same user number as the specified analytic field user number. If not associated on entry, on exit, a pointer to the created analytic field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,fieldUserNumber
    TYPE(DecompositionType), POINTER :: analyticDecomposition,geometricDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,analyticFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_AnalyticCreateStart",err,error,*998)

    CALL EquationsSet_AssertAnalyticNotCreated(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"analytic",analyticFieldUserNumber,analyticField,err,error,*998)
    
    !Initialise the equations set analytic
    CALL EquationsSet_AnalyticInitialise(equationsSet,err,error,*999)
    equationsSet%analytic%analyticFieldAutoCreated=(.NOT.ASSOCIATED(analyticField))
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=analyticFieldUserNumber
    equationsSetSetupInfo%field=>analyticField
    equationsSetSetupInfo%analyticFunctionType=analyticFunctionType
    !Start the equations set specific analytic setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(equationsSet%analytic%analyticFieldAutoCreated) THEN
      analyticField=>equationsSet%analytic%analyticField
    ELSE
      equationsSet%analytic%analyticField=>analyticField
    ENDIF
    
    EXITS("EquationsSet_AnalyticCreateStart")
    RETURN
999 CALL EquationsSet_AnalyticFinalise(equationsSet%analytic,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_AnalyticCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticCreateStart
  
  !
  !================================================================================================================================
  !
  
  !>Destroy the analytic solution for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_AnalyticDestroy
  SUBROUTINE EquationsSet_AnalyticDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy the analytic solutins for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticDestroy",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)

    CALL EquationsSet_AnalyticFinalise(equationsSet%analytic,err,error,*999)
        
    EXITS("EquationsSet_AnalyticDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticDestroy

  !
  !================================================================================================================================
  !

  !>Evaluates the current analytic solution for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_AnalyticEvaluate
  SUBROUTINE EquationsSet_AnalyticEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the current analytic solutins for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentInterpolationType,componentIdx,derivativeIdx,elementIdx,gaussPointIdx, &
      & globalDerivativeIndex,nodeIdx,numberOfAnalyticComponents,numberOfComponents,numberOfDimensions,numberOfElements, &
      & numberOfGauss,numberOfNodes,numberOfNodeDerivatives,numberOfVariables,numberOfVersions,variableIdx,variableType,versionIdx
    REAL(DP) :: analyticTime,normal(3),position(3),tangents(3,3),analyticValue
    REAL(DP) :: analyticDummyValues(1)=0.0_DP
    REAL(DP) :: materialsDummyValues(1)=0.0_DP
    LOGICAL :: reverseNormal=.FALSE.
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: analyticInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: analyticInterpPoint,geometricInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldPhysicalPointType), POINTER :: analyticPhysicalPoint,materialsPhysicalPoint
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,geometricVariable,materialsVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_AnalyticEvaluate",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsFinished(equationsSet,err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    CALL FieldVariable_InterpolationParameterInitialise(geometricVariable,geometricInterpParameters,err,error,*999)
    CALL Field_InterpolatedPointInitialise(geometricInterpParameters,geometricInterpPoint,err,error,*999)
    CALL Field_InterpolatedPointMetricsInitialise(geometricInterpPoint,geometricInterpPointMetrics,err,error,*999)
    NULLIFY(analyticField)
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    IF(ASSOCIATED(analyticField)) THEN
      NULLIFY(analyticVariable)
      CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(analyticVariable,numberOfAnalyticComponents,err,error,*999)
      CALL FieldVariable_InterpolationParameterInitialise(analyticVariable,analyticInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointInitialise(analyticInterpParameters,analyticInterpPoint,err,error,*999)
      CALL Field_PhysicalPointInitialise(analyticInterpPoint,geometricInterpPoint,analyticPhysicalPoint,err,error,*999)
    ENDIF
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    IF(ASSOCIATED(materialsField)) THEN
      NULLIFY(materialsVariable)
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(materialsVariable,numberOfAnalyticComponents,err,error,*999)
      CALL FieldVariable_InterpolationParameterInitialise(materialsVariable,materialsInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointInitialise(materialsInterpParameters,materialsInterpPoint,err,error,*999)
      CALL Field_PhysicalPointInitialise(materialsInterpPoint,geometricInterpPoint,materialsPhysicalPoint,err,error,*999)
    ENDIF
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    CALL EquationsSet_AnalyticTimeGet(equationsSet,analyticTime,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        CALL FieldVariable_ComponentInterpolationGet(dependentVariable,componentIdx,componentInterpolationType,err,error,*999)
        SELECT CASE(componentInterpolationType)
        CASE(FIELD_CONSTANT_INTERPOLATION)
          CALL FlagError("Cannot evaluate an analytic solution for a constant interpolation components.",err,error,*999)
        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          !Loop over the local elements excluding the ghosts
          CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
          DO elementIdx=1,numberOfElements
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
            IF(ASSOCIATED(analyticField)) THEN
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,analyticInterpParameters, &
                & err,error,*999)
            ENDIF
            IF(ASSOCIATED(materialsField)) THEN
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,materialsInterpParameters, &
                & err,error,*999)
            ENDIF
            CALL Field_InterpolateXi(FIRST_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP],geometricInterpPoint,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE,geometricInterpPointMetrics,err,error,*999)
            CALL Field_PositionNormalTangentsCalculateIntPtMetric(geometricInterpPointMetrics,reverseNormal,position,normal, &
              & tangents,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolateXi(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP],analyticInterpPoint, &
              & err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolateXi(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP],materialsInterpPoint, &
              & err,error,*999)
!! \todo Maybe do this with optional arguments?
            IF(ASSOCIATED(analyticField)) THEN
              IF(ASSOCIATED(materialsField)) THEN
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticInterpPoint%values(:,NO_PART_DERIV), &
                  & materialsInterpPoint%values(:,NO_PART_DERIV),analyticValue,err,error,*999)
              ELSE
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticInterpPoint%values(:,NO_PART_DERIV), &
                  & materialsDummyValues,analyticValue,err,error,*999)
              ENDIF
            ELSE
              IF(ASSOCIATED(materialsField)) THEN
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues,materialsInterpPoint% &
                  & values(:,NO_PART_DERIV),analyticValue,err,error,*999)
              ELSE
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues,materialsDummyValues, &
                  & analyticValue,err,error,*999)
              ENDIF
            ENDIF
            CALL FieldVariable_ParameterSetUpdateLocalElement(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,elementIdx, &
              & componentIdx,analyticValue,err,error,*999)
          ENDDO !elementIdx
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL Field_PositionNormalTangentsCalculateNode(dependentField,variableType,componentIdx,nodeIdx,position,normal, &
              & tangents,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolateFieldVariableNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE, &
                & analyticVariable,componentIdx,nodeIdx,analyticPhysicalPoint,err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolateFieldVariableNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE, &
              & materialsVariable,componentIdx,nodeIdx,materialsPhysicalPoint,err,error,*999)
            !Loop over the derivatives
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives                                
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
!! \todo Maybe do this with optional arguments?
              IF(ASSOCIATED(analyticField)) THEN
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticPhysicalPoint%values, &
                    & materialsPhysicalPoint%values,analyticValue,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticPhysicalPoint%values, &
                    & materialsDummyValues,analyticValue,err,error,*999)
                ENDIF
              ELSE
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues, &
                    & materialsPhysicalPoint%values,analyticValue,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues,materialsDummyValues, &
                    & analyticValue,err,error,*999)
                ENDIF
              ENDIF
              !Loop over the versions
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_ParameterSetUpdateLocalNode(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,versionIdx, &
                  & derivativeIdx,nodeIdx,componentIdx,analyticValue,err,error,*999)
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          !Loop over the local elements excluding the ghosts
          CALL DomainElements_NumberOfElementsGet(domainElements,numberOfElements,err,error,*999)
          DO elementIdx=1,numberOfElements
            NULLIFY(basis)
            CALL DomainElements_ElementBasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx,geometricInterpParameters,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & analyticInterpParameters,err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & materialsInterpParameters,err,error,*999)
            !Loop over the Gauss points in the element
            NULLIFY(quadratureScheme)
            CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
            CALL BasisQuadratureScheme_NumberOfGaussGet(quadratureScheme,numberOfGauss,err,error,*999)
            DO gaussPointIdx=1,numberOfGauss
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
                & err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE,geometricInterpPointMetrics,err,error,*999)
              CALL Field_PositionNormalTangentsCalculateIntPtMetric(geometricInterpPointMetrics,reverseNormal,position,normal, &
                & tangents,err,error,*999)
              IF(ASSOCIATED(analyticField)) CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME, &
                & gaussPointIdx,analyticInterpPoint,err,error,*999)
              IF(ASSOCIATED(materialsField)) CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME, &
                & gaussPointIdx,materialsInterpPoint,err,error,*999)
!! \todo Maybe do this with optional arguments?
              IF(ASSOCIATED(analyticField)) THEN
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticInterpPoint%values(:,NO_PART_DERIV), &
                    & materialsInterpPoint%values(:,NO_PART_DERIV),analyticValue,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticInterpPoint%values(:,NO_PART_DERIV), &
                    & materialsDummyValues,analyticValue,err,error,*999)
                ENDIF
              ELSE
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues, &
                    & materialsInterpPoint%values(:,NO_PART_DERIV),analyticValue,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivativeIndex,componentIdx,analyticDummyValues,materialsDummyValues, &
                    & analyticValue,err,error,*999)
                ENDIF
              ENDIF
              CALL FieldVariable_ParameterSetUpdateLocalGaussPoint(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                & gaussPointIdx,elementIdx,componentIdx,analyticValue,err,error,*999)
            ENDDO !gaussPointIdx
          ENDDO !elementIdx
        CASE DEFAULT
          localError="The interpolation type of "//TRIM(NumberToVString(componentInterpolationType,"*",err,error))// &
            & " for component "//TRIM(NumberToVString(componentIdx,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    IF(ASSOCIATED(materialsField)) THEN
      CALL Field_PhysicalPointFinalise(materialsPhysicalPoint,err,error,*999)
      CALL Field_InterpolatedPointFinalise(materialsInterpPoint,err,error,*999)
      CALL FieldVariable_InterpolationParameterFinalise(materialsInterpParameters,err,error,*999)
    ENDIF
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_PhysicalPointFinalise(analyticPhysicalPoint,err,error,*999)
      CALL Field_InterpolatedPointFinalise(analyticInterpPoint,err,error,*999)
      CALL FieldVariable_InterpolationParameterFinalise(analyticInterpParameters,err,error,*999)
    ENDIF
    CALL Field_InterpolatedPointMetricsFinalise(geometricInterpPointMetrics,err,error,*999)
    CALL Field_InterpolatedPointFinalise(geometricInterpPoint,err,error,*999)
    CALL FieldVariable_InterpolationParameterFinalise(geometricInterpParameters,err,error,*999)
           
    EXITS("EquationsSet_AnalyticEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticEvaluate

  !
  !================================================================================================================================
  !

  !>Finalise the analytic solution for an equations set and deallocate all memory.
  SUBROUTINE EquationsSet_AnalyticFinalise(equationsSetAnalytic,err,error,*)

    !Argument variables
    TYPE(EquationsSetAnalyticType), POINTER :: equationsSetAnalytic!<A pointer to the equations set analytic to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_AnalyticFinalise",err,error,*999)

    IF(ASSOCIATED(equationsSetAnalytic)) DEALLOCATE(equationsSetAnalytic)
       
    EXITS("EquationsSet_AnalyticFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticFinalise

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for an equations set.
  SUBROUTINE EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal,time, &
    & variableType,globalDerivative,componentNumber,analyticParameters,materialsParameters,analyticValue,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the analytic for
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: position(:) !<position(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivative !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: analyticValue !<On return, the analtyic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_AnalyticFunctionsEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      IF(SIZE(equationsSet%specification,1)<2) CALL FlagError("Equations set specification must have at least two "// &
        & "entries for a classical field equations set.",err,error,*999)
      CALL ClassicalField_AnalyticFunctionsEvaluate(equationsSet,equationsSet%specification(2),analyticFunctionType,position, &
        & tangents,normal,time,variableType,globalDerivative,componentNumber,analyticParameters,materialsParameters, &
        & analyticValue,err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsSet_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticFunctionsEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Initialises the analytic solution for an equations set.
  SUBROUTINE EquationsSet_AnalyticInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("EquationsSet_AnalyticInitialise",err,error,*998)

    CALL EquationsSet_AssertAnalyticNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%analytic,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set analytic.",err,error,*999)
    equationsSet%analytic%equationsSet=>equationsSet
    equationsSet%analytic%analyticFinished=.FALSE.
    equationsSet%analytic%analyticFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%analytic%analyticField)
    equationsSet%analytic%analyticTime=0.0_DP
       
    EXITS("EquationsSet_AnalyticInitialise")
    RETURN
999 CALL EquationsSet_AnalyticFinalise(equationsSet%analytic,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_AnalyticInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticInitialise

  !
  !================================================================================================================================
  !

  !>Sets the analytic problem user parameter
  SUBROUTINE EquationsSet_AnalyticUserParamSet(equationsSet,parameterIdx,parameter,err,error,*)
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(IN) :: parameterIdx !<Index of the user parameter
    REAL(DP), INTENT(IN) :: parameter !<Value of the parameter
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_AnalyticUserParamSet",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
 
    IF(parameterIdx<1.OR.parameterIdx>SIZE(equationsSet%analytic%analyticUserParams,1)) THEN
      localError="The specified parameter index of "//TRIM(NumberToVString(parameterIdx,"*",err,error))// &
        & " is invalid. The parameter index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(equationsSet%analytic%analyticUserParams,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Set the value
    equationsSet%analytic%analyticUserParams(parameterIdx)=parameter

    EXITS("EquationsSet_AnalyticUserParamSet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticUserParamSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticUserParamSet

  !
  !================================================================================================================================
  !

  !>Sets the analytic problem user parameter
  SUBROUTINE EquationsSet_AnalyticUserParamGet(equationsSet,parameterIdx,parameter,err,error,*)
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(IN) :: parameterIdx !<Index of the user parameter
    REAL(DP), INTENT(OUT) :: parameter !<Value of the parameter
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_AnalyticUserParamGet",err,error,*999)
    
    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
 
    IF(parameterIdx<1.OR.parameterIdx>SIZE(equationsSet%analytic%analyticUserParams,1)) THEN
      localError="The specified parameter index of "//TRIM(NumberToVString(parameterIdx,"*",err,error))// &
        & " is invalid. The parameter index should be >= 1 and <= "// &
        & TRIM(NumberToVString(SIZE(equationsSet%analytic%analyticUserParams,1),"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    !Set the value    
    PARAMETER=equationsSet%analytic%analyticUserParams(parameterIdx)

    EXITS("EquationsSet_AnalyticUserParamGet")
    RETURN
999 ERRORSEXITS("EquationsSet_AnalyticUserParamGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AnalyticUserParamGet

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an equations set.
  SUBROUTINE EquationsSet_Assemble(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearity,timeDependence
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_Assemble",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_AssertIsFinished(equations,err,error,*999)
    
    IF(equationsSet%outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set assemble: ",equationsSet%label,err,error,*999)
    ENDIF

    CALL Equations_TimeDependenceTypeGet(equations,timeDependence,err,error,*999)
    CALL Equations_LinearityTypeGet(equations,linearity,err,error,*999)
    SELECT CASE(timeDependence)
    CASE(EQUATIONS_STATIC)
      SELECT CASE(linearity)
      CASE(EQUATIONS_LINEAR)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleStaticLinearFEM(equationsSet,err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleStaticNonlinearFEM(equationsSet,err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL EquationsSet_AssembleStaticNonlinearNodal(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity of "//TRIM(NumberToVString(linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_QUASISTATIC)
      SELECT CASE(linearity)
      CASE(EQUATIONS_LINEAR)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleQuasistaticLinearFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR)
        CALL EquationsSet_AssembleQuasistaticNonlinearFEM(equationsSet,err,error,*999)
      CASE(EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations linearity of "//TRIM(NumberToVString(linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(linearity)
      CASE(EQUATIONS_LINEAR)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleDynamicLinearFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_NONLINEAR_BCS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations set linearity of "//TRIM(NumberToVString(linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_TIME_STEPPING)
      CALL FlagError("Time stepping equations are not assembled.",err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(timeDependence,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsSet_Assemble")
    RETURN
999 ERRORSEXITS("EquationsSet_Assemble",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_Assemble

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a dynamic linear equations set using the finite element method.
  SUBROUTINE EquationsSet_AssembleDynamicLinearFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalStart,internalFinish,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    
    ENTERS("EquationsSet_AssembleDynamicLinearFEM",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
    
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AssembleDynamicLinearFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_AssembleDynamicLinearFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssembleDynamicLinearFEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a linear static equations set using the finite element method.
  SUBROUTINE EquationsSet_AssembleStaticLinearFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalStart,internalFinish,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    
    ENTERS("EquationsSet_AssembleStaticLinearFEM",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
 
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AssembleStaticLinearFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_AssembleStaticLinearFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssembleStaticLinearFEM

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear static equations set using the finite element method.
  SUBROUTINE EquationsSet_AssembleStaticNonlinearFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalStart,internalFinish,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    
    ENTERS("EquationsSet_AssembleStaticNonlinearFEM",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)

    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AssembleStaticNonlinearFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_AssembleStaticNonlinearFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssembleStaticNonlinearFEM

  !
  !================================================================================================================================
  !

  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear quasistatic equations set using the finite
  !>element method. Currently the same as the static nonlinear case
  SUBROUTINE EquationsSet_AssembleQuasistaticNonlinearFEM(equationsSet,err,error,*)
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_AssembleQuasistaticNonlinearFEM",err,error,*999)

    !Currently no difference
    CALL EquationsSet_AssembleStaticNonlinearFEM(equationsSet,err,error,*999)
    
    EXITS("EquationsSet_AssembleQuasistaticNonlinearFEM")
    RETURN
999 ERRORS("EquationsSet_AssembleQuasistaticNonlinearFEM",err,error)
    EXITS("EquationsSet_AssembleQuasistaticNonlinearFEM")
    RETURN 1
    
  END  SUBROUTINE EquationsSet_AssembleQuasistaticNonlinearFEM

  !
  !================================================================================================================================
  !

  !>Assembles the equations stiffness matrix and rhs for a linear quasistatic equations set using the finite element method.
  SUBROUTINE EquationsSet_AssembleQuasistaticLinearFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalStart,internalFinish,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    
    ENTERS("EquationsSet_AssembleQuasistaticLinearFEM",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
   
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AssembleQuasistaticLinearFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_AssembleQuasistaticLinearFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssembleQuasistaticLinearFEM

  !
  !================================================================================================================================
  !

  !>Backsubstitutes with an equations set to calculate unknown right hand side vectors
  SUBROUTINE EquationsSet_Backsubstitute(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to backsubstitute
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<The boundary conditions to use for the backsubstitution
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dampingMatrixNumber,dynamicMatrixIdx,equationsColumnIdx,equationsColumnNumber,equationsMatrixIdx, &
      & equationsRowNumber,equationsStorageType,lhsBoundaryCondition,lhsVariableDOF,linearMatrixIdx,massMatrixNumber, &
      & numberOfDirichletConditions,numberOfDynamicMatrices,numberOfLinearMatrices,numberOfResiduals,numberOfRows, &
      & numberOfSources,residualIdx,rhsBoundaryCondition,rhsGlobalDOF,rhsVariableDOF,rhsVariableType,rowCondition, &
      & sourceIdx,stiffnessMatrixNumber,variableDOF,variableType
    INTEGER(INTG), POINTER :: equationsRowToLHSDOFMap(:),equationsRowToRHSDOFMap(:)
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:)
    REAL(DP) :: dampingAlpha,dependentValue,lhsSum,linearAlpha,massAlpha,matrixValue,residualAlpha,residualValue,rhsValue, &
      & sourceAlpha,sourceValue,stiffnessAlpha
    TYPE(BoundaryConditionsVariableType), POINTER :: lhsBoundaryConditionsVariable
    TYPE(BoundaryConditionsRowVariableType), POINTER :: lhsBoundaryConditionsRowVariable
    TYPE(DomainMappingType), POINTER :: columnDomainMapping,rhsDomainMapping
    TYPE(DistributedMatrixType), POINTER :: dampingDistributedMatrix,linearDistributedMatrix,massDistributedMatrix, &
      & stiffnessDistributedMatrix
    TYPE(DistributedVectorType), POINTER :: accelerationDistributedVector,displacementDistributedVector,linearDistributedVector, &
      & residualDistributedVector,sourceDistributedVector,velocityDistributedVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,linearMatrix,massMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldParameterSetType), POINTER :: accelerationParameters,displacementParameters,linearParameters,velocityParameters
    TYPE(FieldVariableType), POINTER :: dynamicVariable,lhsVariable,linearVariable,rhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_Backsubstitute",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)
   
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)

    IF(ASSOCIATED(rhsMapping)) THEN
      !Only back substitute if we have a RHS variable
 
      NULLIFY(rhsVariable)
      CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
      NULLIFY(equationsRowToRHSDOFMap)
      CALL EquationsMappingRHS_EquationsRowToRHSDOFMapGet(rhsMapping,equationsRowToRHSDOFMap,err,error,*999)
      
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(lhsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,lhsVariable,err,error,*999)
      CALL EquationsMappingLHS_NumberOfRowsGet(lhsMapping,numberOfRows,err,error,*999)
      NULLIFY(equationsRowToLHSDOFMap)
      CALL EquationsMappingLHS_EquationsRowTOLHSDOFMapGet(lhsMapping,equationsRowToLHSDOFMap,err,error,*999)
      
      NULLIFY(lhsBoundaryConditionsRowVariable)
      CALL BoundaryConditions_RowVariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsRowVariable,err,error,*999)
      NULLIFY(lhsBoundaryConditionsVariable)
      CALL BoundaryConditions_VariableGet(boundaryConditions,lhsVariable,lhsBoundaryConditionsVariable,err,error,*999)
      CALL BoundaryConditionsVariable_NumberOfDirichletConditionsGet(lhsBoundaryConditionsVariable,numberOfDirichletConditions, &
        & err,error,*999)
      
      IF(numberOfDirichletConditions>0) THEN
      
        NULLIFY(vectorMatrices)
        CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)

        NULLIFY(dynamicMapping)
        CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
        IF(ASSOCIATED(dynamicMapping)) THEN
          NULLIFY(dynamicMatrices)
          CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
          CALL EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*999)
          CALL EquationsMappingDynamic_StiffnessMatrixNumberGet(dynamicMapping,stiffnessMatrixNumber,err,error,*999)
          CALL EquationsMappingDynamic_DampingMatrixNumberGet(dynamicMapping,dampingMatrixNumber,err,error,*999)
          CALL EquationsMappingDynamic_MassMatrixNumberGet(dynamicMapping,massMatrixNumber,err,error,*999)
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          IF(stiffnessMatrixNumber/=0) THEN
            NULLIFY(displacementParameters)
            CALL FieldVariable_ParameterSetGet(dynamicVariable,FIELD_VALUES_SET_TYPE,displacementParameters,err,error,*999)
            NULLIFY(displacementDistributedVector)
            CALL FieldParameterSet_ParametersGet(displacementParameters,displacementDistributedVector,err,error,*999)
            NULLIFY(stiffnessMatrix)
            CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,stiffnessMatrixNumber,stiffnessMatrix,err,error,*999)
            NULLIFY(stiffnessDistributedMatrix)
            CALL EquationsMatrix_DistributedMatrixGet(stiffnessMatrix,stiffnessDistributedMatrix,err,error,*999)
            CALL EquationsMatrix_MatrixCoefficientGet(stiffnessMatrix,stiffnessAlpha,err,error,*999)
          ENDIF
          IF(dampingMatrixNumber/=0) THEN
            NULLIFY(velocityParameters)
            CALL FieldVariable_ParameterSetExists(dynamicVariable,FIELD_VELOCITY_VALUES_SET_TYPE,velocityParameters,err,error,*999)
            !If no velocity values then use the incremental alpha
            IF(.NOT.ASSOCIATED(velocityParameters)) &
              & CALL FieldVariable_ParameterSetGet(dynamicVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,velocityParameters, &
              & err,error,*999)
            NULLIFY(velocityDistributedVector)
            CALL FieldParameterSet_ParametersGet(velocityParameters,velocityDistributedVector,err,error,*999)
            NULLIFY(dampingMatrix)
            CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,dampingMatrixNumber,dampingMatrix,err,error,*999)
            NULLIFY(dampingDistributedMatrix)
            CALL EquationsMatrix_DistributedMatrixGet(dampingMatrix,dampingDistributedMatrix,err,error,*999)
            CALL EquationsMatrix_MatrixCoefficientGet(dampingMatrix,dampingAlpha,err,error,*999)
          ENDIF
          IF(massMatrixNumber/=0) THEN
            NULLIFY(accelerationParameters)
            CALL FieldVariable_ParameterSetExists(dynamicVariable,FIELD_ACCELERATION_VALUES_SET_TYPE,accelerationParameters, &
              & err,error,*999)
            !If no accerlation values then use the incremental alpha
            IF(ASSOCIATED(accelerationParameters)) &
              & CALL FieldVariable_ParameterSetGet(dynamicVariable,FIELD_INCREMENTAL_VALUES_SET_TYPE,accelerationParameters, &
              & err,error,*999)
            NULLIFY(accelerationDistributedVector)
            CALL FieldParameterSet_ParametersGet(accelerationParameters,accelerationDistributedVector,err,error,*999)
            NULLIFY(massMatrix)
            CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,massMatrixNumber,massMatrix,err,error,*999)
            NULLIFY(massDistributedMatrix)
            CALL EquationsMatrix_DistributedMatrixGet(massMatrix,massDistributedMatrix,err,error,*999)
            CALL EquationsMatrix_MatrixCoefficientGet(massMatrix,massAlpha,err,error,*999)
          ENDIF
        ENDIF
        NULLIFY(linearMapping)
        CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
        IF(ASSOCIATED(linearMapping)) THEN
          NULLIFY(linearMatrices)
          CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
          CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
        ENDIF
        NULLIFY(nonlinearMapping)
        CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
        IF(ASSOCIATED(nonlinearMapping)) THEN
          NULLIFY(nonlinearMatrices)
          CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
          CALL EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*999)
        ENDIF
        NULLIFY(sourcesMapping)
        CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
        IF(ASSOCIATED(sourcesMapping)) THEN
          NULLIFY(sourceVectors)
          CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
          CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
        ENDIF

        !Loop over the rows in the equations set
        DO equationsRowNumber=1,numberOfRows
          CALL BoundaryConditionsRowVariable_RowConditionTypeGet(lhsBoundaryConditionsRowVariable,equationsRowNumber, &
            & rowCondition,err,error,*999)
          lhsVariableDOF=equationsRowToLHSDOFMap(equationsRowNumber)
          CALL BoundaryConditionsVariable_DOFTypeGet(lhsBoundaryConditionsVariable,lhsVariableDOF,lhsBoundaryCondition, &
            & err,error,*999)

          rhsVariableDOF=equationsRowToRHSDOFMap(equationsRowNumber)

          SELECT CASE(rowCondition)
          CASE(BOUNDARY_CONDITION_FREE_ROW)
            !OK, do nothing
          CASE(BOUNDARY_CONDITION_DIRICHLET_ROW)
            !Backsubtitute to find RHS value
            lhsSum=0.0_DP
            !Dynamic matrices
            IF(ASSOCIATED(dynamicMapping)) THEN
              DO dynamicMatrixIdx=1,numberOfDynamicMatrices
                IF(dynamicMatrixIdx==stiffnessMatrixNumber) THEN
                  CALL DistributedMatrix_MatrixRowByVectorAdd(stiffnessDistributedMatrix,.FALSE.,displacementDistributedVector, &
                    & equationsRowNumber,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,stiffnessAlpha,lhsSum,err,error,*999)
                ELSE IF(dynamicMatrixIdx==dampingMatrixNumber) THEN
                  CALL DistributedMatrix_MatrixRowByVectorAdd(dampingDistributedMatrix,.FALSE.,velocityDistributedVector, &
                    & equationsRowNumber,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,dampingAlpha,lhsSum,err,error,*999)
                ELSE IF(dynamicMatrixIdx==massMatrixNumber) THEN
                  CALL DistributedMatrix_MatrixRowByVectorAdd(massDistributedMatrix,.FALSE.,accelerationDistributedVector, &
                    & equationsRowNumber,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,massAlpha,lhsSum,err,error,*999)
                ELSE
                  localError="The dynamic matrix number of "//TRIM(NumberToVString(dynamicMatrixIdx,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
              ENDDO !dynamicMatrixIdx
            ENDIF
            !Linear matrices
            IF(ASSOCIATED(linearMapping)) THEN
              DO linearMatrixIdx=1,numberOfLinearMatrices
                NULLIFY(linearVariable)
                CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,linearMatrixIdx,linearVariable,err,error,*999)
                !Get the dependent field variable parameters
                NULLIFY(linearParameters)
                CALL FieldVariable_ParameterSetGet(linearVariable,FIELD_VALUES_SET_TYPE,linearParameters,err,error,*999)
                NULLIFY(linearDistributedVector)
                CALL FieldParameterSet_ParametersGet(linearParameters,linearDistributedVector,err,error,*999)
                NULLIFY(linearMatrix)
                CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,linearMatrixIdx,linearMatrix,err,error,*999)
                NULLIFY(linearDistributedMatrix)
                CALL EquationsMatrix_DistributedMatrixGet(linearMatrix,linearDistributedMatrix,err,error,*999)
                CALL EquationsMatrix_MatrixCoefficientGet(linearMatrix,linearAlpha,err,error,*999)
                CALL DistributedMatrix_MatrixRowByVectorAdd(linearDistributedMatrix,.FALSE.,linearDistributedVector, &
                  & equationsRowNumber,DISTRIBUTED_MATRIX_VECTOR_NO_GHOSTS_TYPE,linearAlpha,lhsSum,err,error,*999)
              ENDDO !linearMatrixIdx
            ENDIF
            !Residual vectors
            IF(ASSOCIATED(nonlinearMapping)) THEN
              DO residualIdx=1,numberOfResiduals
                NULLIFY(residualVector)
                CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
                NULLIFY(residualDistributedVector)
                CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & residualDistributedVector,err,error,*999)
                CALL EquationsMatricesResidual_VectorCoefficientGet(residualVector,residualAlpha,err,error,*999)
                CALL DistributedVector_ValuesGet(residualDistributedVector,equationsRowNumber,residualValue,err,error,*999)
                lhsSum=lhsSum+residualAlpha*residualValue
              ENDDO !residualIdx
            ENDIF
            !Source vectors
            IF(ASSOCIATED(sourcesMapping)) THEN
              DO sourceIdx=1,numberOfSources
                NULLIFY(sourceVector)
                CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
                NULLIFY(sourceDistributedVector)
                CALL EquationsMatricesSource_DistributedVectorGet(sourceVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
                  & sourceDistributedVector,err,error,*999)
                CALL EquationsMatricesSource_VectorCoefficientGet(sourceVector,sourceAlpha,err,error,*999)
                CALL DistributedVector_ValuesGet(sourceDistributedVector,equationsRowNumber,sourceValue,err,error,*999)
                lhsSum=lhsSum+sourceAlpha*SourceValue              
              ENDDO !sourceIdx
            ENDIF
            !Set the RHS value
            CALL FieldVariable_ParameterSetUpdateLocalDOF(rhsVariable,FIELD_VALUES_SET_TYPE,rhsVariableDOF,lhsSum,err,error,*999)
          CASE(BOUNDARY_CONDITION_NEUMANN_ROW)
            !OK, do nothing
          CASE(BOUNDARY_CONDITION_ROBIN_ROW)
            !Robin boundary conditions
            CALL FlagError("Robin boundary conditions are not implemented.",err,error,*999)
          CASE(BOUNDARY_CONDITION_CAUCHY_ROW)
            !Cauchy boundary conditions
            CALL FlagError("Cauchy boundary conditions are not implemented.",err,error,*999)
          CASE(BOUNDARY_CONDITION_CONSTRAINED_ROW)
            !Constrained row boundary conditions
            CALL FlagError("Constrained row boundary conditions are not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The LHS boundary condition of "// &
              & TRIM(NumberToVString(lhsBoundaryCondition,"*",err,error))// &
              & " for RHS variable dof number "// &
              & TRIM(NumberToVString(lhsVariableDOF,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDDO !equationsRowNumber
        !Update RHS variable
        CALL FieldVariable_ParameterSetUpdateStart(rhsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(rhsVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDIF !Number of Dirichlet conditions > 0
    ENDIF !RHS mapping
          
    EXITS("EquationsSet_Backsubstitute")
    RETURN
999 ERRORSEXITS("EquationsSet_Backsubstitute",err,error)
    RETURN 1
   
  END SUBROUTINE EquationsSet_Backsubstitute
  
  !
  !================================================================================================================================
  !

  !>Set boundary conditions for an equation set according to the analytic equations. \see OpenCMISS::cmfe_EquationsSet_BoundaryConditionsAnalytic
  SUBROUTINE EquationsSet_BoundaryConditionsAnalytic(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the analyticboundary conditions for.
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_BoundaryConditionsAnalytic",err,error,*999)

    CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
    CALL EquationsSet_AssertAnalyticIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "//TRIM(NumberToVString(equationsSet%specification(1),"*", &
        & err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("EquationsSet_BoundaryConditionsAnalytic")
    RETURN
999 ERRORSEXITS("EquationsSet_BoundaryConditionsAnalytic",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_BoundaryConditionsAnalytic

  !
  !================================================================================================================================
  !

  !>Checks that the fields and region in an equations set are compitable
  SUBROUTINE EquationsSet_FieldRegionSetupCheck(equationsSet,checkString,fieldUserNumber,field,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to check the setup for.
    CHARACTER(LEN=*), INTENT(IN) :: checkString !<The string on the equations set part to check.
    INTEGER(INTG), INTENT(IN) :: fieldUserNumber !<The user number of the field to check
    TYPE(FieldType), POINTER :: field !<A pointer to the field to check the setup for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: actualFieldDecompositionUserNumber,actualFieldRegionUserNumber,actualFieldUserNumber, &
      & equationsSetRegionUserNumber,geometricDecompositionUserNumber
    TYPE(DecompositionType), POINTER :: fieldDecomposition,geometricDecomposition
    TYPE(FieldType), POINTER :: checkField,geometricField
    TYPE(RegionType), POINTER :: actualFieldRegion,equationsSetRegion
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_FieldRegionSetupCheck",err,error,*999)

    NULLIFY(equationsSetRegion)
    CALL EquationsSet_RegionGet(equationsSet,equationsSetRegion,err,error,*999)
    CALL Region_UserNumberGet(equationsSetRegion,equationsSetRegionUserNumber,err,error,*999)
    IF(ASSOCIATED(field)) THEN
      !Check that the specified field and the specified field user number matches
      CALL Field_AssertIsFinished(field,err,error,*999)
      CALL Field_UserNumberGet(field,actualFieldUserNumber,err,error,*999)
      IF(fieldUserNumber/=actualFieldUserNumber) THEN
        localError="The specified "//TRIM(checkString)//" field user number of "// &
          & TRIM(NumberToVString(fieldUserNumber,"*",err,error))//" does not match the user number of the specified "// &
          & TRIM(checkString)//" field of "//TRIM(NumberToVString(actualFieldUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check that the specified field has been created on the same region as the equations set
      NULLIFY(actualFieldRegion)
      CALL Field_RegionGet(field,actualFieldRegion,err,error,*999)
      CALL Region_UserNumberGet(actualFieldRegion,actualFieldRegionUserNumber,err,error,*999)
      IF(equationsSetRegionUserNumber/=actualFieldRegionUserNumber) THEN
        localError="Invalid region setup. The specified "//TRIM(checkString)//" field has been created on region number "// &
          & TRIM(NumberToVString(actualFieldRegionUserNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(equationsSetRegionUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      geometricField=>equationsSet%geometry%geometricField
      IF(ASSOCIATED(geometricField)) THEN
        !Check that the specified field has the same decomposition as the geometric field
        NULLIFY(geometricDecomposition)
        CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
        CALL Decomposition_UserNumberGet(geometricDecomposition,geometricDecompositionUserNumber,err,error,*999)
        NULLIFY(fieldDecomposition)
        CALL Field_DecompositionGet(field,fieldDecomposition,err,error,*999)
        CALL Decomposition_UserNumberGet(fieldDecomposition,actualFieldDecompositionUserNumber,err,error,*999)
        IF(.NOT.ASSOCIATED(geometricDecomposition,fieldDecomposition)) THEN
          localError="Invalid decomposition setup. The specified "//TRIM(checkString)// &
            & " field has been decomposed with decomposition number "// &
            & TRIM(NumberToVString(actualFieldDecompositionUserNumber,"*",err,error))// &
            & " and the equation set geometric field has been decomposed with decomposition number "// &
            & TRIM(NumberToVString(geometricDecompositionUserNumber,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(checkField)
      CALL Field_UserNumberFind(fieldUserNumber,equationsSetRegion,checkField,err,error,*999)
      IF(ASSOCIATED(checkField)) THEN
        localError="The specified "//TRIM(checkString)//" field user number of "// &
          & TRIM(NumberToVString(fieldUserNumber,"*",err,error))// &
          & " has already been used to create a field on region number "// &
          & TRIM(NumberToVString(equationsSetRegionUserNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF       
    ENDIF
    
    EXITS("EquationsSet_FieldRegionSetupCheck")
    RETURN
999 ERRORSEXITS("EquationsSet_FieldRegionSetupCheck",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FieldRegionSetupCheck

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equation set on a region. \see OpenCMISS::Iron::cmfe_EquationsSet_CreateStart
  SUBROUTINE EquationsSet_CreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    
    ENTERS("EquationsSet_CreateFinish",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)
    
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_INITIAL_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    !Finish the equations set specific setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    !Finish the equations set specific geometry setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish the equations set creation
    equationsSet%equationsSetFinished=.TRUE.
  
    EXITS("EquationsSet_CreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_CreateFinish",err,error)
    RETURN 1
   
  END SUBROUTINE EquationsSet_CreateFinish
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an equations set defined by userNumber in the region identified by REGION. \see OpenCMISS::Iron::cmfe_EquationsSet_CreateStart
  !>Default values set for the EQUATIONS_SET's attributes are:
  !>- LINEARITY: 1 (EQUATIONS_SET_LINEAR)
  !>- TIME_DEPENDENCE: 1 (EQUATIONS_SET_STATIC)
  !>- SOLUTION_METHOD: 1 (EQUATIONS_SET_FEM_SOLUTION_METHOD)
  !>- GEOMETRY 
  !>- MATERIALS 
  !>- SOURCE 
  !>- DEPENDENT
  !>- ANALYTIC
  !>- FIXED_CONDITIONS 
  !>- EQUATIONS 
  SUBROUTINE EquationsSet_CreateStart(userNumber,equationsSetRegion,geometricFibreField,equationsSetSpecification,&
      & equationsSetFieldUserNumber,equationsSetField,equationsSet,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: userNumber !<The user number of the equations set
    TYPE(RegionType), POINTER :: equationsSetRegion !<A pointer to the region to create the equations set on
    TYPE(FieldType), POINTER :: geometricFibreField !<A pointer to the either the geometry or, if appropriate, the fibre field for the equation set
    INTEGER(INTG), INTENT(IN) :: equationsSetSpecification(:) !<The equations set specification array to set
    INTEGER(INTG), INTENT(IN) :: equationsSetFieldUserNumber !<The user number of the equations set field
    TYPE(FieldType), POINTER :: equationsSetField !<On return, a pointer to the equations set field
    TYPE(EquationsSetType), POINTER :: equationsSet !<On return, a pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equationsSetIdx,equationsSetRegionUserNumber,geometricFibreFieldRegionUserNumber, &
      & geometricFibreFieldType
    TYPE(DecompositionType), POINTER :: equationsSetFieldDecomposition,geometricFibreFieldDecomposition
    TYPE(EquationsSetType), POINTER :: newEquationsSet
    TYPE(EquationsSetPtrType), POINTER :: newEquationsSets(:)
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field
    TYPE(RegionType), POINTER :: geometricFibreFieldRegion,equationsSetFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_CreateStart",err,error,*997)

    CALL Region_UserNumberGet(equationsSetRegion,equationsSetRegionUserNumber,err,error,*999)
    CALL Field_AssertIsFinished(geometricFibreField,err,error,*999)
    NULLIFY(newEquationsSet)
    CALL EquationsSet_UserNumberFind(userNumber,equationsSetRegion,newEquationsSet,err,error,*997)
    IF(ASSOCIATED(newEquationsSet)) THEN
      localError="Equations set user number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(equationsSetRegionUserNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*997)
    ENDIF
    CALL Field_TypeGet(geometricFibreField,geometricFibreFieldType,err,error,*999)
    IF(geometricFibreFieldType/=FIELD_GEOMETRIC_TYPE.AND.geometricFibreFieldType/=FIELD_FIBRE_TYPE) &
      & CALL FlagError("The specified geometric field is not a geometric or fibre field.",err,error,*997)
    NULLIFY(geometricFibreFieldRegion)
    CALL Field_RegionGet(geometricFibreField,geometricFibreFieldRegion,err,error,*999)
    CALL Region_UserNumberGet(geometricFibreFieldRegion,geometricFibreFieldRegionUserNumber,err,error,*999)
    IF(geometricFibreFieldRegionUserNumber/=equationsSetRegionUserNumber) THEN
      localError="The geometric field region and the specified region do not match. "// &
        & "The geometric field was created on region number "// &
        & TRIM(NumberToVString(geometricFibreFieldRegionUserNumber,"*",err,error))// &
        & " and the specified equation set region number is "// &
        & TRIM(NumberToVString(equationsSetRegionUserNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equationsSetField)) THEN
      !Check the specified equations set field has the same decomposition as the geometric field
      NULLIFY(geometricFibreFieldDecomposition)
      CALL Field_DecompositionGet(geometricFibreField,geometricFibreFieldDecomposition,err,error,*999)
      NULLIFY(equationsSetFieldDecomposition)
      CALL Field_DecompositionGet(equationsSetField,equationsSetFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(geometricFibreFieldDecomposition,equationsSetFieldDecomposition)) THEN
        CALL FlagError("The specified equations set field does not have the same decomposition "// &
          & "as the geometric field for the specified equations set.",err,error,*999)
      ENDIF
    ENDIF
    
    !Initalise equations set
    NULLIFY(newEquationsSet)
    CALL EquationsSet_Initialise(newEquationsSet,err,error,*999)
    !Set default equations set values
    newEquationsSet%userNumber=userNumber
    newEquationsSet%globalNumber=equationsSetRegion%equationsSets%numberOfEquationsSets+1
    newEquationsSet%equationsSets=>equationsSetRegion%equationsSets
    newEquationsSet%label="Equations Set "//TRIM(NumberToVString(userNumber,"*",err,error))
    newEquationsSet%region=>equationsSetRegion
    !Check field
    CALL EquationsSet_FieldRegionSetupCheck(newEquationsSet,"equations set",equationsSetFieldUserNumber,equationsSetField, &
      & err,error,*999)
    !Set the equations set class, type and subtype
    CALL EquationsSet_SpecificationSet(newEquationsSet,equationsSetSpecification,err,error,*999)
    newEquationsSet%equationsSetFinished=.FALSE.
    !Initialise the setup
    newEquationsSet%equationsField%equationsSetFieldAutoCreated=(.NOT.ASSOCIATED(equationsSetField))
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_INITIAL_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    !Here, we get a pointer to the equations_set_field; default is null
    equationsSetSetupInfo%fieldUserNumber=equationsSetFieldUserNumber
    equationsSetSetupInfo%field=>equationsSetField
    !Start equations set specific setup
    CALL EquationsSet_Setup(newEquationsSet,equationsSetSetupInfo,err,error,*999)
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set up the equations set geometric fields
    CALL EquationsSet_GeometryInitialise(newEquationsSet,err,error,*999)
    IF(geometricFibreFieldType==FIELD_GEOMETRIC_TYPE) THEN
      newEquationsSet%geometry%geometricField=>geometricFibreField
      NULLIFY(newEquationsSet%geometry%fibreField)
    ELSE
      newEquationsSet%geometry%geometricField=>geometricFibreField%geometricField
      newEquationsSet%geometry%fibreField=>geometricFibreField
    ENDIF
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=geometricFibreField%userNumber
    equationsSetSetupInfo%field=>geometricFibreField
    !Set up equations set specific geometry
    CALL EquationsSet_Setup(newEquationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Add new equations set into list of equations set in the region
    ALLOCATE(newEquationsSets(equationsSetRegion%equationsSets%numberOfEquationsSets+1),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
    DO equationsSetIdx=1,equationsSetRegion%equationsSets%numberOfEquationsSets
      newEquationsSets(equationsSetIdx)%ptr=>equationsSetRegion%equationsSets%equationsSets(equationsSetIdx)%ptr
    ENDDO !equationsSetIdx
    newEquationsSets(equationsSetRegion%equationsSets%numberOfEquationsSets+1)%ptr=>newEquationsSet
    IF(ASSOCIATED(equationsSetRegion%equationsSets%equationsSets)) DEALLOCATE(equationsSetRegion%equationsSets%equationsSets)
    equationsSetRegion%equationsSets%equationsSets=>newEquationsSets
    equationsSetRegion%equationsSets%numberOfEquationsSets=equationsSetRegion%equationsSets%numberOfEquationsSets+1
    equationsSet=>newEquationsSet
    !\todo check pointer setup
    IF(equationsSet%equationsField%equationsSetFieldAutoCreated) THEN
      equationsSetField=>equationsSet%equationsField%equationsSetField
    ELSE
      equationsSet%equationsField%equationsSetField=>equationsSetField
    ENDIF
   
    EXITS("EquationsSet_CreateStart")
    RETURN
999 IF(ASSOCIATED(newEquationsSet))CALL EquationsSet_Finalise(newEquationsSet,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(newEquationsSets)) DEALLOCATE(newEquationsSets)
997 ERRORSEXITS("EquationsSet_CreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_CreateStart
  
  !
  !================================================================================================================================
  !

  !>Destroys an equations set identified by a pointer and deallocates all memory. \see OpenCMISS::Iron::cmfe_EquationsSet_Destroy
  SUBROUTINE EquationsSet_Destroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,equationsSetPosition
    TYPE(EquationsSetsType), POINTER :: equationsSets
    TYPE(EquationsSetPtrType), POINTER :: newEquationsSets(:)

    NULLIFY(newEquationsSets)

    ENTERS("EquationsSet_Destroy",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*998)
    NULLIFY(equationsSets)
    CALL EquationsSet_EquationsSetsGet(equationsSet,equationsSets,err,error,*999)
    equationsSetPosition=equationsSet%globalNumber

    !Destroy all the equations set components
    CALL EquationsSet_Finalise(equationsSet,err,error,*999)
    
    !Remove the equations set from the list of equations set
    IF(equationsSets%numberOfEquationsSets>1) THEN
      ALLOCATE(newEquationsSets(equationsSets%numberOfEquationsSets-1),STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
      DO equationsSetIdx=1,equationsSets%numberOfEquationsSets
        IF(equationsSetIdx<equationsSetPosition) THEN
          newEquationsSets(equationsSetIdx)%ptr=>equationsSets%equationsSets(equationsSetIdx)%ptr
        ELSE IF(equationsSetIdx>equationsSetPosition) THEN
          equationsSets%equationsSets(equationsSetIdx)%ptr%globalNumber=equationsSets% &
            & equationsSets(equationsSetIdx)%ptr%globalNumber-1
          newEquationsSets(equationsSetIdx-1)%ptr=>equationsSets%equationsSets(equationsSetIdx)%ptr
        ENDIF
      ENDDO !equationsSetIdx
      IF(ASSOCIATED(equationsSets%equationsSets)) DEALLOCATE(equationsSets%equationsSets)
      equationsSets%equationsSets=>newEquationsSets
      equationsSets%numberOfEquationsSets=equationsSets%numberOfEquationsSets-1
    ELSE
      DEALLOCATE(equationsSets%equationsSets)
      equationsSets%numberOfEquationsSets=0
    ENDIF

    EXITS("EquationsSet_Destroy")
    RETURN
999 IF(ASSOCIATED(newEquationsSets)) DEALLOCATE(newEquationsSets)
998 ERRORSEXITS("EquationsSet_Destroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_Destroy
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations set and deallocate all memory.
  SUBROUTINE EquationsSet_Finalise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_Finalise",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      CALL EquationsSet_GeometryFinalise(equationsSet%geometry,err,error,*999)
      CALL EquationsSet_DependentFinalise(equationsSet%dependent,err,error,*999)
      CALL EquationsSet_IndependentFinalise(equationsSet%independent,err,error,*999)
      CALL EquationsSet_MaterialsFinalise(equationsSet%materials,err,error,*999)
      CALL EquationsSet_SourceFinalise(equationsSet%source,err,error,*999)
      CALL EquationsSet_AnalyticFinalise(equationsSet%analytic,err,error,*999)
      CALL EquationsSet_EquationsFieldFinalise(equationsSet%equationsField,err,error,*999)
      CALL EquationsSet_DerivedFinalise(equationsSet%derived,err,error,*999)
      IF(ASSOCIATED(equationsSet%equations)) CALL Equations_Destroy(equationsSet%equations,err,error,*999)
      IF(ALLOCATED(equationsSet%specification)) DEALLOCATE(equationsSet%specification)
      equationsSet%label=""
      DEALLOCATE(equationsSet)
    ENDIF
       
    EXITS("EquationsSet_Finalise")
    RETURN
999 ERRORSEXITS("EquationsSet_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_Finalise

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calcualte
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,numberOfDynamicMatrices,numberOfLinearMatrices,numberOfSources,outputType,sourceIdx
    LOGICAL :: updateMatrix,updateVector
    TYPE(ElementMatrixType), POINTER :: elementMatrix
    TYPE(ElementVectorType), POINTER :: elementVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_FiniteElementCalculate",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)

    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL Fitting_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL Bioelectric_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",err,error,*999)
        CALL EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",numberOfDynamicMatrices,err,error,*999)
        DO matrixIdx=1,numberOfDynamicMatrices
          NULLIFY(dynamicMatrix)
          CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(dynamicMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,dynamicMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",numberOfLinearMatrices,err,error,*999)
        DO matrixIdx=1,numberOfLinearMatrices
          NULLIFY(linearMatrix)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,linearMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(sourceVectors)
      CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
      IF(ASSOCIATED(sourceVectors)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element source vectors:",err,error,*999)
        CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of source vectors = ",numberOfSources,err,error,*999)
        DO sourceIdx=1,numberOfSources
          NULLIFY(sourceVector)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Source vector : ",sourceIdx,err,error,*999)
          CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
          IF(updateVector) CALL EquationsMatricesSource_ElementVectorOutput(GENERAL_OUTPUT_TYPE,sourceVector,err,error,*999)
        ENDDO !sourceIdx
      ENDIF
      NULLIFY(rhsVector)
      CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element RHS vector :",err,error,*999)
        CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
        IF(updateVector) CALL EquationsMatricesRHS_ElementVectorOutput(GENERAL_OUTPUT_TYPE,rhsVector,err,error,*999)
      ENDIF
    ENDIF
       
    EXITS("EquationsSet_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("EquationsSet_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian for the given element number for a finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianCalculationType,matrixIdx,numberOfJacobians,numberOfResiduals,outputType,residualIdx
    LOGICAL :: updateJacobian
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_FiniteElementJacobianEvaluate",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    CALL EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*999)
    
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    DO residualIdx=1,numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
      CALL EquationsMatricesResidual_NumberOfJacobiansGet(residualvector,numberOfJacobians,err,error,*999)
      DO matrixIdx=1,numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
        CALL JacobianMatrix_CalculationTypeGet(jacobianMatrix,jacobianCalculationType,err,error,*999)
        SELECT CASE(jacobianCalculationType)
        CASE(EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
          ! None of these routines currently support calculating off diagonal terms for coupled problems,
          ! but when one does we will have to pass through the matrixIdx parameter
!!TODO: need to pass list of residuals and jacobians to calculate to the FE routines otherwise they will compute all the
!!      jacobians multiple times inside these loops.
          IF(matrixIdx>1) CALL FlagError("Analytic off-diagonal Jacobian calculation not implemented.",err,error,*999)
          SELECT CASE(equationsSet%specification(1))
          CASE(EQUATIONS_SET_ELASTICITY_CLASS)
            CALL Elasticity_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
            CALL FluidMechanics_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
          CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
            CALL ClassicalField_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
          CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MODAL_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
            CALL MultiPhysics_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*999)
          CASE DEFAULT
            localError="The first equations set specification of"// &
              & TRIM(NumberToVString(equationsSet%specification(1),"*", &
              & err,error))//" is not valid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED)
          CALL EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,residualIdx,matrixIdx,elementNumber,err,error,*999)
        CASE DEFAULT
          localError="The Jacobian calculation type of "//TRIM(NumberToVString(jacobianCalculationType,"*",err,error))// &
            & " is not valid for Jacobian matrix index "//TRIM(NumberToVString(matrixIdx,"*",err,error))// &
            & " of residual index "//TRIM(NumberToVString(residualIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !matrixIdx
    ENDDO !residualIdx
    IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element Jacobian matrices:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of residuals = ",numberOfResiduals,err,error,*999)
      DO residualIdx=1,numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Residual number : ",residualIdx,err,error,*999)
        CALL EquationsMatricesResidual_NumberOfJacobiansGet(residualvector,numberOfJacobians,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of Jacobians = ",numberOfJacobians,err,error,*999)             
        DO matrixIdx=1,numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Jacobian matrix : ",matrixIdx,err,error,*999)
          CALL JacobianMatrix_UpdateMatrixGet(jacobianMatrix,updateJacobian,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",updateJacobian,err,error,*999)
          IF(updateJacobian) CALL JacobianMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,jacobianMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF
      
    EXITS("EquationsSet_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_FiniteElementJacobianEvaluate",err,error)
    RETURN 1

  END SUBROUTINE EquationsSet_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the element Jacobian matrix entries using finite differencing for a general finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,residualNumber,jacobianNumber,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet  !<A pointer to the equations set to evaluate the element Jacobian for
    INTEGER(INTG), INTENT(IN) :: residualNumber  !<The residual number to calculate the Jacobian for
    INTEGER(INTG), INTENT(IN) :: jacobianNumber  !<The Jacobian number to calculate the Jacobian for
    INTEGER(INTG), INTENT(IN) :: elementNumber  !<The element number to calculate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err  !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    INTEGER(INTG) :: column,componentIdx,componentInterpolationType,derivative,derivativeIdx,globalDerivativeIndex,localDOF, &
      & node,nodeIdx,numberOfComponents,numberOfLocalNodes,numberOfNodeDerivatives,numberOfResiduals,numberOfRows,residualIdx, &
      & version
    REAL(DP) :: delta,origDepVar
    LOGICAL :: updateMatrix,updateVector
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DistributedVectorType), POINTER :: parameters
    TYPE(ElementVectorType) :: elementVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldParameterSetType), POINTER :: parameterSet
    TYPE(FieldVariableType), POINTER :: rowVariable,columnVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix

    ENTERS("EquationsSet_FiniteElementJacobianEvaluateFD",err,error,*999)

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
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualNumber,residualVector,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)    
    NULLIFY(rowVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowVariable,err,error,*999)
    ! For coupled problems this routine will be called multiple times if multiple Jacobians use finite
    ! differencing, so make sure we only calculate the residual vector once, to save time and because
    ! it would otherwise add together. TODO: sort a better way to do this.
    IF(residualVector%elementResidualCalculated/=elementNumber) &
      & CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    ! Make a temporary copy of the unperturbed residuals
    elementVector=residualVector%elementResidual
    ! For coupled nonlinear problems there will be multiple Jacobians
    ! For this equations set, we calculate the residual for the row variable
    ! while pertubing parameters from the column variable.
    ! For non coupled problems these two variables will be the same
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualNumber,residualMapping,err,error,*999)
    NULLIFY(columnVariable)
    CALL EquationsMappingResidual_VariableGet(residualMapping,jacobianNumber,columnVariable,err,error,*999)
    NULLIFY(parameterSet)
    CALL FieldVariable_ParameterSetGet(columnVariable,FIELD_VALUES_SET_TYPE,parameterSet,err,error,*999)
    NULLIFY(parameters)
    CALL FieldParameterSet_ParametersGet(parameterSet,parameters,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualNumber,residualVector,err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,jacobianNumber,jacobianMatrix,err,error,*999)
    numberOfRows=jacobianMatrix%elementJacobian%numberOfRows
    IF(numberOfRows/=residualVector%elementResidual%numberOfRows) &
      & CALL FlagError("Element matrix number of rows does not match element residual vector size.",err,error,*999)
    ! determine step size
    CALL DistributedVector_L2Norm(parameters,delta,err,error,*999)
    !delta=(1.0_DP+delta)*1.0E-6_DP
    delta=(1.0_DP+delta)*jacobianMatrix%jacobianFiniteDifferenceStepSize
    ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
    ! distributed out, have to use proper field accessing routines..
    ! so let's just loop over component, node/el, derivative
    column=0  ! element jacobian matrix column number
    CALL FieldVariable_NumberOfComponentsGet(columnVariable,numberOfComponents,err,error,*999)
    DO componentIdx=1,numberOfComponents
      NULLIFY(domain)
      CALL FieldVariable_ComponentDomainGet(columnVariable,componentIdx,domain,err,error,*999)
      NULLIFY(domainTopology)
      CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
      NULLIFY(domainElements)
      CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
      CALL FieldVariable_ComponentInterpolationGet(columnVariable,componentIdx,componentInterpolationType,err,error,*999)
      SELECT CASE(componentInterpolationType)
      CASE (FIELD_NODE_BASED_INTERPOLATION)
        NULLIFY(basis)
        CALL DomainElements_ElementBasisGet(domainElements,elementNumber,basis,err,error,*999)
        CALL Basis_NumberOfLocalNodesGet(basis,numberOfLocalNodes,err,error,*999)
        DO nodeIdx=1,numberOfLocalNodes
          CALL DomainElements_ElementNodeGet(domainElements,nodeIdx,elementNumber,node,err,error,*999)
          CALL Basis_NodeNumberOfDerivativesGet(basis,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainElements_ElementDerivativeGet(domainElements,derivativeIdx,nodeIdx,elementNumber,globalDerivativeIndex, &
              & err,error,*999)
            CALL DomainElements_ElementVersionGet(domainElements,derivativeIdx,nodeIdx,elementNumber,version,err,error,*999)
            CALL FieldVariable_LocalNodeDOFGet(columnVariable,version,globalDerivativeIndex,node,componentIdx,localDOF, &
              & err,error,*999)
            ! one-sided finite difference
            CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
            CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
            residualVector%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
            CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
            CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
            column=column+1
            jacobianMatrix%elementJacobian%matrix(1:numberOfRows,column)= &
              & (residualVector%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
        CALL FieldVariable_LocalElementDOFGet(columnVariable,elementNumber,componentIdx,localDOF,err,error,*999)
        ! one-sided finite difference
        CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
        CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
        residualVector%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
        CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
        CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
        column=column+1
        jacobianMatrix%elementJacobian%matrix(1:numberOfRows,column)= &
          & (residualVector%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
      CASE DEFAULT
        CALL FlagError("Unsupported type of interpolation.",err,error,*999)
      END SELECT
    END DO !componentIdx
    ! put the original residual back in
    residualVector%elementResidual=elementVector

    EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN
999 ERRORS("EquationsSet_FiniteElementJacobianEvaluateFD",err,error)
    EXITS("EquationsSet_FiniteElementJacobianEvaluateFD")
    RETURN 1
    
  END SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD

  !
  !================================================================================================================================
  !

  !>Evaluates the element residual and rhs vector for the given element number for a finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,numberOfDynamicMatrices,numberOfLinearMatrices,numberOfResiduals,numberOfSources,outputType, &
      & residualIdx,sourceIdx
    LOGICAL :: updateMatrix,updateVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_FiniteElementResidualEvaluate",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)

    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element matrices and vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      NULLIFY(dynamicMatrices)
      CALL EquationsMatricesVector_DynamicMatricesExists(vectorMatrices,dynamicMatrices,err,error,*999)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",err,error,*999)
        CALL EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of dynamic matrices = ",numberOfDynamicMatrices,err,error,*999)
        DO matrixIdx=1,numberOfDynamicMatrices
          NULLIFY(dynamicMatrix)
          CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Dynamic matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(dynamicMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,dynamicMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of linear matrices = ",numberOfLinearMatrices,err,error,*999)
        DO matrixIdx=1,numberOfLinearMatrices
          NULLIFY(linearMatrix)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Linear matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,linearMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Residual vectors:",err,error,*999)
        CALL EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of residuals = ",numberOfResiduals,err,error,*999)
        DO residualIdx=1,numberOfResiduals
          NULLIFY(residualVector)
          CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Residual vector : ",residualIdx,err,error,*999)
          CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
          IF(updateVector) CALL EquationsMatricesResidual_ElementVectorOutput(GENERAL_OUTPUT_TYPE,residualVector,err,error,*999)
        ENDDO !sourceIdx
      ENDIF
      NULLIFY(sourceVectors)
      CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
      IF(ASSOCIATED(sourceVectors)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Source vectors:",err,error,*999)
        CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of source vectors = ",numberOfSources,err,error,*999)
        DO sourceIdx=1,numberOfSources
          NULLIFY(sourceVector)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Source vector : ",sourceIdx,err,error,*999)
          CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
          IF(updateVector) CALL EquationsMatricesSource_ElementVectorOutput(GENERAL_OUTPUT_TYPE,sourceVector,err,error,*999)
        ENDDO !sourceIdx
      ENDIF
      NULLIFY(rhsVector)
      CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"RHS vector :",err,error,*999)
        CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
        IF(updateVector) CALL EquationsMatricesRHS_ElementVectorOutput(GENERAL_OUTPUT_TYPE,rhsVector,err,error,*999)
      ENDIF
    ENDIF
       
    EXITS("EquationsSet_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_FiniteElementResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Finish the creation of independent variables for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_IndependentCreateFinish
  SUBROUTINE EquationsSet_IndependentCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish the creation of the independent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: independentField

    ENTERS("EquationsSet_IndependentCreateFinish",err,error,*999)

    CALL EquationsSet_AssertIndependentNotFinished(equationsSet,err,error,*999)
    
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    NULLIFY(independentField)
    CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
    independentField=>equationsSet%INDEPENDENT%independentField
    equationsSetSetupInfo%fieldUserNumber=independentField%userNumber
    equationsSetSetupInfo%field=>independentField
    !Finish equations set specific startup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish independent creation
    equationsSet%independent%independentFinished=.TRUE.
       
    EXITS("EquationsSet_IndependentCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_IndependentCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of independent variables for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_IndependentCreateStart
  SUBROUTINE EquationSet_IndependentCreateStart(equationsSet,independentFieldUserNumber,independentField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: independentFieldUserNumber !<The user specified independent field number
    TYPE(FieldType), POINTER :: independentField !<If associated on entry, a pointer to the user created independent field which has the same user number as the specified independent field user number. If not associated on entry, on exit, a pointer to the created independent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: geometricDecomposition,independentDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,independentFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationSet_IndependentCreateStart",err,error,*998)

    CALL EquationsSet_AssertIndependentNotCreated(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"independent",independentFieldUserNumber,independentField,err,error,*998)
    
    !Initialise the equations set independent
    CALL EquationsSet_IndependentInitialise(equationsSet,err,error,*999)
    equationsSet%independent%independentFieldAutoCreated=(.NOT.ASSOCIATED(independentField))
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=independentFieldUserNumber
    equationsSetSetupInfo%field=>independentField
    !Start equations set specific startup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(equationsSet%independent%independentFieldAutoCreated) THEN            
      independentField=>equationsSet%independent%independentField
    ELSE
      equationsSet%independent%independentField=>independentField
    ENDIF
       
    EXITS("EquationSet_IndependentCreateStart")
    RETURN
999 CALL EquationsSet_IndependentFinalise(equationsSet%independent,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationSet_IndependentCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationSet_IndependentCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the independent field for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_IndependentDestroy
  SUBROUTINE EquationsSet_IndependentDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy the independent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_IndependentDestroy",err,error,*999)

    CALL EquationsSet_AssertIndependentIsCreated(equationsSet,err,error,*999)
    
    CALL EquationsSet_IndependentFinalise(equationsSet%independent,err,error,*999)
       
    EXITS("EquationsSet_IndependentDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_IndependentDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentDestroy

  !
  !================================================================================================================================
  !

  !>Finalise the independent field for an equations set.
  SUBROUTINE EquationsSet_IndependentFinalise(equationsSetIndependent,err,error,*)

    !Argument variables
    TYPE(EquationsSetIndependentType), POINTER :: equationsSetIndependent !<A pointer to the equations set independent to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_IndependentFinalise",err,error,*999)

    IF(ASSOCIATED(equationsSetIndependent)) DEALLOCATE(equationsSetIndependent)
       
    EXITS("EquationsSet_IndependentFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_IndependentFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the independent field for an equations set.
  SUBROUTINE EquationsSet_IndependentInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the independent for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsSet_IndependentInitialise",err,error,*998)

    CALL EquationsSet_AssertAnalyticNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%independent,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set independent field.",err,error,*999)
    equationsSet%independent%equationsSet=>equationsSet
    equationsSet%independent%independentFinished=.FALSE.
    equationsSet%independent%independentFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%independent%independentField)
      
    EXITS("EquationsSet_IndependentInitialise")
    RETURN
999 CALL EquationsSet_IndependentFinalise(equationsSet%independent,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_IndependentInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_IndependentInitialise

  !
  !================================================================================================================================
  !

  !>Initialises an equations set.
  SUBROUTINE EquationsSet_Initialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The pointer to the equations set to initialise. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("EquationsSet_Initialise",err,error,*998)

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is already associated.",err,error,*998)
    
    ALLOCATE(equationsSet,STAT=err)    
    IF(err/=0) CALL FlagError("Could not allocate equations set.",err,error,*999)
    equationsSet%userNumber=0
    equationsSet%globalNumber=0
    equationsSet%equationsSetFinished=.FALSE.
    NULLIFY(equationsSet%equationsSets)
    equationsSet%label=""
    NULLIFY(equationsSet%region)
    equationsSet%specificationLength=0
    equationsSet%currentTime=0.0_DP
    equationsSet%deltaTime=0.0_DP
    equationsSet%outputType=EQUATIONS_SET_NO_OUTPUT
    equationsSet%solutionMethod=0
    CALL EquationsSet_GeometryInitialise(equationsSet,err,error,*999)
    CALL EquationsSet_DependentInitialise(equationsSet,err,error,*999)
    CALL EquationsSet_EquationsFieldInitialise(equationsSet,err,error,*999)
    NULLIFY(equationsSet%independent)
    NULLIFY(equationsSet%materials)
    NULLIFY(equationsSet%source)
    NULLIFY(equationsSet%analytic)
    NULLIFY(equationsSet%derived)
    NULLIFY(equationsSet%equations)
    NULLIFY(equationsSet%boundaryConditions)
       
    EXITS("EquationsSet_Initialise")
    RETURN
999 CALL EquationsSet_Finalise(equationsSet,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_Initialise

  !
  !================================================================================================================================
  !

  !>Finalise the geometry for an equations set
  SUBROUTINE EquationsSet_GeometryFinalise(equationsSetGeometry,err,error,*)

    !Argument variables
    TYPE(EquationsSetGeometryType) :: equationsSetGeometry !<A pointer to the equations set geometry to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_GeometryFinalise",err,error,*999)
    
    NULLIFY(equationsSetGeometry%geometricField)
    NULLIFY(equationsSetGeometry%fibreField)
       
    EXITS("EquationsSet_GeometryFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_GeometryFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_GeometryFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the geometry for an equation set
  SUBROUTINE EquationsSet_GeometryInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the geometry for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_GeometryInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    equationsSet%geometry%equationsSet=>equationsSet
    NULLIFY(equationsSet%geometry%geometricField)
    NULLIFY(equationsSet%geometry%fibreField)
        
    EXITS("EquationsSet_GeometryInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_GeometryInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_GeometryInitialise
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of materials for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_MaterialsCreateFinish
  SUBROUTINE EquationsSet_MaterialsCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish the creation of the materials field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: materialsField

    ENTERS("EquationsSet_MaterialsCreateFinish",err,error,*999)

    CALL EquationsSet_AssertMaterialsNotFinished(equationsSet,err,error,*999)
    
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_MATERIALS_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    materialsField=>equationsSet%materials%materialsField
    IF(ASSOCIATED(materialsField)) THEN
      equationsSetSetupInfo%fieldUserNumber=materialsField%userNumber
      equationsSetSetupInfo%field=>materialsField
    ENDIF
    !Finish equations set specific startup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish materials creation
    equationsSet%materials%materialsFinished=.TRUE.
       
    EXITS("EquationsSet_MaterialsCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_MaterialsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of materials for a problem. \see OpenCMISS::Iron::cmfe_EquationsSet_MaterialsCreateStart
  SUBROUTINE EquationsSet_MaterialsCreateStart(equationsSet,materialsFieldUserNumber,materialsField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: materialsFieldUserNumber !<The user specified materials field number
    TYPE(FieldType), POINTER :: materialsField !<If associated on entry, a pointer to the user created materials field which has the same user number as the specified materials field user number. If not associated on entry, on exit, a pointer to the created materials field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: geometricDecomposition,materialsDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,materialsFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_MaterialsCreateStart",err,error,*998)

    CALL EquationsSet_AssertMaterialsNotCreated(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"materials",materialsFieldUserNumber,materialsField,err,error,*998)
    
    !Initialise the equations set materials
    CALL EquationsSet_MaterialsInitialise(equationsSet,err,error,*999)
    equationsSet%materials%materialsFieldAutoCreated=(.NOT.ASSOCIATED(materialsField))
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_MATERIALS_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=materialsFieldUserNumber
    equationsSetSetupInfo%field=>materialsField
    !Start equations set specific startup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(equationsSet%materials%materialsFieldAutoCreated) THEN            
      materialsField=>equationsSet%materials%materialsField
    ELSE
      equationsSet%materials%materialsField=>materialsField
    ENDIF
       
    EXITS("EquationsSet_MaterialsCreateStart")
    RETURN
999 CALL EquationsSet_MaterialsFinalise(equationsSet%MATERIALS,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_MaterialsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the materials for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_MaterialsDestroy
  SUBROUTINE EquationsSet_MaterialsDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy the materials for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_MaterialsDestroy",err,error,*999)

    CALL EquationsSet_AssertMaterialsIsCreated(equationsSet,err,error,*999)
    
    CALL EquationsSet_MaterialsFinalise(equationsSet%materials,err,error,*999)
       
    EXITS("EquationsSet_MaterialsDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_MaterialsDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsDestroy

  !
  !================================================================================================================================
  !

  !>Finalise the materials for an equations set.
  SUBROUTINE EquationsSet_MaterialsFinalise(equationsSetMaterials,err,error,*)

    !Argument variables
    TYPE(EquationsSetMaterialsType), POINTER :: equationsSetMaterials !<A pointer to the equations set materials to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_MaterialsFinalise",err,error,*999)

    IF(ASSOCIATED(equationsSetMaterials)) DEALLOCATE(equationsSetMaterials)
       
    EXITS("EquationsSet_MaterialsFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_MaterialsFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the materials for an equations set.
  SUBROUTINE EquationsSet_MaterialsInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the materials for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsSet_MaterialsInitialise",err,error,*998)

    CALL EquationsSet_AssertMaterialsNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%materials,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set materials.",err,error,*999)
    equationsSet%materials%equationsSet=>equationsSet
    equationsSet%materials%materialsFinished=.FALSE.
    equationsSet%materials%materialsFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%materials%materialsField)
       
    EXITS("EquationsSet_MaterialsInitialise")
    RETURN
999 CALL EquationsSet_MaterialsFinalise(equationsSet%materials,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_MaterialsInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_MaterialsInitialise

  !
  !================================================================================================================================
  !

  !>Finish the creation of a dependent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_DependentCreateFinish
  SUBROUTINE EquationsSet_DependentCreateFinish(equationsSet,err,error,*)
    
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: dependentField

    ENTERS("EquationsSet_DependentCreateFinish",err,error,*999)

    CALL EquationsSet_AssertDependentNotFinished(equationsSet,err,error,*999)
    
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    equationsSetSetupInfo%fieldUserNumber=dependentField%userNumber
    equationsSetSetupInfo%field=>dependentField
    !Finish equations set specific setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish the equations set creation
    equationsSet%dependent%dependentFinished=.TRUE.
       
    EXITS("EquationsSet_DependentCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_DependentCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of dependent variables for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DependentCreateStart
  SUBROUTINE EquationsSet_DependentCreateStart(equationsSet,dependentFieldUserNumber,dependentField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of a dependent field on
    INTEGER(INTG), INTENT(IN) :: dependentFieldUserNumber !<The user specified dependent field number
    TYPE(FieldType), POINTER :: dependentField !<If associated on entry, a pointer to the user created dependent field which has the same user number as the specified dependent field user number. If not associated on entry, on exit, a pointer to the created dependent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,dependentFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EquationsSet_DependentCreateStart",err,error,*998)

    CALL EquationsSet_AssertDependentNotFinished(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"dependent",dependentFieldUserNumber,dependentField,err,error,*998)
    
    !Initialise the setup
    equationsSet%dependent%dependentFieldAutoCreated=(.NOT.ASSOCIATED(dependentField))
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=dependentFieldUserNumber
    equationsSetSetupInfo%field=>dependentField
    !Start the equations set specfic solution setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
      dependentField=>equationsSet%dependent%dependentField
    ELSE
      equationsSet%dependent%dependentField=>dependentField
    ENDIF
       
    EXITS("EquationsSet_DependentCreateStart")
    RETURN
999 CALL EquationsSet_DependentFinalise(equationsSet%dependent,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_DependentCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentCreateStart

  !
  !================================================================================================================================
  !
  
  !>Destroy the dependent variables for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DependentDestroy
  SUBROUTINE EquationsSet_DependentDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_DependentDestroy",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated",err,error,*999)

    CALL EquationsSet_DependentFinalise(equationsSet%dependent,err,error,*999)
    
    EXITS("EquationsSet_DependentDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_DependentDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentDestroy
  
  !
  !================================================================================================================================
  !

  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE EquationsSet_DependentFinalise(equationsSetDependent,err,error,*)

    !Argument variables
    TYPE(EquationsSetDependentType) :: equationsSetDependent !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_DependentFinalise",err,error,*999)

    NULLIFY(equationsSetDependent%equationsSet)
    equationsSetDependent%dependentFinished=.FALSE.
    equationsSetDependent%dependentFieldAutoCreated=.FALSE.
    NULLIFY(equationsSetDependent%dependentField)
    
    EXITS("EquationsSet_DependentFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_DependentFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialises the dependent variables for a equations set.
  SUBROUTINE EquationsSet_DependentInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_DependentInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    equationsSet%dependent%equationsSet=>equationsSet
    equationsSet%dependent%dependentFinished=.FALSE.
    equationsSet%dependent%dependentFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%dependent%dependentField)
       
    EXITS("EquationsSet_DependentInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_DependentInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DependentInitialise

  !
  !================================================================================================================================
  !

  !>Finish the creation of a derived variables field for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedCreateFinish
  SUBROUTINE EquationsSet_DerivedCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish the derived variable creation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) :: equationsSetSetupInfo
    TYPE(FieldType), POINTER :: derivedField

    ENTERS("EquationsSet_DerivedCreateFinish",err,error,*999)

    CALL EquationsSet_AssertDerivedNotFinished(equationsSet,err,error,*999)

    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_DERIVED_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    NULLIFY(derivedField)
    CALL EquationsSet_DerivedFieldGet(equationsSet,derivedField,err,error,*999)
    equationsSetSetupInfo%fieldUserNumber=derivedField%userNumber
    equationsSetSetupInfo%field=>derivedField
    !Finish equations set specific setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish the equations set derived creation
    equationsSet%derived%derivedFinished=.TRUE.

    EXITS("EquationsSet_DerivedCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of derived variables field for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedCreateStart
  SUBROUTINE EquationsSet_DerivedCreateStart(equationsSet,derivedFieldUserNumber,derivedField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of a derived field on
    INTEGER(INTG), INTENT(IN) :: derivedFieldUserNumber !<The user specified derived field number
    TYPE(FieldType), POINTER :: derivedField !<If associated on entry, a pointer to the user created derived field which has the same user number as the specified derived field user number. If not associated on entry, on exit, a pointer to the created derived field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: derivedDecomposition,geometricDecomposition
    TYPE(EquationsSetSetupType) :: equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,derivedFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_DerivedCreateStart",err,error,*998)

    CALL EquationsSet_AssertDerivedNotCreated(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"derived",derivedFieldUserNumber,derivedField,err,error,*998)

    CALL EquationsSet_DerivedInitialise(equationsSet,err,error,*999)
    equationsSet%derived%derivedFieldAutoCreated=(.NOT.ASSOCIATED(derivedField))
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_DERIVED_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=derivedFieldUserNumber
    equationsSetSetupInfo%field=>derivedField
    !Start the equations set specfic solution setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(.NOT.equationsSet%derived%derivedFieldAutoCreated) THEN
      equationsSet%derived%derivedField=>derivedField
    ELSE
      equationsSet%derived%derivedField=>derivedField
    ENDIF

    EXITS("EquationsSet_DerivedCreateStart")
    RETURN
999 CALL EquationsSet_DerivedFinalise(equationsSet%derived,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_DerivedCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the derived variables for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedDestroy
  SUBROUTINE EquationsSet_DerivedDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The pointer to the equations set to destroy the derived fields for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_DerivedDestroy",err,error,*999)

    CALL EquationsSet_AssertDerivedIsCreated(equationsSet,err,error,*999)
    
    CALL EquationsSet_DerivedFinalise(equationsSet%derived,err,error,*999)
 
    EXITS("EquationsSet_DerivedDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedDestroy

  !
  !================================================================================================================================
  !

  !>Finalises the derived variables for an equation set and deallocates all memory.
  SUBROUTINE EquationsSet_DerivedFinalise(equationsSetDerived,err,error,*)

    !Argument variables
    TYPE(EquationsSetDerivedType), POINTER :: equationsSetDerived !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_DerivedFinalise",err,error,*999)

    IF(ASSOCIATED(equationsSetDerived)) THEN
      IF(ALLOCATED(equationsSetDerived%variableTypes)) DEALLOCATE(equationsSetDerived%variableTypes)
      DEALLOCATE(equationsSetDerived)
    ENDIF

    EXITS("EquationsSet_DerivedFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the derived variables for a equations set.
  SUBROUTINE EquationsSet_DerivedInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the derived field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError

    ENTERS("EquationsSet_DerivedInitialise",err,error,*999)

    CALL EquationsSet_AssertDerivedNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%derived,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set derived information.",err,error,*998)
    ALLOCATE(equationsSet%derived%variableTypes(EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES),STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set derived variable types.",err,error,*999)
    equationsSet%derived%variableTypes=0
    equationsSet%derived%numberOfVariables=0
    equationsSet%derived%equationsSet=>equationsSet
    equationsSet%derived%derivedFinished=.FALSE.
    equationsSet%derived%derivedFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%derived%derivedField)

    EXITS("EquationsSet_DerivedInitialise")
    RETURN
999 CALL EquationsSet_DerivedFinalise(equationsSet%derived,dummyErr,dummyError,*999)
998 ERRORSEXITS("EquationsSet_DerivedInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the equation set field variables for an equation set and deallocates all memory.
  SUBROUTINE EquationsSet_EquationsFieldFinalise(equationsField,err,error,*)

    !Argument variables
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField !<The pointer to the equations set equations field
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_EquationsFieldFinalise",err,error,*999)

    IF(ASSOCIATED(equationsField)) THEN
      NULLIFY(equationsField%equationsSet)
      equationsField%equationsSetFieldFinished=.FALSE.
      equationsField%equationsSetFieldAutoCreated=.FALSE.
      NULLIFY(equationsField%equationsSetField)
      DEALLOCATE(equationsField)
    ENDIF
    
    EXITS("EquationsSet_EquationsFieldFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsFieldFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsFieldFinalise
  
  !
  !================================================================================================================================
  !
  !>Initialises the equations set field for a equations set.
  SUBROUTINE EquationsSet_EquationsFieldInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the equations field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_EquationsFieldInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    CALL EquationsSet_AssertEquationsFieldNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%equationsField,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set equations field information.",err,error,*999)
    equationsSet%equationsField%equationsSet=>equationsSet
    equationsSet%equationsField%equationsSetFieldFinished=.FALSE.
    equationsSet%equationsField%equationsSetFieldAutoCreated=.TRUE.
    NULLIFY(equationsSet%equationsField%equationsSetField)
        
    EXITS("EquationsSet_EquationsFieldInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsFieldInitialise

  !
  !================================================================================================================================
  !



  !>Sets up the specifices for an equation set.
  SUBROUTINE EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the setup on
    TYPE(EquationsSetSetupType), INTENT(INOUT) ::equationsSetSetupInfo !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_Setup",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
       & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL Bioelectric_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL Fitting_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsSet_Setup")
    RETURN
999 ERRORSEXITS("EquationsSet_Setup",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_Setup

  !
  !================================================================================================================================
  !

 !>Finish the creation of equations for the equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_EquationsCreateFinish
  SUBROUTINE EquationsSet_EquationsCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to finish the creation of the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    
    ENTERS("EquationsSet_EquationsCreateFinish",err,error,*999)

    CALL EquationsSet_AssertEquationsNotFinished(equationsSet,err,error,*999)

    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    !Finish the equations specific solution setup.
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
       
    EXITS("EquationsSet_EquationsCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set. \see OpenCMISS::Iron::cmfe_EquationsSet_EquationsCreateStart
  !>Default values set for the EQUATIONS's attributes are:
  !>- OUTPUT_TYPE: 0 (EQUATIONS_SET_NO_OUTPUT)
  !>- SPARSITY_TYPE: 1 (EQUATIONS_SET_SPARSE_MATRICES)
  !>- NONLINEAR_JACOBIAN_TYPE: 0
  !>- INTERPOLATION: null
  !>- LINEAR_DATA: null 
  !>- NONLINEAR_DATA: null
  !>- TIME_DATA: null
  !>- vectorMapping:  
  !>- vectorMatrices:  
  SUBROUTINE EquationsSet_EquationsCreateStart(equationsSet,equations,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to create equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo

    ENTERS("EquationsSet_EquationsCreateStart",err,error,*999)

    CALL EquationsSet_AssertEquationsNotCreated(equationsSet,err,error,*999)
    IF(ASSOCIATED(equations)) CALL FlagError("Equations is already associated.",err,error,*999)
    
    !Initialise the setup    
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    !Start the equations set specific solution setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Return the pointer
    equations=>equationsSet%equations
       
    EXITS("EquationsSet_EquationsCreateStart")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the equations for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_EquationsDestroy
  SUBROUTINE EquationsSet_EquationsDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_EquationsDestroy",err,error,*999)

    CALL EquationsSet_AssertEquationsIsCreated(equationsSet,err,error,*999)
    
    CALL Equations_Finalise(equationsSet%equations,err,error,*999)
        
    EXITS("EquationsSet_EquationsDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsDestroy

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear equations set.
  SUBROUTINE EquationsSet_JacobianEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_JacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)

    IF(equationsSet%outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set Jacobian evaluate: ",equationsSet%label,err,error,*999)
    ENDIF
    
    SELECT CASE(equations%linearity)
    CASE(EQUATIONS_LINEAR)
      SELECT CASE(equations%timeDependence)
      CASE(EQUATIONS_STATIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleStaticLinearFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_QUASISTATIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleQuasistaticLinearFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "// &
            & TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_AssembleDynamicLinearFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_NONLINEAR)
      SELECT CASE(equations%timeDependence)
      CASE(EQUATIONS_STATIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_JacobianEvaluateStaticFEM(equationsSet,err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL EquationsSet_JacobianEvaluateStaticNodal(equationsSet,err,error,*999)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_QUASISTATIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_JacobianEvaluateStaticFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        ! sebk 15/09/09
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_JacobianEvaluateDynamicFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_TIME_STEPPING)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The equations set time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_NONLINEAR_BCS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsSet_JacobianEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_JacobianEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_JacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an static equations set using the finite element method
  SUBROUTINE EquationsSet_JacobianEvaluateStaticFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalStart,internalFinish,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
  
    ENTERS("EquationsSet_JacobianEvaluateStaticFEM",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(dependentDomain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)

    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx                  
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_JacobianEvaluateStaticFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_JacobianEvaluateStaticFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_JacobianEvaluateStaticFEM

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an dynamic equations set using the finite element method
  SUBROUTINE EquationsSet_JacobianEvaluateDynamicFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalFinish,internalStart,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
  
    ENTERS("EquationsSet_JacobianEvaluateDynamicFEM",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(dependentDomain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
 
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
      
    EXITS("EquationsSet_JacobianEvaluateDynamicFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_JacobianEvaluateDynamicFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_JacobianEvaluateDynamicFEM

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an equations set.
  SUBROUTINE EquationsSet_ResidualEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: linearity,numberOfResiduals,numberOfResidualVariables,outputType,residualIdx,residualVariableIdx, &
      & timeDependence
    TYPE(DistributedVectorType), POINTER :: residualDistributedVector,residualVariableDistributedVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldParameterSetType), POINTER :: residualParameterSet
    TYPE(FieldVariableType), POINTER :: residualVariable
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_ResidualEvaluate",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
        
    IF(outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set residual evaluate: ",equationsSet%label,err,error,*999)
    ENDIF

    CALL Equations_LinearityTypeGet(equations,linearity,err,error,*999)
    CALL Equations_TimeDependenceTypeGet(equations,timeDependence,err,error,*999)
    SELECT CASE(linearity)
    CASE(EQUATIONS_LINEAR)
      CALL FlagError("Can not evaluate a residual for linear equations.",err,error,*999)
    CASE(EQUATIONS_NONLINEAR)
      SELECT CASE(timeDependence)
      CASE(EQUATIONS_STATIC,EQUATIONS_QUASISTATIC)!Quasistatic handled like static
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_ResidualEvaluateStaticFEM(equationsSet,err,error,*999)
        CASE(EQUATIONS_SET_NODAL_SOLUTION_METHOD)
          CALL EquationsSet_ResidualEvaluateStaticNodal(equationsSet,err,error,*999)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(equationsSet%solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          CALL EquationsSet_ResidualEvaluateDynamicFEM(equationsSet,err,error,*999)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%solutionMethod,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set time dependence type of "//TRIM(NumberToVString(timeDependence,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_NONLINEAR_BCS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The equations linearity of "//TRIM(NumberToVString(linearity,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    !Update the residual parameter set if it exists
!!TODO: This is wrong. There should be no assumption that the residual vector (based on lhsVariable) should have any connection
!!      to a residual variable (based on the residualVariable).
    CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)
    DO residualIdx=1,numberOfResiduals
      NULLIFY(residualMapping)
      CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,residualIdx,residualMapping,err,error,*999)
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
      NULLIFY(residualDistributedVector)
      CALL EquationsMatricesResidual_DistributedVectorGet(residualVector,EQUATIONS_MATRICES_CURRENT_VECTOR, &
        & residualDistributedVector,err,error,*999)
      CALL EquationsMappingResidual_NumberOfResidualVariablesGet(residualMapping,numberOfResidualVariables,err,error,*999)
      DO residualVariableIdx=1,numberOfResidualVariables
        NULLIFY(residualVariable)
        CALL EquationsMappingResidual_VariableGet(residualMapping,residualVariableIdx,residualVariable,err,error,*999)
        NULLIFY(residualParameterSet)
        CALL FieldVariable_ParameterSetExists(residualVariable,FIELD_RESIDUAL_SET_TYPE,residualParameterSet,err,error,*999)
        IF(ASSOCIATED(residualParameterSet)) THEN
          NULLIFY(residualVariableDistributedVector)
          CALL FieldParameterSet_ParametersGet(residualParameterSet,residualVariableDistributedVector,err,error,*999)
          !Residual parameter set exists. Copy the residual vector to the residuals parameter set.
          CALL DistributedVector_Copy(residualDistributedVector,residualVariableDistributedVector,1.0_DP,err,error,*999)
        ENDIF
      ENDDO !residualVariableIdx
    ENDDO !residualIdx
       
    EXITS("EquationsSet_ResidualEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_ResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_ResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an dynamic equations set using the finite element method
  SUBROUTINE EquationsSet_ResidualEvaluateDynamicFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalFinish,internalStart,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
 
    ENTERS("EquationsSet_ResidualEvaluateDynamicFEM",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(dependentDomain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)
 
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx                  
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
      
    EXITS("EquationsSet_ResidualEvaluateDynamicFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_ResidualEvaluateDynamicFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_ResidualEvaluateDynamicFEM

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an static equations set using the finite element method
  SUBROUTINE EquationsSet_ResidualEvaluateStaticFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,elementIdx,element,ghostFinish,internalFinish,internalStart,numberOfTimes,outputType
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: dependentDecomposition
    TYPE(DomainType), POINTER :: dependentDomain
    TYPE(DomainMappingType), POINTER :: elementsMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
 
    ENTERS("EquationsSet_ResidualEvaluateStaticFEM",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(dependentDecomposition)
    CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
    NULLIFY(dependentDomain)
    CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(dependentDomain,domainMappings,err,error,*999)
    NULLIFY(elementsMapping)
    CALL DomainMappings_ElementsMappingGet(domainMappings,elementsMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(elementsMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(elementsMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(elementsMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(elementsMapping,ghostFinish,err,error,*999)

    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      elementUserElapsed=0.0_SP
      elementSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal elements
    DO elementIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(elementsMapping,elementIdx,element,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_ResidualEvaluateStaticFEM")
    RETURN
999 ERRORSEXITS("EquationsSet_ResidualEvaluateStaticFEM",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_ResidualEvaluateStaticFEM

  !
  !================================================================================================================================
  !

  !>Finalises the equations set setup and deallocates all memory
  SUBROUTINE EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*)

    !Argument variables
    TYPE(EquationsSetSetupType), INTENT(OUT) :: equationsSetSetupInfo !<The equations set setup to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SetupFinalise",err,error,*999)

    equationsSetSetupInfo%setupType=0
    equationsSetSetupInfo%actionType=0
    equationsSetSetupInfo%fieldUserNumber=0
    NULLIFY(equationsSetSetupInfo%field)
    equationsSetSetupInfo%analyticFunctionType=0
    
    EXITS("EquationsSet_SetupFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_SetupFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SetupFinalise
  
  !
  !================================================================================================================================
  !

  !>Initialise the equations set setup.
  SUBROUTINE EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*)

    !Argument variables
    TYPE(EquationsSetSetupType), INTENT(OUT) :: equationsSetSetupInfo !<The equations set setup to be initialised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SetupInitialise",err,error,*999)

    equationsSetSetupInfo%setupType=0
    equationsSetSetupInfo%actionType=0
    equationsSetSetupInfo%fieldUserNumber=0
    NULLIFY(equationsSetSetupInfo%field)
    equationsSetSetupInfo%analyticFunctionType=0
    
    EXITS("EquationsSet_SetupInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_SetupInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SetupInitialise
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for an equations set.
  SUBROUTINE EquationsSet_OutputTypeSet(equationsSet,outputType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see EquationsSetRoutines_OutputTypes,EquationsSetRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_OutputTypeSet",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)

    SELECT CASE(outputType)
    CASE(EQUATIONS_SET_NO_OUTPUT)
      equationsSet%outputType=EQUATIONS_SET_NO_OUTPUT
    CASE(EQUATIONS_SET_PROGRESS_OUTPUT)
      equationsSet%outputType=EQUATIONS_SET_PROGRESS_OUTPUT
    CASE DEFAULT
      localError="The specified output type of "//TRIM(NumberToVString(outputType,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("EquationsSet_OutputTypeSet")
    RETURN
999 ERRORSEXITS("EquationsSet_OutputTypeSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_OutputTypeSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SolutionMethodSet
  SUBROUTINE EquationsSet_SolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The equations set solution method to set \see EquationsSetRoutines_SolutionMethods,EquationsSetRoutines
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_SolutionMethodSet",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)
    
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))      
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      IF(SIZE(equationsSet%specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation set.", &
          & err,error,*999)
      END IF
      IF(equationsSet%specification(2) == EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
        CALL Monodomain_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
      ELSE
        CALL Bioelectric_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*999)
      END IF
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_EquationsSetSolnMethodSet(equationsSet,solutionMethod,err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("EquationsSet_SolutionMethodSet")
    RETURN
999 ERRORSEXITS("EquationsSet_SolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SolutionMethodSet
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of a source for an equation set. \see OpenCMISS::Iron::cmfe_EquationsSet_SourceCreateFinish
  SUBROUTINE EquationsSet_SourceCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of a souce for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: sourceField

    ENTERS("EquationsSet_SourceCreateFinish",err,error,*999)

    CALL EquationsSet_AssertSourceNotFinished(equationsSet,err,error,*999)
    
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_SOURCE_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_FINISH_ACTION
    NULLIFY(sourceField)
    CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
    equationsSetSetupInfo%fieldUserNumber=sourceField%userNumber
    equationsSetSetupInfo%field=>sourceField
    !Finish the equation set specific source setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Finish the source creation
    equationsSet%source%sourceFinished=.TRUE.
       
    EXITS("EquationsSet_SourceCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_SourceCreateFinish",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of a source for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SourceCreateStart
  SUBROUTINE EquationsSet_SourceCreateStart(equationsSet,sourceFieldUserNumber,sourceField,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to start the creation of a source for
    INTEGER(INTG), INTENT(IN) :: sourceFieldUserNumber !<The user specified source field number
    TYPE(FieldType), POINTER :: sourceField !<If associated on entry, a pointer to the user created source field which has the same user number as the specified source field user number. If not associated on entry, on exit, a pointer to the created source field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: geometricDecomposition,sourceDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,sourceFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_SourceCreateStart",err,error,*998)

    CALL EquationsSet_AssertSourceNotCreated(equationsSet,err,error,*998)
    CALL EquationsSet_FieldRegionSetupCheck(equationsSet,"source",sourceFieldUserNumber,sourceField,err,error,*998)
    
    !Initialise the equations set source
    CALL EquationsSet_SourceInitialise(equationsSet,err,error,*999)
    equationsSet%source%sourceFieldAutoCreated=(.NOT.ASSOCIATED(sourceField))
    !Initialise the setup
    CALL EquationsSet_SetupInitialise(equationsSetSetupInfo,err,error,*999)
    equationsSetSetupInfo%setupType=EQUATIONS_SET_SETUP_SOURCE_TYPE
    equationsSetSetupInfo%actionType=EQUATIONS_SET_SETUP_START_ACTION
    equationsSetSetupInfo%fieldUserNumber=sourceFieldUserNumber
    equationsSetSetupInfo%field=>sourceField
    !Start the equation set specific source setup
    CALL EquationsSet_Setup(equationsSet,equationsSetSetupInfo,err,error,*999)
    !Finalise the setup
    CALL EquationsSet_SetupFinalise(equationsSetSetupInfo,err,error,*999)
    !Set pointers
    IF(equationsSet%source%sourceFieldAutoCreated) THEN            
      sourceField=>equationsSet%source%sourceField
    ELSE
      equationsSet%source%sourceField=>sourceField
    ENDIF
       
    EXITS("EquationsSet_SourceCreateStart")
    RETURN
999 CALL EquationsSet_SourceFinalise(equationsSet%source,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_SourceCreateStart",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the source for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SourceDestroy
  SUBROUTINE EquationsSet_SourceDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to destroy the source for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SourceDestroy",err,error,*999)

    CALL EquationsSet_AssertSourceIsCreated(equationsSet,err,error,*999)
    
    CALL EquationsSet_SourceFinalise(equationsSet%source,err,error,*999)
       
    EXITS("EquationsSet_SourceDestroy")
    RETURN
999 ERRORSEXITS("EquationsSet_SourceDestroy",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceDestroy

  !
  !================================================================================================================================
  !

  !>Finalise the source for a equations set and deallocate all memory.
  SUBROUTINE EquationsSet_SourceFinalise(equationsSetSource,err,error,*)

    !Argument variables
    TYPE(EquationsSetSourceType), POINTER :: equationsSetSource !<A pointer to the equations set source to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SourceFinalise",err,error,*999)

    IF(ASSOCIATED(equationsSetSource)) DEALLOCATE(equationsSetSource)
       
    EXITS("EquationsSet_SourceFinalise")
    RETURN
999 ERRORSEXITS("EquationsSet_SourceFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceFinalise

  !
  !================================================================================================================================
  !

  !>Initialises the source for an equations set.
  SUBROUTINE EquationsSet_SourceInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the source field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EquationsSet_SourceInitialise",err,error,*998)

    CALL EquationsSet_AssertSourceNotCreated(equationsSet,err,error,*999)
    
    ALLOCATE(equationsSet%source,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set source.",err,error,*999)
    equationsSet%source%equationsSet=>equationsSet
    equationsSet%source%sourceFinished=.FALSE.
    equationsSet%source%sourceFieldAutoCreated=.FALSE.
    NULLIFY(equationsSet%source%sourceField)
        
    EXITS("EquationsSet_SourceInitialise")
    RETURN
999 CALL EquationsSet_SourceFinalise(equationsSet%source,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_SourceInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SourceInitialise

  !
  !================================================================================================================================
  !

  !>Sets/changes the current times for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_TimesSet
  SUBROUTINE EquationsSet_TimesSet(equationsSet,currentTime,deltaTime,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the times for
    REAL(DP), INTENT(IN) :: currentTime !<The current time for the equations set to set.
    REAL(DP), INTENT(IN) :: deltaTime !<The current time incremenet for the equations set to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_TimesSet",err,error,*999) 

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(currentTime<0.0_DP) CALL FlagError("Invalid current time. The time must be >= zero.",err,error,*999)

    equationsSet%currentTime=currentTime
    equationsSet%deltaTime=deltaTime
      
    EXITS("EquationsSet_TimesSet")
    RETURN
999 ERRORSEXITS("EquationsSet_TimesSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_TimesSet
  
  !
  !================================================================================================================================
  !

  !>Calculates a derived variable value for the equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedVariableCalculate
  SUBROUTINE EquationsSet_DerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to calculate. \see EquationsSetRoutines_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_DerivedVariableCalculate",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))  
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_EquationsSetDerivedVariableCalculate(equationsSet,derivedType,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "//TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("EquationsSet_DerivedVariableCalculate")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedVariableCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedVariableCalculate

  !
  !================================================================================================================================
  !

  !>Sets the field variable type of the derived field to be used to store a derived variable. \see OpenCMISS::Iron::cmfe_EquationsSet_DerivedVariableSet
  SUBROUTINE EquationsSet_DerivedVariableSet(equationsSet,derivedType,variableType,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate a derived field for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to calculate. \see EquationsSetRoutines_DerivedTypes.
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type used to store the calculated derived value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FieldType), POINTER :: derivedField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_DerivedVariableSet",err,error,*999)

    !Check pointers and finished state
    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    NULLIFY(derivedField)
    CALL EquationsSet_DerivedFieldGet(equationsSet,derivedField,err,error,*999)
    IF(equationsSet%derived%derivedFinished) CALL FlagError("Equations set derived information is already finished.",err,error,*999)
    IF(derivedType<1.OR.derivedType>EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES) THEN
      localError="The specified derived variable type of "//TRIM(NumberToVString(derivedType,"*",err,error))// &
        & " is invalid. It should be between >= 1 and <= "//TRIM(NumberToVString(EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES,"*", &
        & err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    NULLIFY(fieldVariable)
    CALL Field_VariableGet(derivedField,variableType,fieldVariable,err,error,*999)
    
    IF(equationsSet%derived%variableTypes(derivedType)==0) &
      & equationsSet%derived%numberOfVariables=equationsSet%derived%numberOfVariables+1
    equationsSet%derived%variableTypes(derivedType)=variableType

    EXITS("EquationsSet_DerivedVariableSet")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedVariableSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedVariableSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set. \see OpenCMISS::Iron::cmfe_EquationsSet_SpecificationSet
  SUBROUTINE EquationsSet_SpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification array to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: specificationIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_SpecificationSet",err,error,*999)

    CALL EquationsSet_AssertNotFinished(equationsSet,err,error,*999)
    IF(SIZE(specification,1)<1) CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL ClassicalField_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      IF(SIZE(specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation class.", &
          & err,error,*999)
      END IF
      IF(specification(2)==EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
        CALL Monodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      ELSE
        CALL Bioelectric_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
      END IF
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MultiPhysics_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL Fitting_EquationsSetSpecificationSet(equationsSet,specification,err,error,*999)
    CASE(EQUATIONS_SET_OPTIMISATION_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVstring(specification(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set the specification length
    equationsSet%specificationLength=0
    DO specificationIdx=1,SIZE(equationsSet%specification,1)
      IF(equationsSet%specification(specificationIdx)>0) THEN
        equationsSet%specificationLength=specificationIdx
      ENDIF
    ENDDO !specificationIdx
    
    EXITS("EquationsSet_SpecificationSet")
    RETURN
999 ERRORS("EquationsSet_SpecificationSet",err,error)
    EXITS("EquationsSet_SpecificationSet")
    RETURN 1
    
  END SUBROUTINE EquationsSet_SpecificationSet
  
  !
  !================================================================================================================================
  !

  !>Evaluate a tensor at a given element Gauss point.
  SUBROUTINE EquationsSet_TensorInterpolateGaussPoint(equationsSet,tensorEvaluateType,gaussPointNumber,userElementNumber,values, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number of the field to interpolate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_TensorInterpolateGaussPoint",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) CALL FlagError("Equations set specification must have at least one entry.", &
      & err,error,*999)

    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_TensorInterpolateGaussPoint(equationsSet,tensorEvaluateType,gaussPointNumber,userElementNumber,values, &
        & err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_OPTIMISATION_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVstring(equationsSet%specification(1),"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)      
    END SELECT

    EXITS("EquationsSet_TensorInterpolateGaussPoint")
    RETURN
999 ERRORSEXITS("EquationsSet_TensorInterpolateGaussPoint",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_TensorInterpolateGaussPoint

  !
  !================================================================================================================================
  !

  !>Evaluate a tensor at a given element xi location.
  SUBROUTINE EquationsSet_TensorInterpolateXi(equationsSet,tensorEvaluateType,userElementNumber,xi,values,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(IN) :: xi(:) !<The element xi to interpolate the field at.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_TensorInterpolateXi",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) CALL FlagError("Equations set specification must have at least one entry.", &
      & err,error,*999)

    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_TensorInterpolateXi(equationsSet,tensorEvaluateType,userElementNumber,xi,values,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_OPTIMISATION_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      CALL FlagError("The first equations set specification of "// &
        & TRIM(NumberToVstring(equationsSet%specification(1),"*",err,error))// &
        & " is not valid.",err,error,*999)
    END SELECT

    EXITS("EquationsSet_TensorInterpolateXi")
    RETURN
999 ERRORSEXITS("EquationsSet_TensorInterpolateXi",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_TensorInterpolateXi

  !
  !================================================================================================================================
  !

  !>Finalises all equations sets on a region and deallocates all memory.
  SUBROUTINE EquationsSets_Finalise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to finalise the problems for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSets_Finalise",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*999)
    
    IF(ASSOCIATED(region%equationsSets)) THEN
      DO WHILE(region%equationsSets%numberOfEquationsSets>0)
        CALL EquationsSet_Destroy(region%equationsSets%equationsSets(1)%ptr,err,error,*999)
      ENDDO 
      DEALLOCATE(region%equationsSets)
    ENDIF

    EXITS("EquationsSets_Finalise")
    RETURN
999 ERRORSEXITS("EquationsSets_Finalise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSets_Finalise

  !
  !================================================================================================================================
  !

  !>Intialises all equations sets on a region.
  SUBROUTINE EquationsSets_Initialise(region,err,error,*)

    !Argument variables
    TYPE(RegionType), POINTER :: region !<A pointer to the region to initialise the equations sets for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSets_Initialise",err,error,*999)

    IF(.NOT.ASSOCIATED(region)) CALL FlagError("Region is not associated.",err,error,*998)
    IF(ASSOCIATED(region%equationsSets)) CALL FlagError("Region already has associated equations sets",err,error,*998)
    
!!TODO: Inherit any equations sets from the parent region???
    ALLOCATE(region%equationsSets,STAT=err)
    IF(err/=0) CALL FlagError("Could not allocate region equations sets",err,error,*999)
    region%equationsSets%region=>region
    region%equationsSets%numberOfEquationsSets=0
    NULLIFY(region%equationsSets%equationsSets)
 
    EXITS("EquationsSets_Initialise")
    RETURN
999 IF(ASSOCIATED(region%equationsSets)) DEALLOCATE(region%equationsSets)
998 ERRORSEXITS("EquationsSets_Initialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSets_Initialise
  
  !
  !================================================================================================================================
  !

  !> Apply the boundary condition load increment to dependent field
  SUBROUTINE EquationsSet_BoundaryConditionsIncrement(equationsSet,boundaryConditions,iterationNumber, &
    & maximumNumberOfIterations,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<The equations set to increment the boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<The boundary conditions to apply the increment to
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: conditionGlobalDOF,conditionIdx,conditionLocalDOF,conditionType,dirichletIdx,dirichletDOFIdx, &
      & dirichletDomainNumber,fixedIncrementedCount,ghostStart,globalDirichletDOFIdx,globalNeumannDOFIdx,globalPressureIncDOFIdx, &
      & localDirichletDOFIdx,localNeumannDOFIdx,localNeumannPointDOFIdx,localPressureIncDOFIdx,movedWallIncrementedCount, &
      & myGroupComputationNodeNumber,neumannDomainNumber,neumannIdx,neumannPointCount,neumannPointIncrementedCount, &
      & numberOfDirichletConditions,numberOfVariables,pressureIncDomainNumber,pressureIncIdx,pressureIncrementedCount, &
      & variableIdx,variableType
    REAL(DP) :: fullLoad, currentLoad, newLoad, prevLoad
    REAL(DP), POINTER :: fullLoads(:),currentLoads(:), prevLoads(:)
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions
    TYPE(BoundaryConditionsVariableType),   POINTER :: boundaryConditionsVariable
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DistributedVectorType), POINTER :: neumannPointDOFValues
    TYPE(DomainMappingType), POINTER :: domainMapping,neumannPointDOFMapping
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("EquationsSet_BoundaryConditionsIncrement",err,error,*999)
   
    !Take the stored load, scale it down appropriately then apply to the unknown variables    
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    
    IF(diagnostics1) THEN
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set BC increment: ",equationsSet%label,err,error,*999)
    ENDIF

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    !Loop over the variables associated with this equations set
    !\todo: Looping over all field variables is not safe when volume-coupled problem is solved. Look at matrix and rhs mapping instead?
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(dependentVariable,domainMapping,err,error,*999)
        CALL DomainMapping_GhostStartGet(domainMapping,ghostStart,err,error,*999)
        ! Check if there are any incremented conditions applied for this boundary conditions variable
        CALL BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,BOUNDARY_CONDITION_FIXED_INCREMENTED, &
          & fixedIncrementedCount,err,error,*999)
        CALL BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED, &
          & movedWallIncrementedCount,err,error,*999)
        IF(fixedIncrementedCount>0.OR.movedWallIncrementedCount>0) THEN
          NULLIFY(dirichletBoundaryConditions)
          CALL BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,dirichletBoundaryConditions, &
            & err,error,*999)
          !Get the pointer to vector holding the full and current loads
          !   full load: FIELD_BOUNDARY_CONDITIONS_SET_TYPE - holds the target load values
          !   current load: FIELD_VALUES_SET_TYPE - holds the current increment values
          NULLIFY(fullLoads)
          CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,fullLoads,err,error,*999)
          !chrm 22/06/2010: 'FIELD_BOUNDARY_CONDITIONS_SET_TYPE' does not get updated with time (update_BCs)
          !\ToDo: How can this be achieved ???
          NULLIFY(currentLoads)
          CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_VALUES_SET_TYPE,currentLoads,err,error,*999)
          !Get full increment, calculate new load, then apply to dependent field
          CALL BoundaryConditionsVariable_NumberOfDirichletConditionsGet(boundaryConditionsVariable,numberOfDirichletConditions, &
            & err,error,*999)
          NULLIFY(dirichletBoundaryConditions)
          IF(numberOfDirichletConditions>0) CALL BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable, &
            & dirichletBoundaryConditions,err,error,*999)
          DO dirichletIdx=1,numberOfDirichletConditions
            CALL BoundaryConditionsDirichlet_DirichletDOFIndexGet(dirichletBoundaryConditions,dirichletIdx,globalDirichletDOFIdx, &
              & err,error,*999)
            !Check whether we have an incremented boundary condition type
            CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalDirichletDOFIdx,conditionType, &
              & err,error,*999)
            SELECT CASE(conditionType)
            CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
              !Convert dof index to local index
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,globalDirichletDOFIdx,1,dirichletDomainNumber, &
                & err,error,*999)
              IF(dirichletDomainNumber==myGroupComputationNodeNumber) THEN
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,globalDirichletDOFIdx,1,localDirichletDOFIdx, &
                  & err,error,*999)
                IF(localDirichletDOFIdx>=1.AND.localDirichletDOFIdx<ghostStart) THEN
                  fullLoad=fullLoads(localDirichletDOFIdx)
                  ! Apply full load if last step, or fixed BC
                  IF(iterationNumber==maximumNumberOfIterations) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDirichletDOFIdx, &
                      & fullLoad,err,error,*999)
                  ELSE
                    !Calculate new load and apply to dependent field
                    currentLoad=currentLoads(localDirichletDOFIdx)
                    newLoad=currentLoad+(fullLoad-currentLoad)/(maximumNumberOfIterations-iterationNumber+1)
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDirichletDOFIdx, &
                      & newLoad,err,error,*999)
                    IF(diagnostics1) THEN
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF idx : ",localDirichletDOFIdx,err,error,*999)
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Current load = ",currentLoad,err,error,*999)
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    New load = ",newLoad,err,error,*999)
                    ENDIF
                  ENDIF !Full or intermediate load
                ENDIF !non-ghost dof
              ENDIF !current domain
            CASE DEFAULT
              !Do nothing for non-incremented boundary conditions
            END SELECT
          ENDDO !dirichletIdx
          !\ToDo: What happens if the call below is issued
          !without actually that the dependent field has been modified in above conditional ?
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Restore the vector handles
          CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE,fullLoads,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_VALUES_SET_TYPE,currentLoads,err,error,*999)
        ENDIF
        ! Also increment any incremented Neumann point conditions
        CALL BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED, &
          & neumannPointIncrementedCount,err,error,*999)
        CALL BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,BOUNDARY_CONDITION_NEUMANN_POINT, &
          & neumannPointCount,err,error,*999)
        IF(neumannPointIncrementedCount>0) THEN
          NULLIFY(neumannBoundaryConditions)
          CALL BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,neumannBoundaryConditions,err,error,*999)
          NULLIFY(neumannPointDOFMapping)
          CALL BoundaryConditionsNeumann_PointDOFMappingGet(neumannBoundaryConditions,neumannPointDOFMapping,err,error,*999)
          NULLIFY(neumannPointDOFValues)
          CALL BoundaryConditionsNeumann_PointDOFValuesGet(neumannBoundaryConditions,neumannPointDOFValues,err,error,*999)
          ! The boundary conditions parameter set contains the full values and the
          ! current incremented values are transferred to the point values vector
          DO neumannIdx=1,neumannPointIncrementedCount+neumannPointCount
            CALL BoundaryConditionsNeumann_NeumannDOFIndexGet(neumannBoundaryConditions,neumannIdx,globalNeumannDOFIdx, &
              & err,error,*999)
            CALL BoundaryConditionsVariable_ConditionTypeGet(boundaryConditionsVariable,globalNeumannDOFIdx,conditionType, &
             & err,error,*999)
            ! conditionGlobalDOF could be for non-incremented point Neumann condition
            IF(conditionType/=BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) CYCLE
            CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,globalDirichletDOFIdx,1,dirichletDomainNumber, &
              & err,error,*999)
            IF(neumannDomainNumber==myGroupComputationNodeNumber) THEN
              CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,globalNeumannDOFIdx,1,localNeumannDOFIdx, &
                & err,error,*999)
              CALL DomainMapping_LocalNumberFromGlobalGet(neumannPointDOFMapping,neumannIdx,1,localNeumannPointDOFIdx, &
                & err,error,*999)
              CALL FieldVariable_ParameterSetGetLocalDOF(dependentVariable,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                & localNeumannDOFIdx,fullLoad,err,error,*999)
              CALL DistributedVector_ValuesSet(neumannPointDOFValues,localNeumannPointDOFIdx,fullLoad*(REAL(iterationNumber)/ &
                & REAL(maximumNumberOfIterations)),err,error,*999)
            ENDIF
          ENDDO !conditionIdx
        ENDIF

        !There might also be pressure incremented conditions
        CALL BoundaryConditionsVariable_DOFCountGet(boundaryConditionsVariable,BOUNDARY_CONDITION_PRESSURE_INCREMENTED, &
          & pressureIncrementedCount,err,error,*999)
        IF(pressureIncrementedCount>0) THEN
          ! handle pressure incremented boundary conditions
          NULLIFY(pressureIncrementedBoundaryConditions)
          CALL BoundaryConditionsVariable_PressureIncConditionsGet(boundaryConditionsVariable, &
            & pressureIncrementedBoundaryConditions,err,error,*999)
          !Due to a variety of reasons, the pressure incremented type is setup differently to dirichlet conditions.
          !We store two sets of vectors, the current and previous values
          !   current: FIELD_PRESSURE_VALUES_SET_TYPE - always holds the current increment, even if not incremented
          !   previous: FIELD_PREVIOUS_PRESSURE_SET_TYPE - holds the previously applied increment
          !Grab the pointers for both
          NULLIFY(prevLoads)
          CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE,prevLoads,err,error,*999)
          NULLIFY(currentLoads)
          CALL FieldVariable_ParameterSetDataGet(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,currentLoads,err,error,*999)
          !Calculate the new load, update the old load
          IF(iterationNumber==1) THEN
            !On the first iteration, FIELD_PRESSURE_VALUES_SET_TYPE actually contains the full load
            DO pressureIncIdx=1,pressureIncrementedCount
              !Global dof index
              CALL BoundaryConditionsPressureInc_PressureIncDOFIndexGet(pressureIncrementedBoundaryConditions,pressureIncIdx, &
                & globalPressureIncDOFIdx,err,error,*999)
              !Must convert into local dof index
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,globalPressureIncDOFIdx,1,pressureIncDomainNumber, &
                & err,error,*999)
              IF(pressureIncDomainNumber==myGroupComputationNodeNumber) THEN
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,globalPressureIncDOFIdx,1,localPressureIncDOFIdx, &
                  & err,error,*999)
                IF(localPressureIncDOFIdx>=1.AND.localPressureIncDOFIdx<ghostStart) THEN
                  newLoad=currentLoads(localPressureIncDOFIdx)
                  newLoad=newLoad/maximumNumberOfIterations
                  !Update current and previous loads
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE, &
                    & localPressureIncDOFIdx,newLoad,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                    & localPressureIncDOFIdx,0.0_dp,err,error,*999)
                  IF(diagnostics1) THEN
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF Idx : ",localPressureIncDOFIdx,err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Current load = ",currentLoads(localPressureIncDOFIdx), &
                      & err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    New load = ",newLoad,err,error,*999)
                  ENDIF
                ENDIF !Non-ghost dof
              ENDIF !Current domain
            ENDDO !conditionIdx
          ELSE
            !Calculate the new load, keep the current load
            DO conditionIdx=1,pressureIncrementedCount
              !This is global dof idx
              CALL BoundaryConditionsPressureInc_PressureIncDOFIndexGet(pressureIncrementedBoundaryConditions,pressureIncIdx, &
                & globalPressureIncDOFIdx,err,error,*999)
              !Must convert into local dof index
              CALL DomainMapping_DomainNumberFromGlobalGet(domainMapping,globalPressureIncDOFIdx,1,pressureIncDomainNumber, &
                & err,error,*999)
              IF(pressureIncDomainNumber==myGroupComputationNodeNumber) THEN
                CALL DomainMapping_LocalNumberFromGlobalGet(domainMapping,globalPressureIncDOFIdx,1,localPressureIncDOFIdx, &
                  & err,error,*999)
                IF(localPressureIncDOFIdx>=1.AND.localPressureIncDOFIdx<ghostStart) THEN
                  prevLoad=prevLoads(localPressureIncDOFIdx)
                  currentLoad=currentLoads(localPressureIncDOFIdx)
                  newLoad=currentLoad+(currentLoad-prevLoad)  !This may be subject to numerical errors...
                  !if (conditionIdx==1) write(*,*) "new load=",new_load
                  !Update current and previous loads
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE, &
                    & localPressureIncDOFIdx,newLoad,err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                    & localPressureIncDOFIdx,currentLoad,err,error,*999)
                  IF(diagnostics1) THEN
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  DOF idx : ",localPressureIncDOFIdx,err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    Current load = ",currentLoads(localPressureIncDOFIdx), &
                      & err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    New load = ",newLoad,err,error,*999)
                  ENDIF
                ENDIF !Non-ghost dof
              ENDIF !Current domain
            ENDDO !conditionIdx
          ENDIF
          !Start transfer of dofs to neighbouring domains
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
          !Restore the vector handles
          CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE,prevLoads,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,currentLoads,err,error,*999)
          !Finish transfer of dofs to neighbouring domains
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
        ENDIF !Pressure incremented bc block                
      ELSE
        ! do nothing - no boundary conditions variable type associated?
      ENDIF
    ENDDO !variableIdx

    EXITS("EquationsSet_BoundaryConditionsIncrement")
    RETURN
999 ERRORSEXITS("EquationsSet_BoundaryConditionsIncrement",err,error)
    RETURN 1

  END SUBROUTINE EquationsSet_BoundaryConditionsIncrement

  !
  !================================================================================================================================
  !

  !> Apply load increments for equations sets
  SUBROUTINE EquationsSet_LoadIncrementApply(equationsSet,boundaryConditions,iterationNumber,maximumNumberOfIterations, &
    & err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<The boundary conditions to apply the increment to
    INTEGER(INTG), INTENT(IN) :: iterationNumber !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: maximumNumberOfIterations !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_LoadIncrementApply",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
   
    !Increment boundary conditions
    CALL EquationsSet_BoundaryConditionsIncrement(equationsSet,boundaryConditions,iterationNumber,maximumNumberOfIterations, &
      & err,error,*999)

    !Apply any other equation set specific increments
    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL Elasticity_LoadIncrementApply(equationsSet,iterationNumber,maximumNumberOfIterations,err,error,*999)
    CASE DEFAULT
      !Do nothing
    END SELECT

    EXITS("EquationsSet_LoadIncrementApply")
    RETURN
999 ERRORSEXITS("EquationsSet_LoadIncrementApply",err,error)
    RETURN 1

  END SUBROUTINE EquationsSet_LoadIncrementApply

  !
  !================================================================================================================================
  !

  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear static equations set using a nodal method
  SUBROUTINE EquationsSet_AssembleStaticNonlinearNodal(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,ghostFinish,internalFinish,internalStart,nodeIdx,nodeNumber,numberOfTimes,outputType
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: nodalMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    
    ENTERS("EquationsSet_AssembleStaticNonlinearNodal",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(nodalMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodalMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(nodalMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(nodalMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(nodalMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(nodalMapping,ghostFinish,err,error,*999)

    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      nodeUserElapsed=0.0_SP
      nodeSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal nodes
    DO nodeIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(nodalMapping,nodeIdx,nodeNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      nodeUserElapsed=userElapsed
      nodeSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost nodes
    DO nodeIdx=nodalMapping%boundaryStart,nodalMapping%ghostFinish
      nodeNumber=nodalMapping%domainList(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average nodes equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_AssembleStaticNonlinearNodal")
    RETURN
999 ERRORSEXITS("EquationsSet_AssembleStaticNonlinearNodal",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_AssembleStaticNonlinearNodal

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal Jacobian for the given node number for a nodal equations set.
  SUBROUTINE EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The node number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: jacobianCalculationType,matrixIdx,numberOfJacobians,numberOfResiduals,outputType,residualIdx
    LOGICAL :: updateJacobian
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(NodalMatrixType), POINTER :: nodalMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_NodalJacobianEvaluate",err,error,*999)

   
    NULLIFY(equations)    
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    CALL EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*999)

    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
     DO residualIdx=1,numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)    
      CALL EquationsMatricesResidual_NumberOfJacobiansGet(residualvector,numberOfJacobians,err,error,*999)
      DO matrixIdx=1,numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
        CALL JacobianMatrix_CalculationTypeGet(jacobianMatrix,jacobianCalculationType,err,error,*999)
        SELECT CASE(jacobianCalculationType)
        CASE(EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
          ! None of these routines currently support calculating off diagonal terms for coupled problems,
          ! but when one does we will have to pass through the matrixIdx parameter
          IF(matrixIdx>1) CALL FlagError("Analytic off-diagonal Jacobian calculation not implemented.",err,error,*999)
          SELECT CASE(equationsSet%specification(1))
          CASE(EQUATIONS_SET_ELASTICITY_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
            CALL FluidMechanics_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
          CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MODAL_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The first equations set specification of "// &
              & TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))//" is not valid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The Jacobian calculation type of "//TRIM(NumberToVString(jacobianCalculationType,"*",err,error))// &
            & " is not valid for matrix index number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END DO !matrixIdx
    ENDDO !residualIdx
    IF(outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal Jacobian matrices:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node number = ",nodeNumber,err,error,*999)      
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of residuals = ",numberOfResiduals,err,error,*999)      
      DO residualIdx=1,numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Residual number : ",residualIdx,err,error,*999)
        CALL EquationsMatricesResidual_NumberOfJacobiansGet(residualvector,numberOfJacobians,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of Jacobians = ",numberOfJacobians,err,error,*999)             
        DO matrixIdx=1,numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Jacobian matrix : ",matrixIdx,err,error,*999)
          CALL JacobianMatrix_UpdateMatrixGet(jacobianMatrix,updateJacobian,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",updateJacobian,err,error,*999)
          IF(updateJacobian) CALL JacobianMatrix_NodalMatrixOutput(GENERAL_OUTPUT_TYPE,jacobianMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDDO !residualIdx
    ENDIF

    EXITS("EquationsSet_NodalJacobianEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_NodalJacobianEvaluate",err,error)
    RETURN 1

  END SUBROUTINE EquationsSet_NodalJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the nodal residual and rhs vector for the given node number for a nodal equations set.
  SUBROUTINE EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The nodal number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx,numberOfDynamicMatrices,numberOfLinearMatrices,numberOfResiduals,numberOfSources,outputType, &
      & residualIdx,sourceIdx
    LOGICAL :: updateMatrix,updateVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: dynamicMatrix,linearMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(NodalMatrixType), POINTER :: nodalMatrix
    TYPE(NodalVectorType), POINTER :: nodalVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_NodalResidualEvaluate",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    CALL EquationsMappingNonlinear_NumberOfResidualsGet(nonlinearMapping,numberOfResiduals,err,error,*999)

    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    SELECT CASE(equationsSet%specification(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FluidMechanics_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "//TRIM(NumberToVString(equationsSet%specification(1),"*",err,error))// &
        & " is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    DO residualIdx=1,numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
      residualVector%nodalResidualCalculated=nodeNumber
    ENDDO !residualIdx
    
    IF(outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal matrices and vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node number = ",nodeNumber,err,error,*999)
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",err,error,*999)
        CALL EquationsMatricesDynamic_NumberOfDynamicMatricesGet(dynamicMatrices,numberOfDynamicMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of dynamic matrices = ",numberOfDynamicMatrices,err,error,*999)
        DO matrixIdx=1,numberOfDynamicMatrices
          NULLIFY(dynamicMatrix)
          CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,matrixIdx,dynamicMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Dynamic matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(dynamicMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_NodalMatrixOutput(GENERAL_OUTPUT_TYPE,dynamicMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesExists(vectorMatrices,linearMatrices,err,error,*999)
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL EquationsMatricesLinear_NumberOfLinearMatricesGet(linearMatrices,numberOfLinearMatrices,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of linear matrices = ",numberOfLinearMatrices,err,error,*999)
        DO matrixIdx=1,numberOfLinearMatrices
          NULLIFY(linearMatrix)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,matrixIdx,linearMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Linear matrix : ",matrixIdx,err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(linearMatrix,updateMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",updateMatrix,err,error,*999)
          IF(updateMatrix) CALL EquationsMatrix_ElementMatrixOutput(GENERAL_OUTPUT_TYPE,linearMatrix,err,error,*999)
        ENDDO !matrixIdx
      ENDIF
      NULLIFY(nonlinearMatrices)
      CALL EquationsMatricesVector_NonlinearMatricesExists(vectorMatrices,nonlinearMatrices,err,error,*999)
      IF(ASSOCIATED(nonlinearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Residual vectors:",err,error,*999)
        CALL EquationsMatricesNonlinear_NumberOfResidualsGet(nonlinearMatrices,numberOfResiduals,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of residuals = ",numberOfResiduals,err,error,*999)
        DO residualIdx=1,numberOfResiduals
          NULLIFY(residualVector)
          CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Residual vector : ",residualIdx,err,error,*999)
          CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
          IF(updateVector) CALL EquationsMatricesResidual_NodalVectorOutput(GENERAL_OUTPUT_TYPE,residualVector,err,error,*999)
        ENDDO !residualIdx
      ENDIF
      NULLIFY(sourceVectors)
      CALL EquationsMatricesVector_SourceVectorsExists(vectorMatrices,sourceVectors,err,error,*999)
      IF(ASSOCIATED(sourceVectors)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Source vectors:",err,error,*999)
        CALL EquationsMatricesSources_NumberOfSourcesGet(sourceVectors,numberOfSources,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of source vectors = ",numberOfSources,err,error,*999)
        DO sourceIdx=1,numberOfSources
          NULLIFY(sourceVector)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Source vector : ",sourceIdx,err,error,*999)
          CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
          IF(updateVector) CALL EquationsMatricesSource_NodalVectorOutput(GENERAL_OUTPUT_TYPE,sourceVector,err,error,*999)
        ENDDO !sourceIdx
      ENDIF
      NULLIFY(rhsVector)
      CALL EquationsMatricesVector_RHSVectorExists(vectorMatrices,rhsVector,err,error,*999)
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"RHS vector :",err,error,*999)
        CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",updateVector,err,error,*999)
        IF(updateVector) CALL EquationsMatricesRHS_NodalVectorOutput(GENERAL_OUTPUT_TYPE,rhsVector,err,error,*999)
      ENDIF
    ENDIF
       
    EXITS("EquationsSet_NodalResidualEvaluate")
    RETURN
999 ERRORSEXITS("EquationsSet_NodalResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_NodalResidualEvaluate


  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for an static equations set using the finite nodal method
  SUBROUTINE EquationsSet_JacobianEvaluateStaticNodal(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,ghostFinish,internalFinish,internalStart,nodeIdx,nodeNumber,numberOfTimes,outputType
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: nodalMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
  
    ENTERS("EquationsSet_JacobianEvaluateStaticNodal",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(nodalMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodalMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(nodalMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(nodalMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(nodalMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(nodalMapping,ghostFinish,err,error,*999)
    
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the nodes
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    numberOfTimes=0
    !Loop over the internal nodes
    DO nodeIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(nodalMapping,nodeIdx,nodeNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_JacobianNodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      nodeUserElapsed=userElapsed
      nodeSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost nodes
    DO nodeIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(nodalMapping,nodeIdx,nodeNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_JacobianNodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average node equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations Jacobian if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) &
      & CALL EquationsMatricesVector_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
       
    EXITS("EquationsSet_JacobianEvaluateStaticNodal")
    RETURN
999 ERRORSEXITS("EquationsSet_JacobianEvaluateStaticNodal",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_JacobianEvaluateStaticNodal

  !
  !================================================================================================================================
  !

  !>Evaluates the residual for an static equations set using the nodal method
  SUBROUTINE EquationsSet_ResidualEvaluateStaticNodal(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: boundaryStart,ghostFinish,internalFinish,internalStart,nodeIdx,nodeNumber,numberOfTimes,outputType
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: nodalMapping
    TYPE(DomainMappingsType), POINTER :: domainMappings
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField
 
    ENTERS("EquationsSet_ResidualEvaluateStaticNodal",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_OutputTypeGet(equations,outputType,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(domain)
    CALL Decomposition_DomainGet(decomposition,0,domain,err,error,*999)
    NULLIFY(domainMappings)
    CALL Domain_DomainMappingsGet(domain,domainMappings,err,error,*999)
    NULLIFY(nodalMapping)
    CALL DomainMappings_NodesMappingGet(domainMappings,nodalMapping,err,error,*999)
    CALL DomainMapping_InternalStartGet(nodalMapping,internalStart,err,error,*999)
    CALL DomainMapping_InternalFinishGet(nodalMapping,internalFinish,err,error,*999)
    CALL DomainMapping_BoundaryStartGet(nodalMapping,boundaryStart,err,error,*999)
    CALL DomainMapping_GhostFinishGet(nodalMapping,ghostFinish,err,error,*999)
    
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
      nodeUserElapsed=0.0_SP
      nodeSystemElapsed=0.0_SP
    ENDIF
    numberOfTimes=0
    !Loop over the internal nodes
    DO nodeIdx=internalStart,internalFinish
      CALL DomainMapping_NumberGet(nodalMapping,nodeIdx,nodeNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      nodeUserElapsed=userElapsed
      nodeSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost nodes
    DO nodeIdx=boundaryStart,ghostFinish
      CALL DomainMapping_NumberGet(nodalMapping,nodeIdx,nodeNumber,err,error,*999)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average node equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations residual vector if required
    IF(outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatricesVector_Output(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_ResidualEvaluateStaticNodal")
    RETURN
999 ERRORSEXITS("EquationsSet_ResidualEvaluateStaticNodal",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_ResidualEvaluateStaticNodal

  !
  !================================================================================================================================
  !

END MODULE EquationsSetRoutines
