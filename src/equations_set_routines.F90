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

!> This module handles all equations set routines.
MODULE EquationsSetRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BIOELECTRIC_ROUTINES
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
  USE MONODOMAIN_EQUATIONS_ROUTINES
#ifndef NOMPIMOD
  USE MPI
#endif
  USE MULTI_PHYSICS_ROUTINES
  USE ProfilingRoutines
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
  
  PUBLIC EquationsSet_Backsubstitute,EquationsSet_NonlinearRHSUpdate
  
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
      equationsSetSetupInfo%fieldUserNumber=analyticField%userNumber
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
    INTEGER(INTG) :: dummyErr
    TYPE(DecompositionType), POINTER :: analyticDecomposition,geometricDecomposition
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(FieldType), POINTER :: field,geometricField
    TYPE(RegionType), POINTER :: region,analyticFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_AnalyticCreateStart",err,error,*998)

    CALL EquationsSet_AssertAnalyticNotCreated(equationsSet,err,error,*998)

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_AssertIsFinished(analyticField,err,error,*998)
      !Check the user numbers match
      IF(analyticFieldUserNumber/=analyticField%userNumber) THEN
        localError="The specified analytic field user number of "// &
          & TRIM(NumberToVString(analyticFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified analytic field of "// &
          & TRIM(NumberToVString(analyticField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      NULLIFY(analyticFieldRegion)
      CALL Field_RegionGet(analyticField,analyticFieldRegion,err,error,*998)
      !Check the field is defined on the same region as the equations set
      IF(analyticFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified analytic field has been created on region number "// &
          & TRIM(NumberToVString(analyticFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
      !Check the specified analytic field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(analyticDecomposition)
      CALL Field_DecompositionGet(analyticField,analyticDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,analyticDecomposition)) THEN
        CALL FlagError("The specified analytic field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*998)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(analyticFieldUserNumber,region,field,err,error,*998)
      IF(ASSOCIATED(field)) THEN
        localError="The specified analytic field user number of "// &
          & TRIM(NumberToVString(analyticFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*998)
      ENDIF
    ENDIF
    !Initialise the equations set analytic
    CALL EquationsSet_AnalyticInitialise(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(analyticField)) equationsSet%analytic%analyticFieldAutoCreated=.TRUE.
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
      & globalDerivIndex,nodeIdx,numberOfAnalyticComponents,numberOfDimensions,variableIdx,variableType,versionIdx
    REAL(DP) :: analyticTime,normal(3),position(3),tangents(3,3),VALUE
    REAL(DP) :: analyticDummyValues(1)=0.0_DP
    REAL(DP) :: materialsDummyValues(1)=0.0_DP
    LOGICAL :: reverseNormal=.FALSE.
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FieldInterpolationParametersPtrType), POINTER :: analyticInterpParameters(:),geometricInterpParameters(:), &
      & materialsInterpParameters(:)
    TYPE(FieldInterpolatedPointPtrType), POINTER :: analyticInterpPoint(:),geometricInterpPoint(:), &
      & materialsInterpPoint(:)
    TYPE(FieldInterpolatedPointMetricsPtrType), POINTER :: geometricInterpPointMetrics(:)
    TYPE(FieldPhysicalPointPtrType), POINTER :: analyticPhysicalPoint(:),materialsPhysicalPoint(:)
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_AnalyticEvaluate",err,error,*999)

    CALL EquationsSet_AssertAnalyticIsFinished(equationsSet,err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
    CALL Field_InterpolationParametersInitialise(geometricField,geometricInterpParameters,err,error,*999)
    CALL Field_InterpolatedPointsInitialise(geometricInterpParameters,geometricInterpPoint,err,error,*999)
    CALL Field_InterpolatedPointsMetricsInitialise(geometricInterpPoint,geometricInterpPointMetrics, &
      & err,error,*999)
    NULLIFY(analyticField)
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_NumberOfComponentsGet(analyticField,FIELD_U_VARIABLE_TYPE,numberOfAnalyticComponents, &
        & err,error,*999)
      CALL Field_InterpolationParametersInitialise(analyticField,analyticInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(analyticInterpParameters,analyticInterpPoint,err,error,*999)
      CALL Field_PhysicalPointsInitialise(analyticInterpPoint,geometricInterpPoint,analyticPhysicalPoint, &
        & err,error,*999)
    ENDIF
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    IF(ASSOCIATED(materialsField)) THEN
      CALL Field_NumberOfComponentsGet(materialsField,FIELD_U_VARIABLE_TYPE,numberOfAnalyticComponents, &
        & err,error,*999)
      CALL Field_InterpolationParametersInitialise(materialsField,materialsInterpParameters,err,error,*999)
      CALL Field_InterpolatedPointsInitialise(materialsInterpParameters,materialsInterpPoint,err,error,*999)
      CALL Field_PhysicalPointsInitialise(materialsInterpPoint,geometricInterpPoint,materialsPhysicalPoint, &
        & err,error,*999)
    ENDIF
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    CALL EquationsSet_AnalyticTimeGet(equationsSet,analyticTime,err,error,*999)
    DO variableIdx=1,dependentField%numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      DO componentIdx=1,dependentVariable%numberOfComponents
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
          DO elementIdx=1,domainElements%numberOfElements
            NULLIFY(basis)
            CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(ASSOCIATED(analyticField)) THEN
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
                & analyticInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            IF(ASSOCIATED(materialsField)) THEN
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
                & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            CALL Field_InterpolateXi(FIRST_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE, &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_PositionNormalTangentsCalculateIntPtMetric( &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,reverseNormal, &
              & position,normal,tangents,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolateXi(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
                & analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolateXi(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
!! \todo Maybe do this with optional arguments?
            IF(ASSOCIATED(analyticField)) THEN
              IF(ASSOCIATED(materialsField)) THEN
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivIndex,componentIdx, &
                  & analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                  & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                  & VALUE,err,error,*999)
              ELSE
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivIndex,componentIdx, &
                  & analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                  & materialsDummyValues,VALUE,err,error,*999)
              ENDIF
            ELSE
              IF(ASSOCIATED(materialsField)) THEN
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                  & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                  & VALUE,err,error,*999)
              ELSE
                CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                  & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                  & materialsDummyValues,VALUE,err,error,*999)
              ENDIF
            ENDIF
            CALL Field_ParameterSetUpdateLocalElement(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
              & elementIdx,componentIdx,VALUE,err,error,*999)
          ENDDO !elementIdx
        CASE(FIELD_NODE_BASED_INTERPOLATION)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          DO nodeIdx=1,domainNodes%numberOfNodes
            CALL Field_PositionNormalTangentsCalculateNode(dependentField,variableType,componentIdx,nodeIdx, &
              & position,normal,tangents,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolateFieldNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE, &
                & analyticField,FIELD_U_VARIABLE_TYPE,componentIdx,nodeIdx,analyticPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr, &
                & err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolateFieldNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE, &
              & materialsField,FIELD_U_VARIABLE_TYPE,componentIdx,nodeIdx,materialsPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr, &
              & err,error,*999)
            !Loop over the derivatives
            DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%numberOfDerivatives                                
              globalDerivIndex=domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%globalDerivativeIndex
!! \todo Maybe do this with optional arguments?
              IF(ASSOCIATED(analyticField)) THEN
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionTYpe,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx, &
                    & analyticPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES, &
                    & materialsPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES,VALUE,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx, &
                    & analyticPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES, &
                    & materialsDummyValues,VALUE,err,error,*999)
                ENDIF
              ELSE
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                    & materialsPhysicalPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES,VALUE,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                    & materialsDummyValues,VALUE,err,error,*999)
                ENDIF
              ENDIF
              !Loop over the versions
              DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                CALL Field_ParameterSetUpdateLocalNode(dependentField,variableType, &
                  & FIELD_ANALYTIC_VALUES_SET_TYPE,versionIdx,derivativeIdx,nodeIdx, &
                  & componentIdx,VALUE,err,error,*999)
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
          NULLIFY(domainElements)
          CALL DomainTopology_DomainElementsGet(domainTopology,domainElements,err,error,*999)
          !Loop over the local elements excluding the ghosts
          DO elementIdx=1,domainElements%numberOfElements
            NULLIFY(basis)
            CALL DomainElements_BasisGet(domainElements,elementIdx,basis,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(ASSOCIATED(analyticField)) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & analyticInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(ASSOCIATED(materialsField)) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementIdx, &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Loop over the Gauss points in the element
            NULLIFY(quadratureScheme)
            CALL Basis_QuadratureSchemeGet(basis,BASIS_DEFAULT_QUADRATURE_SCHEME,quadratureScheme,err,error,*999)
            DO gaussPointIdx=1,quadratureScheme%numberOfGauss
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE, &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_PositionNormalTangentsCalculateIntPtMetric(geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr, &
                & reverseNormal,position,normal,tangents,err,error,*999)
              IF(ASSOCIATED(analyticField)) CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME, &
                & gaussPointIdx,analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              IF(ASSOCIATED(materialsField)) CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME, &
                & gaussPointIdx,materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
!! \todo Maybe do this with optional arguments?
              IF(ASSOCIATED(analyticField)) THEN
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx, &
                    & analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                    & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                    & VALUE,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx, &
                    & analyticInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                    & materialsDummyValues,VALUE,err,error,*999)
                ENDIF
              ELSE
                IF(ASSOCIATED(materialsField)) THEN
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                    & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                    & VALUE,err,error,*999)
                ELSE
                  CALL EquationsSet_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,position,tangents,normal, &
                    & analyticTime,variableType,globalDerivIndex,componentIdx,analyticDummyValues, &
                    & materialsDummyValues,VALUE,err,error,*999)
                ENDIF
              ENDIF
              CALL Field_ParameterSetUpdateLocalGaussPoint(dependentField,variableType, &
                & FIELD_ANALYTIC_VALUES_SET_TYPE,gaussPointIdx,elementIdx,componentIdx, &
                & VALUE,err,error,*999)
            ENDDO !gaussPointIdx
          ENDDO !elementIdx
        CASE DEFAULT
          localError="The interpolation type of "//TRIM(NumberToVString(componentInterpolationType,"*",err,error))// &
            & " for component "//TRIM(NumberToVString(componentIdx,"*",err,error))//" of variable type "// &
            & TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !componentIdx
      CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    IF(ASSOCIATED(materialsField)) THEN
      CALL Field_PhysicalPointsFinalise(materialsPhysicalPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(materialsInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersFinalise(materialsInterpParameters,err,error,*999)
    ENDIF
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_PhysicalPointsFinalise(analyticPhysicalPoint,err,error,*999)
      CALL Field_InterpolatedPointsFinalise(analyticInterpPoint,err,error,*999)
      CALL Field_InterpolationParametersFinalise(analyticInterpParameters,err,error,*999)
    ENDIF
    CALL Field_InterpolatedPointsMetricsFinalise(geometricInterpPointMetrics,err,error,*999)
    CALL Field_InterpolatedPointsFinalise(geometricInterpPoint,err,error,*999)
    CALL Field_InterpolationParametersFinalise(geometricInterpParameters,err,error,*999)
           
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
    & variableType,globalDerivative,componentNumber,analyticParameters,materialsParameters,value,err,error,*)

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
    REAL(DP), INTENT(OUT) :: value !<On return, the analtyic function value.
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
    
    SELECT CASE(equationsSet%SPECIFICATION(1))
    CASE(EQUATIONS_SET_ELASTICITY_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      IF(SIZE(equationsSet%specification,1)<2) CALL FlagError("Equations set specification must have at least two "// &
        & "entries for a classical field equations set.",err,error,*999)
      CALL ClassicalField_AnalyticFunctionsEvaluate(equationsSet,equationsSet%specification(2), &
        & analyticFunctionType,position,tangents,normal,time,variableType,globalDerivative, &
        & componentNumber,analyticParameters,materialsParameters,value,err,error,*999)
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
      CALL FLAG_ERROR(localError,err,error,*999)
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
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_Assemble",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    CALL Equations_AssertIsFinished(equations,err,error,*999)
    
    IF(equationsSet%outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set assemble: ",equationsSet%label,err,error,*999)
    ENDIF
    
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      SELECT CASE(equations%linearity)
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
        localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_QUASISTATIC)
      SELECT CASE(equations%linearity)
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
        localError="The equations linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
      SELECT CASE(equations%linearity)
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
        localError="The equations set linearity of "//TRIM(NumberToVString(equations%linearity,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_TIME_STEPPING)
      CALL FlagError("Time stepping equations are not assembled.",err,error,*999)
    CASE DEFAULT
      localError="The equations time dependence type of "//TRIM(NumberToVString(equations%timeDependence,"*",err,error))// &
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
 
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatricesVector_ElementInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: equationsColumnIdx,equationsColumnNumber,equationsMatrixIdx,equationsRowNumber, &
      & equationsStorageType,rhsBoundaryCondition,rhsGlobalDOF,rhsVariableDOF,rhsVariableType,variableDOF,variableType
    INTEGER(INTG), POINTER :: columnIndices(:),rowIndices(:)
    REAL(DP) :: dependentValue,matrixValue,rhsValue,sourceValue
    REAL(DP), POINTER :: dependentParameters(:),equationsMatrixData(:),sourceVectorData(:)
    TYPE(BoundaryConditionVariableType), POINTER :: rhsBoundaryConditions
    TYPE(DomainMappingType), POINTER :: columnDomainMapping,rhsDomainMapping
    TYPE(DistributedMatrixType), POINTER :: equationsDistributedMatrix
    TYPE(DistributedVectorType), POINTER :: sourceDistributedVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable,rhsVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_Backsubstitute",err,error,*999)

    CALL EquationsSet_AssertIsFinished(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)

    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
 
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)

    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingExists(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingExists(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingExists(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)

    dynamicMatrices=>vectorMatrices%dynamicMatrices
    IF(ASSOCIATED(dynamicMatrices)) THEN
      !CALL FlagError("Not implemented.",err,error,*999)
    ELSE
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(rhsMapping)
      CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
      NULLIFY(sourcesMapping)
      CALL EquationsMappingVector_SourcesMappingGet(vectorMapping,sourcesMapping,err,error,*999)
      
      IF(ASSOCIATED(sourcesMapping)) THEN
        DO sourceIdx=1,
        NULLIFY(sourceVector)
        CALL EquationsMatricesVector_SourceVectorGet(vectorMatrices,sourceVector,err,error,*999)
        sourceDistributedVector=>sourceVector%vector
        IF(.NOT.ASSOCIATED(sourceDistributedVector)) &
          & CALL FlagError("Source distributed vector is not associated.",err,error,*999)
        CALL DistributedVector_DataGet(sourceDistributedVector,sourceVectorData,err,error,*999)
      ENDIF
      IF(ASSOCIATED(rhsMapping)) THEN
        NULLIFY(rhsVariable)
        CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
        rhsVariableType=rhsVariable%variableType
        NULLIFY(rhsDomainMapping)
        CALL FieldVariable_DomainMappingGet(rhsVariable,rhsDomainMapping,err,error,*999)
        CALL BoundaryConditions_VariableExists(boundaryConditions,rhsVariable,rhsBoundaryConditions,err,error,*999)
        IF(ASSOCIATED(rhsBoundaryConditions)) THEN
          !Loop over the equations matrices
          DO equationsMatrixIdx=1,linearMatrices%numberOfLinearMatrices
            NULLIFY(dependentVariable)
            CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,equationsMatrixIdx,dependentVariable,err,error,*999)
            variableType=dependentVariable%variableType
            !Get the dependent field variable parameters
            CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_VALUES_SET_TYPE,dependentParameters,err,error,*999)
            NULLIFY(equationsMatrix)
            CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,equationsMatrixIdx,equationsMatrix,err,error,*999)
            columnDomainMapping=>linearMapping%equationsMatrixToVarMaps(equationsMatrixIdx)%columnDOFSMapping
            IF(ASSOCIATED(columnDomainMapping)) &
              & CALL FlagError("Equations column domain mapping is not associated.",err,error,*999)
            equationsDistributedMatrix=>equationsMatrix%matrix
            IF(.NOT.ASSOCIATED(equationsDistributedMatrix)) &
              & CALL FlagError("Equations matrix distributed matrix is not associated.",err,error,*999)
            CALL DistributedMatrix_StorageTypeGet(equationsDistributedMatrix,equationsStorageType,err,error,*999)
            CALL DistributedMatrix_DataGet(equationsDistributedMatrix,equationsMatrixData,err,error,*999)
            SELECT CASE(equationsStorageType)
            CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
              !Loop over the non ghosted rows in the equations set
              DO equationsRowNumber=1,vectorMapping%numberOfRows
                rhsValue=0.0_DP
                rhsVariableDOF=rhsMapping%equationsRowToRHSDOFMap(equationsRowNumber)
                rhsGlobalDOF=rhsDomainMapping%localToGlobalMap(rhsVariableDOF)
                rhsBoundaryCondition=rhsBoundaryConditions%DOFTypes(rhsGlobalDOF)
                !For free RHS DOFs, set the right hand side field values by multiplying the
                !row by the dependent variable value
                SELECT CASE(rhsBoundaryCondition)
                CASE(BOUNDARY_CONDITION_DOF_FREE)
                  !Back substitute
                  !Loop over the local columns of the equations matrix
                  DO equationsColumnIdx=1,columnDomainMapping%totalNumberOfLocal
                    equationsColumnNumber=columnDomainMapping%localToGlobalMap(equationsColumnIdx)
                    variableDOF=equationsColumnIdx
                    matrixValue=equationsMatrixData(equationsRowNumber+ &
                      & (equationsColumnNumber-1)*vectorMatrices%totalNumberOfRows)
                    dependentValue=dependentParameters(variableDOF)
                    rhsValue=rhsValue+matrixValue*dependentValue
                  ENDDO !equationsColumnIdx
                CASE(BOUNDARY_CONDITION_DOF_FIXED)
                  !Do nothing
                CASE(BOUNDARY_CONDITION_DOF_MIXED)
                  !Robin or is it Cauchy??? boundary conditions
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The RHS variable boundary condition of "// &
                    & TRIM(NumberToVString(rhsBoundaryCondition,"*",err,error))// &
                    & " for RHS variable dof number "// &
                    & TRIM(NumberToVString(rhsVariableDOF,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                IF(ASSOCIATED(sourceMapping)) THEN
                  sourceValue=sourceVectorData(equationsRowNumber)
                  rhsValue=rhsValue-sourceValue
                ENDIF
                CALL Field_ParameterSetUpdateLocalDOF(dependentField,rhsVariableType, &
                  & FIELD_VALUES_SET_TYPE,rhsVariableDOF,rhsValue,err,error,*999)
              ENDDO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
              NULLIFY(rowIndices)
              NULLIFY(columnIndices)
              CALL DistributedMatrix_StorageLocationsGet(equationsDistributedMatrix,rowIndices,columnIndices,err,error,*999)
              !Loop over the non-ghosted rows in the equations set
              DO equationsRowNumber=1,vectorMapping%numberOfRows
                rhsValue=0.0_DP
                rhsVariableDOF=rhsMapping%equationsRowToRHSDOFMap(equationsRowNumber)
                rhsGlobalDOF=rhsDomainMapping%localToGlobalMap(rhsVariableDOF)
                rhsBoundaryCondition=rhsBoundaryConditions%DOFTypes(rhsGlobalDOF)
                SELECT CASE(rhsBoundaryCondition)
                CASE(BOUNDARY_CONDITION_DOF_FREE)
                  !Back substitute
                  !Loop over the local columns of the equations matrix
                  DO equationsColumnIdx=rowIndices(equationsRowNumber),rowIndices(equationsRowNumber+1)-1
                    equationsColumnNumber=columnIndices(equationsColumnIdx)
                    variableDOF=columnDomainMapping%globalToLocalMap(equationsColumnNumber)%localNumber(1)
                    matrixValue=equationsMatrixData(equationsColumnIdx)
                    dependentValue=dependentParameters(variableDOF)
                    rhsValue=rhsValue+matrixValue*dependentValue
                  ENDDO !equationsColumnIdx
                CASE(BOUNDARY_CONDITION_DOF_FIXED)
                  !Do nothing
                CASE(BOUNDARY_CONDITION_DOF_MIXED)
                  !Robin or is it Cauchy??? boundary conditions
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global boundary condition of "// &
                    & TRIM(NumberToVString(rhsBoundaryCondition,"*",err,error))// &
                    & " for RHS variable dof number "// &
                    & TRIM(NumberToVString(rhsVariableDOF,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                IF(ASSOCIATED(sourceMapping)) THEN
                  sourceValue=sourceVectorData(equationsRowNumber)
                  rhsValue=rhsValue-sourceValue
                ENDIF
                CALL Field_ParameterSetUpdateLocalDOF(dependentField,rhsVariableType, &
                  & FIELD_VALUES_SET_TYPE,rhsVariableDOF,rhsValue,err,error,*999)
              ENDDO !equationsRowNumber
            CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The matrix storage type of "// &
                & TRIM(NumberToVString(equationsStorageType,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL DistributedMatrix_DataRestore(equationsDistributedMatrix,equationsMatrixData,err,error,*999)
            !Restore the dependent field variable parameters
            CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
              & dependentParameters,err,error,*999)
          ENDDO !equationsMatrixIdx
          
          !Start the update of the field parameters
          CALL Field_ParameterSetUpdateStart(dependentField,rhsVariableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Finish the update of the field parameters
          CALL Field_ParameterSetUpdateFinish(dependentField,rhsVariableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(ASSOCIATED(sourceMapping)) &
            & CALL DistributedVector_DataRestore(sourceDistributedVector,sourceVectorData,err,error,*999)
        ENDIF
      ENDIF
    ENDIF
          
    EXITS("EquationsSet_Backsubstitute")
    RETURN
999 ERRORSEXITS("EquationsSet_Backsubstitute",err,error)
    RETURN 1
   
  END SUBROUTINE EquationsSet_Backsubstitute
  
  !
  !================================================================================================================================
  !

  !>Updates the right hand side variable from the equations residual vector
  SUBROUTINE EquationsSet_NonlinearRHSUpdate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<Boundary conditions to use for the RHS update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variableDOF,rowIdx,variableType,rhsGlobalDOF,rhsBoundaryCondition,equationsMatrixIdx
    REAL(DP) :: VALUE
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(DistributedVectorType), POINTER :: residualVector
    TYPE(FieldType), POINTER :: rhsField
    TYPE(FieldVariableType), POINTER :: rhsVariable,residualVariable
    TYPE(BoundaryConditionVariableType), POINTER :: rhsBoundaryConditions
    TYPE(DomainMappingType), POINTER :: rhsDomainMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_NonlinearRHSUpdate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)
      
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingGet(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(rhsVariable)
    CALL EquationsMappingRHS_RHSVariableGet(rhsMapping,rhsVariable,err,error,*999)
    NULLIFY(rhsField)
    CALL FieldVariable_FieldGet(rhsVariable,rhsField,err,error,*999)
    variableType=rhsVariable%variableType
    NULLIFY(rhsDomainMapping)
    CALL FieldVariable_DomainMappingGet(rhsVariable,rhsDomainMapping,err,error,*999)
    NULLIFY(rhsBoundaryConditions)
    CALL BoundaryConditions_VariableGet(boundaryConditions,rhsVariable,rhsBoundaryConditions,err,error,*999)
    !Get the equations residual vector
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    residualVector=>nonlinearMatrices%residual
    IF(.NOT.ASSOCIATED(residualVector)) CALL FlagError("Residual vector is not associated.",err,error,*999)
    DO equationsMatrixIdx=1,nonlinearMapping%numberOfResidualVariables
      NULLIFY(residualVariable)
      CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,1,equationsMatrixIdx,residualVariable,err,error,*999)
      DO rowIdx=1,vectorMapping%numberOfRows
        variableDOF=rhsMapping%equationsRowToRHSDOFMap(rowIdx)
        rhsGlobalDOF=rhsDomainMapping%localToGlobalMap(variableDOF)
        rhsBoundaryCondition=rhsBoundaryConditions%DOFTypes(rhsGlobalDOF)
        SELECT CASE(rhsBoundaryCondition)
        CASE(BOUNDARY_CONDITION_DOF_FREE)
          !Add residual to field value
          CALL DistributedVector_ValuesGet(residualVector,rowIdx,VALUE,err,error,*999)
          CALL Field_ParameterSetUpdateLocalDOF(rhsField,variableType,FIELD_VALUES_SET_TYPE,variableDOF,VALUE,err,error,*999)
        CASE(BOUNDARY_CONDITION_DOF_FIXED)
          !Do nothing
        CASE(BOUNDARY_CONDITION_DOF_MIXED)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The RHS variable boundary condition of "// &
            & TRIM(NumberToVString(rhsBoundaryCondition,"*",err,error))// &
            & " for RHS variable dof number "// &
            & TRIM(NumberToVString(variableDOF,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !rowIdx
    ENDDO !equationsMatrixIdx
    CALL Field_ParameterSetUpdateStart(rhsField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
    CALL Field_ParameterSetUpdateFinish(rhsField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)

    EXITS("EquationsSet_NonlinearRHSUpdate")
    RETURN
999 ERRORSEXITS("EquationsSet_NonlinearRHSUpdate",err,error)
    RETURN 1

  END SUBROUTINE EquationsSet_NonlinearRHSUpdate

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
    INTEGER(INTG) :: dummyErr,equationsSetIdx
    TYPE(DecompositionType), POINTER :: equationsSetFieldDecomposition,geometricFibreFieldDecomposition
    TYPE(EquationsSetType), POINTER :: newEquationsSet
    TYPE(EquationsSetPtrType), POINTER :: newEquationsSets(:)
    TYPE(EquationsSetSetupType) equationsSetSetupInfo
    TYPE(RegionType), POINTER :: geometricFibreFieldRegion,equationsSetFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(EquationsSetEquationsSetFieldType), POINTER :: equationsEquationsSetField
    TYPE(FieldType), POINTER :: FIELD

    NULLIFY(newEquationsSet)
    NULLIFY(newEquationsSets)
    NULLIFY(equationsEquationsSetField)

    ENTERS("EquationsSet_CreateStart",err,error,*997)

    IF(.NOT.ASSOCIATED(equationsSetRegion)) CALL FlagError("Region is not associated.",err,error,*997)
    IF(.NOT.ASSOCIATED(equationsSetRegion%equationsSets)) THEN
      localError="The equations sets on region number "//TRIM(NumberToVString(equationsSetRegion%userNumber,"*",err,error))// &
        & " are not associated."
      CALL FlagError(localError,err,error,*997)
    ENDIF
    CALL Field_AssertIsFinished(geometricFibreField,err,error,*999)
    CALL EquationsSet_UserNumberFind(userNumber,equationsSetRegion,newEquationsSet,err,error,*997)
    IF(ASSOCIATED(newEquationsSet)) THEN
      localError="Equations set user number "//TRIM(NumberToVString(userNumber,"*",err,error))// &
        & " has already been created on region number "//TRIM(NumberToVString(equationsSetRegion%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*997)
    ENDIF
    NULLIFY(newEquationsSet)
    IF(geometricFibreField%TYPE/=FIELD_GEOMETRIC_TYPE.AND.geometricFibreField%TYPE/=FIELD_FIBRE_TYPE) &
      & CALL FlagError("The specified geometric field is not a geometric or fibre field.",err,error,*997)
    NULLIFY(geometricFibreFieldRegion)
    CALL Field_RegionGet(geometricFibreField,geometricFibreFieldRegion,err,error,*999)
    IF(geometricFibreFieldRegion%userNumber/=equationsSetRegion%userNumber) THEN
      localError="The geometric field region and the specified region do not match. "// &
        & "The geometric field was created on region number "// &
        & TRIM(NumberToVString(geometricFibreFieldRegion%userNumber,"*",err,error))// &
        & " and the specified equation set region number is "// &
        & TRIM(NumberToVString(equationsSetRegion%userNumber,"*",err,error))//"."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    IF(ASSOCIATED(equationsSetField)) THEN
      !Check the equations set field has been finished
      CALL Field_AssertIsFinished(equationsSetField,err,error,*999)
      !Check the user numbers match
      IF(equationsSetFieldUserNumber/=equationsSetField%userNumber) THEN
        localError="The specified equations set field user number of "// &
          & TRIM(NumberToVString(equationsSetFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified equations set field of "// &
          & TRIM(NumberToVString(equationsSetField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(equationsSetFieldRegion)
      CALL Field_RegionGet(equationsSetField,equationsSetFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the equations set
      IF(equationsSetFieldRegion%userNumber/=equationsSetRegion%userNumber) THEN
        localError="Invalid region setup. The specified equations set field was created on region no. "// &
          & TRIM(NumberToVString(equationsSetFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(equationsSetRegion%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified equations set field has the same decomposition as the geometric field
      NULLIFY(geometricFibreFieldDecomposition)
      CALL Field_DecompositionGet(geometricFibreField,geometricFibreFieldDecomposition,err,error,*999)
      NULLIFY(equationsSetFieldDecomposition)
      CALL Field_DecompositionGet(equationsSetField,equationsSetFieldDecomposition,err,error,*999)
      IF(.NOT.ASSOCIATED(geometricFibreFieldDecomposition,equationsSetFieldDecomposition)) THEN
        CALL FlagError("The specified equations set field does not have the same decomposition "// &
          & "as the geometric field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(equationsSetFieldUserNumber,equationsSetRegion,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified equations set field user number of "// &
          & TRIM(NumberToVString(equationsSetFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(equationsSetRegion%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    
    !Initalise equations set
    CALL EquationsSet_Initialise(newEquationsSet,err,error,*999)
    !Set default equations set values
    newEquationsSet%userNumber=userNumber
    newEquationsSet%globalNumber=equationsSetRegion%equationsSets%numberOfEquationsSets+1
    newEquationsSet%equationsSets=>equationsSetRegion%equationsSets
    newEquationsSet%label="Equations Set "//TRIM(NumberToVString(userNumber,"*",err,error))
    newEquationsSet%REGION=>equationsSetRegion
    !Set the equations set class, type and subtype
    CALL EquationsSet_SpecificationSet(newEquationsSet,equationsSetSpecification,err,error,*999)
    newEquationsSet%equationsSetFinished=.FALSE.
    !Initialise the setup
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
    IF(geometricFibreField%type==FIELD_GEOMETRIC_TYPE) THEN
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
    equationsEquationsSetField=>equationsSet%equationsSetField
    !\todo check pointer setup
    IF(equationsEquationsSetField%equationsSetFieldAutoCreated) THEN
      equationsSetField=>equationsSet%equationsSetField%equationsSetFieldField
    ELSE
      equationsSet%equationsSetField%equationsSetFieldField=>equationsSetField
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
      CALL EqutionsSet_EquationsSetFieldFinalise(equationsSet%equationsSetField,err,error,*999)
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
    INTEGER(INTG) :: matrixIdx
    TYPE(ElementMatrixType), POINTER :: elementMatrix
    TYPE(ElementVectorType), POINTER :: elementVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_FiniteElementCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)    
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
      IF(SIZE(equationsSet%specification,1)<2) &
        & CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation class.", &
        & err,error,*999)
      IF(equationsSet%specification(2) == EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
        CALL Monodomain_FiniteElementCalculate(equationsSet,elementNumber,err,error,*999)
      ELSE
        CALL BIOELECTRIC_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
      END IF
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
    CASE DEFAULT
      localError="The first equations set specification of "// &
        & TRIM(NumberToVString(equationsSet%SPECIFICATION(1),"*",err,error))//" is not valid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMatrices)
      CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element stiffness matrices:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      dynamicMatrices=>vectorMatrices%dynamicMatrices
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Dynamic matrices:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",dynamicMatrices%numberOfDynamicMatrices, &
          & err,error,*999)
        DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",dynamicMatrices%matrices(matrixIdx)%ptr%updateMatrix, &
            & err,error,*999)
          IF(dynamicMatrices%matrices(matrixIdx)%ptr%updateMatrix) THEN
            elementMatrix=>dynamicMatrices%matrices(matrixIdx)%ptr%elementMatrix
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",elementMatrix%numberOfColumns,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",elementMatrix%maxNumberOfColumns, &
              & err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix% &
              & columnDOFS,'("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,1,1,elementMatrix% &
              & numberOfColumns,8,8,elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)','     :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      linearMatrices=>vectorMatrices%linearMatrices
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",linearMatrices%numberOfLinearMatrices, &
          & err,error,*999)
        DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",linearMatrices%matrices(matrixIdx)%ptr%updateMatrix, &
            & err,error,*999)
          IF(linearMatrices%matrices(matrixIdx)%ptr%updateMatrix) THEN
            elementMatrix=>linearMatrices%matrices(matrixIdx)%ptr%elementMatrix
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",elementMatrix%numberOfColumns,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",elementMatrix%maxNumberOfColumns, &
              & err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix% &
              & columnDOFS,'("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,1,1,elementMatrix% &
              & numberOfColumns,8,8,elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)','     :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      rhsVector=>vectorMatrices%rhsVector
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element RHS vector :",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",rhsVector%updateVector,err,error,*999)
        IF(rhsVector%updateVector) THEN
          elementVector=>rhsVector%elementVector
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementVector%maxNumberOfRows, &
            & err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
            & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
            & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
        ENDIF
      ENDIF
      sourceVectors=>vectorMatrices%sourceVectors
      IF(ASSOCIATED(sourceVectors)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element source vectors:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of source vectors = ",sourceVectors%numberOfSources,err,error,*999)
        DO sourceIdx=1,sourceVectors%numberOfSources
          NULLIFY(sourceVector)
          CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,sourceIdx,sourceVector,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Source vector : ",sourceIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",sourceVector%updateVector,err,error,*999)
          IF(sourceVector%updateVector) THEN
            elementVector=>sourceVector%elementVector
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementVector%maxNumberOfRows, &
              & err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
              & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !sourceIdx
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
    INTEGER(INTG) :: matrixIdx
    TYPE(ElementMatrixType), POINTER :: elementMatrix
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_FiniteElementJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)

    DO residualIdx=1,nonlinearMatrices%numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
      DO matrixIdx=1,residualVector%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
        SELECT CASE(jacobianMatrix%jacobianCalculationType)
        CASE(EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
          ! None of these routines currently support calculating off diagonal terms for coupled problems,
          ! but when one does we will have to pass through the matrixIdx parameter
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
              & TRIM(NumberToVString(equationsSet%SPECIFICATION(1),"*", &
              & err,error))//" is not valid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_JACOBIAN_FINITE_DIFFERENCE_CALCULATED)
          CALL EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,eresidualIdx,matrixIdx,lementNumber,err,error,*999)
        CASE DEFAULT
          localError="The Jacobian calculation type of "// &
            & TRIM(NumberToVString(jacobianMatrix%jacobianCalculationType,"*",err,error))//" is not valid for matrix index "// &
            & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDDO !matrixIdx
    ENDDO !residualIdx
    IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element Jacobian matrix:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Element Jacobian:",err,error,*999)
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Jacobian number = ",matrixIdx,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",nonlinearMatrices%jacobians(matrixIdx)%ptr% &
          & updateJacobian,err,error,*999)
        IF(nonlinearMatrices%jacobians(matrixIdx)%ptr%updateJacobian) THEN
          elementMatrix=>nonlinearMatrices%jacobians(matrixIdx)%ptr%elementJacobian
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",elementMatrix%numberOfColumns, &
            & err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows, &
            & err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",elementMatrix% &
            & maxNumberOfColumns,err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
            & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix% &
            & columnDOFS,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',err,error,*999)
          CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,1,1,elementMatrix% &
            & numberOfColumns,8,8,elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix% &
            & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
            & '(16X,8(X,E13.6))',err,error,*999)
!!TODO: Write out the element residual???
        END IF
      END DO !matrixIdx
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
    INTEGER(INTG), INTENT(IN) :: residualNumber  !<The residual number to calculate when there are coupled problems
    INTEGER(INTG), INTENT(IN) :: jacobianNumber  !<The Jacobian number to calculate when there are coupled problems
    INTEGER(INTG), INTENT(IN) :: elementNumber  !<The element number to calculate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err  !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    TYPE(BasisType), POINTER :: basis
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainElementsType), POINTER :: domainElements
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(DistributedVectorType), POINTER :: parameters
    TYPE(EquationsType), POINTER :: equations
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldParameterSetType), POINTER :: parameterSet
    TYPE(FieldVariableType), POINTER :: rowVariable,columnVariable
    TYPE(ElementVectorType) :: elementVector
    INTEGER(INTG) :: componentIdx,localDOF,version,derivativeIdx,derivative,nodeIdx,node,column
    INTEGER(INTG) :: componentInterpolationType
    INTEGER(INTG) :: numberOfRows
    REAL(DP) :: delta,origDepVar

    ENTERS("EquationsSet_FiniteElementJacobianEvaluateFD",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
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
    
    ! The first residual variable is always the row variable, which is the variable the
    ! residual is calculated for
    NULLIFY(rowVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowVariable,err,error,*999)
    ! For coupled problems this routine will be called multiple times if multiple Jacobians use finite
    ! differencing, so make sure we only calculate the residual vector once, to save time and because
    ! it would otherwise add together
    IF(nonlinearMatrices%elementResidualCalculated/=elementNumber) &
      & CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,residualNumber,elementNumber,err,error,*999)
    ! Make a temporary copy of the unperturbed residuals
    elementVector=residualVector%elementResidual
    ! For coupled nonlinear problems there will be multiple Jacobians
    ! For this equations set, we calculate the residual for the row variable
    ! while pertubing parameters from the column variable.
    ! For non coupled problems these two variables will be the same
    NULLIFY(columnVariable)
    CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,jacobianNumber,residualNumber,columnVariable,err,error,*999)
    NULLIFY(parameterSet)
    CALL FieldVariable_ParameterSetGet(columnVariable,FIELD_VALUES_SET_TYPE,parameterSet,err,error,*999)
    NULLIFY(parameters)
    CALL FieldParameterSet_ParametersGet(parameterSet,parameters,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualNumber,residualVector,err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,jacobianNumber,jacobianMatrix,err,error,*999)
    numberOfRows=jacobianMatrix%elementJacobian%numberOfRows
    IF(numberOfRows/=nonlinearMatrices%elementResidual%numberOfRows) &
      & CALL FlagError("Element matrix number of rows does not match element residual vector size.",err,error,*999)
    ! determine step size
    CALL DistributedVector_L2Norm(parameters,delta,err,error,*999)
    !delta=(1.0_DP+delta)*1.0E-6_DP
    delta=(1.0_DP+delta)*jacobianMatrix%jacobianFiniteDifferenceStepSize
    ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
    ! distributed out, have to use proper field accessing routines..
    ! so let's just loop over component, node/el, derivative
    column=0  ! element jacobian matrix column number
    DO componentIdx=1,columnVariable%numberOfComponents
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
        CALL DomainElements_BasisGet(domainElements,elementNumber,basis,err,error,*999)
        DO nodeIdx=1,basis%numberOfNodes
          node=domainElements%elements(elementNumber)%elementNodes(nodeIdx)
          DO derivativeIdx=1,basis%numberOfDerivatives(nodeIdx)
            derivative=domainElements%elements(elementNumber)%elementDerivatives(derivativeIdx,nodeIdx)
            version=domainElements%elements(elementNumber)%elementVersions(derivativeIdx,nodeIdx)
            CALL FieldVariable_LocalNodeDOFGet(columnVariable,version,derivative,node,componentIdx,localDOF, &
              & err,error,*999)
            ! one-sided finite difference
            CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
            CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
            nonlinearMatrices%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
            CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,residualNumber,elementNumber,err,error,*999)
            CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
            column=column+1
            jacobianMatrix%elementJacobian%matrix(1:numberOfRows,column)= &
              & (nonlinearMatrices%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
        CALL FieldVariable_LocalElementDOFGet(columnVariable,elementNumber,componentIdx,localDOF,err,error,*999)
        ! one-sided finite difference
        CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
        CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
        nonlinearMatrices%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
        CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
        CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
        column=column+1
        jacobianMatrix%elementJacobian%matrix(1:numberOfRows,column)= &
          & (nonlinearMatrices%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
      CASE DEFAULT
        CALL FlagError("Unsupported type of interpolation.",err,error,*999)
      END SELECT
    END DO !componentIdx
    ! put the original residual back in
    nonlinearMatrices%elementResidual=elementVector

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
    INTEGER(INTG) :: matrixIdx
    TYPE(ElementMatrixType), POINTER :: elementMatrix
    TYPE(ElementVectorType), POINTER :: elementVector
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_FiniteElementResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Finite element residual matrices and vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element number = ",elementNumber,err,error,*999)
      linearMatrices=>vectorMatrices%linearMatrices
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",linearMatrices% &
          & numberOfLinearMatrices,err,error,*999)
        DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",linearMatrices%matrices(matrixIdx)%ptr% &
            & updateMatrix,err,error,*999)
          IF(linearMatrices%matrices(matrixIdx)%ptr%updateMatrix) THEN
            elementMatrix=>linearMatrices%matrices(matrixIdx)%ptr%elementMatrix
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",elementMatrix%numberOfColumns, &
              & err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows, &
              & err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",elementMatrix% &
              & maxNumberOfColumns,err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix% &
              & columnDOFS,'("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,1,1,elementMatrix% &
              & numberOfColumns,8,8,elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)','     :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      dynamicMatrices=>vectorMatrices%dynamicMatrices
      IF(ASSOCIATED(dynamicMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Dynamnic matrices:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of element matrices = ",dynamicMatrices% &
          & numberOfDynamicMatrices,err,error,*999)
        DO matrixIdx=1,dynamicMatrices%numberOfDynamicMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Element matrix : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",dynamicMatrices%matrices(matrixIdx)%ptr% &
            & updateMatrix,err,error,*999)
          IF(dynamicMatrices%matrices(matrixIdx)%ptr%updateMatrix) THEN
            elementMatrix=>dynamicMatrices%matrices(matrixIdx)%ptr%elementMatrix
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",elementMatrix%numberOfColumns, &
              & err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementMatrix%maxNumberOfRows, &
              & err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",elementMatrix% &
              & maxNumberOfColumns,err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,8,8,elementMatrix%rowDOFS, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfColumns,8,8,elementMatrix% &
              & columnDOFS,'("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,elementMatrix%numberOfRows,1,1,elementMatrix% &
              & numberOfColumns,8,8,elementMatrix%matrix(1:elementMatrix%numberOfRows,1:elementMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)','     :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Element residual vector:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",nonlinearMatrices%updateResidual,err,error,*999)
      IF(nonlinearMatrices%updateResidual) THEN
        elementVector=>nonlinearMatrices%elementResidual
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementVector%maxNumberOfRows, &
          & err,error,*999)
        CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
          & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
          & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
      ENDIF
      rhsVector=>vectorMatrices%rhsVector
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element RHS vector :",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",rhsVector%updateVector,err,error,*999)
        IF(rhsVector%updateVector) THEN
          elementVector=>rhsVector%elementVector
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementVector%maxNumberOfRows, &
            & err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
            & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
            & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
        ENDIF
      ENDIF
      sourceVector=>vectorMatrices%sourceVector
      IF(ASSOCIATED(sourceVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Element source vector :",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",sourceVector%updateVector,err,error,*999)
        IF(sourceVector%updateVector) THEN
          elementVector=>sourceVector%elementVector
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",elementVector%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",elementVector%maxNumberOfRows, &
            & err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%rowDOFS, &
            & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,elementVector%numberOfRows,8,8,elementVector%vector, &
            & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
        ENDIF
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
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(independentField)) THEN
      CALL Field_AssertIsFinished(independentField,err,error,*999)
      !Check the user numbers match
      IF(independentFieldUserNumber/=independentField%userNumber) THEN
        localError="The specified independent field user number of "// &
          & TRIM(NumberToVString(independentFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified independent field of "// &
          & TRIM(NumberToVString(independentField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(independentFieldRegion)
      CALL Field_RegionGet(independentField,independentFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the equations set
      IF(independentFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified independent field has been created on region number "// &
          & TRIM(NumberToVString(independentFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified independent field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(independentDecomposition)
      CALL Field_DecompositionGet(independentField,independentDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,independentDecomposition)) THEN
        CALL FlagError("The specified independent field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(independentFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified independent field user number of "// &
          & TRIM(NumberToVString(independentFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    !Initialise the equations set independent
    CALL EquationsSet_IndependentInitialise(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(independentField)) equationsSet%independent%independentFieldAutoCreated=.TRUE.
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
    CALL EquationsSet_EquationsSetFieldInitialise(equationsSet,err,error,*999)
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

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(materialsField)) THEN
      !Check the materials field has been finished
      CALL Field_AssertIsFinished(materialsField,err,error,*999)
      !Check the user numbers match
      IF(materialsFieldUserNumber/=materialsField%userNumber) THEN
        localError="The specified materials field user number of "// &
          & TRIM(NumberToVString(materialsFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified materials field of "// &
          & TRIM(NumberToVString(materialsField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(materialsFieldRegion)
      CALL Field_RegionGet(materialsField,materialsFieldRegion,err,error,*998)
      !Check the field is defined on the same region as the equations set
      IF(materialsFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified materials field has been created on region number "// &
          & TRIM(NumberToVString(materialsFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified materials field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(materialsDecomposition)
      CALL Field_DecompositionGet(materialsField,materialsDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,materialsDecomposition)) THEN
        CALL FlagError("The specified materials field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(materialsFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified materials field user number of "// &
          & TRIM(NumberToVString(materialsFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    !Initialise the equations set materials
    CALL EquationsSet_MaterialsInitialise(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(materialsField)) equationsSet%materials%materialsFieldAutoCreated=.TRUE.
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
    
    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(dependentField)) THEN
      !Check the dependent field has been finished
      CALL Field_AssertIsFinished(dependentField,err,error,*999)
      !Check the user numbers match
      IF(dependentFieldUserNumber/=dependentField%userNumber) THEN
        localError="The specified dependent field user number of "// &
          & TRIM(NumberToVString(dependentFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified dependent field of "// &
          & TRIM(NumberToVString(dependentField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(dependentFieldRegion)
      CALL Field_RegionGet(dependentField,dependentFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the equations set
      IF(dependentFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified dependent field has been created on region number "// &
          & TRIM(NumberToVString(dependentFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified dependent field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,dependentDecomposition)) THEN
        CALL FlagError("The specified dependent field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(dependentFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified dependent field user number of "// &
          & TRIM(NumberToVString(dependentFieldUserNumber,"*",err,error))// &
          & " has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      equationsSet%dependent%dependentFieldAutoCreated=.TRUE.
    ENDIF
    !Initialise the setup
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

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(derivedField)) THEN
      !Check the derived field has been finished
      CALL Field_AssertIsFinished(derivedField,err,error,*999)
      !Check the user numbers match
      IF(derivedFieldUserNumber/=derivedField%userNumber) THEN
        localError="The specified derived field user number of "// &
          & TRIM(NumberToVString(derivedFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified derived field of "// &
          & TRIM(NumberToVString(derivedField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END IF
      NULLIFY(derivedFieldRegion)
      CALL Field_RegionGet(derivedField,derivedFieldRegion,err,error,*999)
      !Check the field is defined on the same region as the equations set
      IF(derivedFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified derived field has been created on region number "// &
          & TRIM(NumberToVString(derivedFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END IF
      !Check the specified derived field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(derivedDecomposition)
      CALL Field_DecompositionGet(derivedField,derivedDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,derivedDecomposition)) THEN
        CALL FlagError("The specified derived field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(derivedFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified derived field user number of "// &
          & TRIM(NumberToVString(derivedFieldUserNumber,"*",err,error))// &
          & " has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      equationsSet%derived%derivedFieldAutoCreated=.TRUE.
    ENDIF
    CALL EquationsSet_DerivedInitialise(equationsSet,err,error,*999)
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

  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE EqutionsSet_EquationsSetFieldFinalise(equationsSetField,err,error,*)

    !Argument variables
    TYPE(EquationsSetEquationsSetFieldType) :: equationsSetField !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EqutionsSet_EquationsSetFieldFinalise",err,error,*999)

    NULLIFY(equationsSetField%equationsSet)
    equationsSetField%equationsSetFieldFinished=.FALSE.
    equationsSetField%equationsSetFieldAutoCreated=.FALSE.
    NULLIFY(equationsSetField%equationsSetFieldField)
    
    EXITS("EqutionsSet_EquationsSetFieldFinalise")
    RETURN
999 ERRORSEXITS("EqutionsSet_EquationsSetFieldFinalise",err,error)
    RETURN 1
    
  END SUBROUTINE EqutionsSet_EquationsSetFieldFinalise
  
  !
  !================================================================================================================================
  !
  !>Initialises the equations set field for a equations set.
  SUBROUTINE EquationsSet_EquationsSetFieldInitialise(equationsSet,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_EquationsSetFieldInitialise",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    equationsSet%equationsSetField%equationsSet=>equationsSet
    equationsSet%equationsSetField%equationsSetFieldFinished=.FALSE.
    equationsSet%equationsSetField%equationsSetFieldAutoCreated=.TRUE.
    NULLIFY(equationsSet%equationsSetField%equationsSetFieldField)
        
    EXITS("EquationsSet_EquationsSetFieldInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsSetFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsSetFieldInitialise

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
      IF(SIZE(equationsSet%specification,1)<2) THEN
        CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation class.", &
          & err,error,*999)
      END IF
      IF(equationsSet%specification(2) == EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
        CALL MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP(equationsSet,equationsSetSetupInfo,err,error,*999)
      ELSE
        CALL BIOELECTRIC_EQUATIONS_SET_SETUP(equationsSet,equationsSetSetupInfo,err,error,*999)
      END IF
    CASE(EQUATIONS_SET_FITTING_CLASS)
      CALL Fitting_EquationsSetSetup(equationsSet,equationsSetSetupInfo,err,error,*999)
    CASE(EQUATIONS_SET_MODAL_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
      CALL MULTI_PHYSICS_EQUATIONS_SET_SETUP(equationsSet,equationsSetSetupInfo,err,error,*999)
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,residualNumber,jacobianNumber,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx                  
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,residualNumber,jacobianNumber,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatrices_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,residualNumber,jacobianNumber,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,residualNumber,jacobianNumber,element,err,error,*999)
      CALL EquationsMatricesVector_JacobianElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
      CALL EquationsMatrices_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    INTEGER(INTG) :: residualIdx,residualVariableIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldParameterSetType), POINTER :: residualParameterSet
    TYPE(FieldVariableType), POINTER :: residualVariable
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_ResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
        
    IF(equationsSet%outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set residual evaluate: ",equationsSet%label,err,error,*999)
    ENDIF
    
    SELECT CASE(equations%linearity)
    CASE(EQUATIONS_LINEAR)
      CALL FlagError("Can not evaluate a residual for linear equations.",err,error,*999)
    CASE(EQUATIONS_NONLINEAR)
      SELECT CASE(equations%timeDependence)
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
    
    !Update the residual parameter set if it exists
    DO residualIdx=1,nonlinearMapping%numberOfResiduals
      DO residualVariableIdx=1,nonlinearMapping%numberOfResidualVariables
        NULLIFY(residualVariable)
        CALL EquationsMappingNonlinear_ResidualVariableGet(nonlinearMapping,residualIdx,residualVariableIdx,residualVariable, &
          & err,error,*999)
        NULLIFY(residualParameterSet)
        CALL FieldVariable_ParameterSetCheck(residualVariable,FIELD_RESIDUAL_SET_TYPE,residualParameterSet,err,error,*999)
        IF(ASSOCIATED(residualParameterSet)) THEN
          !Residual parameter set exists. Copy the residual vector to the residuals parameter set.
          CALL DistributedVector_Copy(nonlinearMatrices%residual,residualParameterSet%parameters,1.0_DP,err,error,*999)
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx                  
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO elementIdx=elementsMapping%internalStart,elementsMapping%internalFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime3,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime3,err,error,*999)
      userElapsed=userTime3(1)-userTime2(1)
      systemElapsed=systemTime3(1)-systemTime2(1)
      elementUserElapsed=userElapsed
      elementSystemElapsed=systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Internal elements equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Loop over the boundary and ghost elements
    DO elementIdx=elementsMapping%boundaryStart,elementsMapping%ghostFinish
      element=elementsMapping%domainList(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatricesVector_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      elementUserElapsed=elementUserElapsed+userElapsed
      elementSystemElapsed=elementSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_ELEMENT_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost elements equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average element equations assembly", &
        & elementUserElapsed/numberOfTimes,elementSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the element matrices
    CALL EquationsMatricesVector_ElementFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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

    NULLIFY(region)
    CALL EquationsSet_RegionGet(equationsSet,region,err,error,*998)
    IF(ASSOCIATED(sourceField)) THEN
      !Check the source field has been finished
      CALL Field_AssertIsFinished(sourceField,err,error,*999)
      !Check the user numbers match
      IF(sourceFieldUserNumber/=sourceField%userNumber) THEN
        localError="The specified source field user number of "// &
          & TRIM(NumberToVString(sourceFieldUserNumber,"*",err,error))// &
          & " does not match the user number of the specified source field of "// &
          & TRIM(NumberToVString(sourceField%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(sourceFieldRegion)
      CALL Field_RegionGet(sourceField,sourceFieldRegion,err,error,*998)
      !Check the field is defined on the same region as the equations set
      IF(sourceFieldRegion%userNumber/=region%userNumber) THEN
        localError="Invalid region setup. The specified source field has been created on region number "// &
          & TRIM(NumberToVString(sourceFieldRegion%userNumber,"*",err,error))// &
          & " and the specified equations set has been created on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      !Check the specified source field has the same decomposition as the geometric field
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*998)
      NULLIFY(geometricDecomposition)
      CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*998)
      NULLIFY(sourceDecomposition)
      CALL Field_DecompositionGet(sourceField,sourceDecomposition,err,error,*998)
      IF(.NOT.ASSOCIATED(geometricDecomposition,sourceDecomposition)) THEN
        CALL FlagError("The specified source field does not have the same decomposition as the geometric "// &
          & "field for the specified equations set.",err,error,*999)
      ENDIF
    ELSE
      !Check the user number has not already been used for a field in this region.
      NULLIFY(field)
      CALL Field_UserNumberFind(sourceFieldUserNumber,region,field,err,error,*999)
      IF(ASSOCIATED(field)) THEN
        localError="The specified source field user number of "// &
          & TRIM(NumberToVString(sourceFieldUserNumber,"*",err,error))// &
          & "has already been used to create a field on region number "// &
          & TRIM(NumberToVString(region%userNumber,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ENDIF
    !Initialise the equations set source
    CALL EquationsSet_SourceInitialise(equationsSet,err,error,*999)
    IF(.NOT.ASSOCIATED(sourceField)) equationsSet%source%sourceFieldAutoCreated=.TRUE.
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
    INTEGER(INTG) :: variableIdx,variableType,dirichletIdx,dirichletDOFIdx,neumannPointDOF
    INTEGER(INTG) :: conditionIdx, conditionGlobalDOF, conditionLocalDOF, myGroupComputationNodeNumber
    REAL(DP), POINTER :: fullLoads(:),currentLoads(:), prevLoads(:)
    REAL(DP) :: fullLoad, currentLoad, newLoad, prevLoad
    TYPE(BoundaryConditionsDirichletType), POINTER :: dirichletBoundaryConditions
    TYPE(BoundaryConditionsNeumannType), POINTER :: neumannBoundaryConditions
    TYPE(BoundaryConditionsPressureIncrementedType), POINTER :: pressureIncrementedBoundaryConditions
    TYPE(BoundaryConditionVariableType),   POINTER :: boundaryConditionsVariable
    TYPE(DecompositionType), POINTER :: decomposition
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldVariableType), POINTER :: dependentVariable
    TYPE(WorkGroupType), POINTER :: workGroup

    ENTERS("EquationsSet_BoundaryConditionsIncrement",err,error,*999)
   
    !Take the stored load, scale it down appropriately then apply to the unknown variables    
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions are not associated.",err,error,*999)
    
    IF(DIAGNOSTICS1) THEN
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set BC increment: ",equationsSet%label,err,error,*999)
    ENDIF

    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(decomposition)
    CALL Field_DecompositionGet(dependentField,decomposition,err,error,*999)
    NULLIFY(workGroup)
    CALL Decomposition_WorkGroupGet(decomposition,workGroup,err,error,*999)
    CALL WorkGroup_GroupNodeNumberGet(workGroup,myGroupComputationNodeNumber,err,error,*999)
    IF(.NOT.ALLOCATED(dependentField%variables)) CALL FlagError("Dependent field variables are not allocated.",err,error,*999)
    !Loop over the variables associated with this equations set
    !\todo: Looping over all field variables is not safe when volume-coupled problem is solved. Look at matrix and rhs mapping instead?
    DO variableIdx=1,dependentField%numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      NULLIFY(boundaryConditionsVariable)
      CALL BoundaryConditions_VariableExists(boundaryConditions,dependentVariable,boundaryConditionsVariable,err,error,*999)
      IF(ASSOCIATED(boundaryConditionsVariable)) THEN
        NULLIFY(domainMapping)
        CALL FieldVariable_DomainMappingGet(dependentVariable,domainMapping,err,error,*999)
        ! Check if there are any incremented conditions applied for this boundary conditions variable
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_FIXED_INCREMENTED)>0.OR. &
          & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)>0) THEN
          NULLIFY(dirichletBoundaryConditions)
          CALL BoundaryConditionsVariable_DirichletConditionsGet(boundaryConditionsVariable,dirichletBoundaryConditions, &
            & err,error,*999)
          !Get the pointer to vector holding the full and current loads
          !   full load: FIELD_BOUNDARY_CONDITIONS_SET_TYPE - holds the target load values
          !   current load: FIELD_VALUES_SET_TYPE - holds the current increment values
          NULLIFY(fullLoads)
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
            & fullLoads,err,error,*999)
          !chrm 22/06/2010: 'FIELD_BOUNDARY_CONDITIONS_SET_TYPE' does not get updated with time (update_BCs)
          !\ToDo: How can this be achieved ???
          NULLIFY(currentLoads)
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & currentLoads,err,error,*999)
          !Get full increment, calculate new load, then apply to dependent field
          DO dirichletIdx=1,boundaryConditionsVariable%NUMBER_OF_DIRICHLET_CONDITIONS
            dirichletDOFIdx=dirichletBoundaryConditions%DIRICHLET_DOF_INDICES(dirichletIdx)
            !Check whether we have an incremented boundary condition type
            SELECT CASE(boundaryConditionsVariable%conditionTypes(dirichletDOFIdx))
            CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED, &
              & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
              !Convert dof index to local index
              IF(domainMapping%globalToLocalMap(dirichletDOFIdx)%domainNumber(1)==myGroupComputationNodeNumber) THEN
                dirichletDOFIdx=domainMapping%globalToLocalMap(dirichletDOFIdx)%localNumber(1)
                IF(0<dirichletDOFIdx.AND.dirichletDOFIdx<domainMapping%ghostStart) THEN
                  fullLoad=fullLoads(dirichletDOFIdx)
                  ! Apply full load if last step, or fixed BC
                  IF(iterationNumber==maximumNumberOfIterations) THEN
                    CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
                      & dirichletDOFIdx,fullLoad,err,error,*999)
                  ELSE
                    !Calculate new load and apply to dependent field
                    currentLoad=currentLoads(dirichletDOFIdx)
                    newLoad=currentLoad+(fullLoad-currentLoad)/(maximumNumberOfIterations-iterationNumber+1)
                    CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
                      & dirichletDOFIdx,newLoad,err,error,*999)
                    IF(diagnostics1) THEN
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx",dirichletDOFIdx,err,error,*999)
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load",currentLoad,err,error,*999)
                      CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",newLoad,err,error,*999)
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
          CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          !Restore the vector handles
          CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
            & fullLoads,err,error,*999)
          CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
            & currentLoads,err,error,*999)
        ENDIF
        ! Also increment any incremented Neumann point conditions
        IF(boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
          NULLIFY(neumannBoundaryConditions)
          CALL BoundaryConditionsVariable_NeumannConditionsGet(boundaryConditionsVariable,neumannBoundaryConditions, &
            & err,error,*999)
          ! The boundary conditions parameter set contains the full values and the
          ! current incremented values are transferred to the point values vector
          DO conditionIdx=1,boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)+ &
              & boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_NEUMANN_POINT)
            conditionGlobalDOF=neumannBoundaryConditions%setDofs(conditionIdx)
            ! conditionGlobalDOF could be for non-incremented point Neumann condition
            IF(boundaryConditionsVariable%conditionTypes(conditionGlobalDOF)/= &
              & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) CYCLE
            IF(domainMapping%globalToLocalMap(conditionGlobalDOF)%domainNumber(1)== &
              & myGroupComputationNodeNumber) THEN
              conditionLocalDOF=domainMapping%globalToLocalMap(conditionGlobalDOF)%localNumber(1)
              neumannPointDOF=boundaryConditionsVariable%neumannBoundaryConditions%pointDofMapping% &
                & globalToLocalMap(conditionIdx)%localNumber(1)
              CALL Field_ParameterSetGetLocalDOF(dependentField,variableType, &
                & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,conditionLocalDOF,fullLoad,err,error,*999)
              CALL DistributedVector_ValuesSet(neumannBoundaryConditions%pointValues,neumannPointDOF, &
                & fullLoad*(REAL(iterationNumber)/REAL(maximumNumberOfIterations)), &
                & err,error,*999)
            ENDIF
          ENDDO !conditionIdx
        ENDIF

        !There might also be pressure incremented conditions
        IF (boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0) THEN
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
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_PREVIOUS_PRESSURE_SET_TYPE,prevLoads,err,error,*999)
          NULLIFY(currentLoads)
          CALL Field_ParameterSetDataGet(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE,currentLoads,err,error,*999)
          !Calculate the new load, update the old load
          IF(iterationNumber==1) THEN
            !On the first iteration, FIELD_PRESSURE_VALUES_SET_TYPE actually contains the full load
            DO conditionIdx=1,boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
              !Global dof index
              conditionGlobalDOF=pressureIncrementedBoundaryConditions%pressureIncrementedDOFIndices(conditionIdx)
              !Must convert into local dof index
              IF(domainMapping%globalToLocalMap(conditionGlobalDOF)%domainNumber(1)==myGroupComputationNodeNumber) THEN
                conditionLocalDOF=domainMapping%globalToLocalMap(conditionGlobalDOF)%localNumber(1)
                IF(0<conditionLocalDOF.AND.conditionLocalDOF<domainMapping%ghostStart) THEN
                  newLoad=currentLoads(conditionLocalDOF)
                  newLoad=newLoad/maximumNumberOfIterations
                  !Update current and previous loads
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,conditionLocalDOF,newLoad,err,error,*999)
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType, &
                    & FIELD_PREVIOUS_PRESSURE_SET_TYPE,conditionLocalDOF,0.0_dp,err,error,*999)
                  IF(diagnostics1) THEN
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx", &
                      & conditionLocalDOF,err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load", &
                      & currentLoads(conditionLocalDOF),err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",newLoad,err,error,*999)
                  ENDIF
                ENDIF !Non-ghost dof
              ENDIF !Current domain
            ENDDO !conditionIdx
          ELSE
            !Calculate the new load, keep the current load
            DO conditionIdx=1,boundaryConditionsVariable%dofCounts(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
              !This is global dof idx
              conditionGlobalDOF=pressureIncrementedBoundaryConditions%pressureIncrementedDOFIndices(conditionIdx)
              !Must convert into local dof index
              IF(domainMapping%globalToLocalMap(conditionGlobalDOF)%domainNumber(1)==myGroupComputationNodeNumber) THEN
                conditionLocalDOF=domainMapping%globalToLocalMap(conditionGlobalDOF)%localNumber(1)
                IF(0<conditionLocalDOF.AND.conditionLocalDOF<domainMapping%ghostStart) THEN
                  prevLoad=prevLoads(conditionLocalDOF)
                  currentLoad=currentLoads(conditionLocalDOF)
                  newLoad=currentLoad+(currentLoad-prevLoad)  !This may be subject to numerical errors...
                  !if (conditionIdx==1) write(*,*) "new load=",new_load
                  !Update current and previous loads
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType, &
                    & FIELD_PRESSURE_VALUES_SET_TYPE,conditionLocalDOF,newLoad,err,error,*999)
                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,variableType, &
                    & FIELD_PREVIOUS_PRESSURE_SET_TYPE,conditionLocalDOF,currentLoad,err,error,*999)
                  IF(diagnostics1) THEN
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx", &
                      & conditionLocalDOF,err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load", &
                      & currentLoads(conditionLocalDOF),err,error,*999)
                    CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",newLoad,err,error,*999)
                  ENDIF
                ENDIF !Non-ghost dof
              ENDIF !Current domain
            ENDDO !conditionIdx
          ENDIF
          !Start transfer of dofs to neighbouring domains
          CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
          !Restore the vector handles
          CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_PREVIOUS_PRESSURE_SET_TYPE,prevLoads,err,error,*999)
          CALL Field_ParameterSetDataRestore(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE,currentLoads, &
            & err,error,*999)
          !Finish transfer of dofs to neighbouring domains
          CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_PREVIOUS_PRESSURE_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_PRESSURE_VALUES_SET_TYPE,err,error,*999)
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

    IF(ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
   
    !Increment boundary conditions
    CALL EquationsSet_BoundaryConditionsIncrement(equationsSet,boundaryConditions,iterationNumber, &
      & maximumNumberOfIterations,err,error,*999)

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
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO nodeIdx=nodalMapping%internalStart,nodalMapping%internalFinish
      nodeNumber=nodalMapping%domainList(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average nodes equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations matrices and RHS vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(NodalMatrixType), POINTER :: nodalMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_NodalJacobianEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<1) &
      & CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
    
    NULLIFY(equations)    
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)

    DO residualIdx=1,nonlinearMatrices%numberOfResiduals
      NULLIFY(residualVector)
      CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)    
      DO matrixIdx=1,residualVector%numberOfJacobians
        NULLIFY(jacobianMatrix)
        CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
        SELECT CASE(jacobianMatrix%jacobianCalculationType)
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
          localError="The Jacobian calculation type of "// &
            & TRIM(NumberToVString(jacobianMatrix%jacobianCalculationType,"*",err,error))// &
            & " is not valid for matrix index number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END DO !matrixIdx
    ENDDO !residualIdx
    IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal Jacobian matrix:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node number = ",nodeNumber,err,error,*999)      
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal Jacobians:",err,error,*999)
      DO residualIdx=1,nonlinearlinearMatrices%numberOfResiduals
        NULLIFY(residualVector)
        CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,residualIdx,residualVector,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Residual number = ",residualIdx,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update residual = ",residualVector%updateNodalJacobain,err,error,*999)
        DO matrixIdx=1,residualVector%numberOfJacobians
          NULLIFY(jacobianMatrix)
          CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,matrixIdx,jacobianMatrix,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"    Jacobian number = ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"      Update Jacobian = ",jacobianMatrix%updateJacobian,err,error,*999)
          IF(jacobianMatrix%updateJacobian) THEN
            NULLIFY(nodalMatrix)
            CALL JacobianMatrix_NodalJacobianGet(jacobianMatrix,nodalMatrix,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"      Number of rows = ",nodalMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"      Number of columns = ",nodalMatrix%numberOfColumns,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"      Maximum number of rows = ",nodalMatrix%maxNumberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"      Maximum number of columns = ",nodalMatrix%maxNumberOfColumns, &
              & err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,8,8,nodalMatrix%rowDofs, &
              & '("      Row dofs     :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfColumns,8,8,nodalMatrix% &
              & columnDofs,'("      Column dofs  :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,1,1,nodalMatrix% &
              & numberOfColumns,8,8,nodalMatrix%matrix(1:nodalMatrix%numberOfRows,1:nodalMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("      Matrix','(",I2,",:)',' :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
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
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(NodalMatrixType), POINTER :: nodalMatrix
    TYPE(NodalVectorType), POINTER :: nodalVector
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_NodalResidualEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    nonlinearMatrices%nodalResidualCalculated=nodeNumber
    
    IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal residual matrices and vectors:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node number = ",nodeNumber,err,error,*999)
      linearMatrices=>vectorMatrices%linearMatrices
      IF(ASSOCIATED(linearMatrices)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Linear matrices:",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Number of node matrices = ",linearMatrices%numberOfLinearMatrices,err,error,*999)
        DO matrixIdx=1,linearMatrices%numberOfLinearMatrices
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node matrix : ",matrixIdx,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update matrix = ",linearMatrices%matrices(matrixIdx)%ptr%updateMatrix, &
            & err,error,*999)
          IF(linearMatrices%matrices(matrixIdx)%ptr%updateMatrix) THEN
            nodalMatrix=>linearMatrices%matrices(matrixIdx)%ptr%nodalMatrix
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",nodalMatrix%numberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",nodalMatrix%numberOfColumns,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",nodalMatrix%maxNumberOfRows,err,error,*999)
            CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",nodalMatrix%maxNumberOfColumns, &
              & err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,8,8,nodalMatrix%rowDofs, &
              & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfColumns,8,8,nodalMatrix% &
              & columnDofs,'("  Column dofs      :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
            CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,1,1,nodalMatrix% &
              & numberOfColumns,8,8,nodalMatrix%matrix(1:nodalMatrix%numberOfRows,1:nodalMatrix% &
              & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)','     :",8(X,E13.6))', &
              & '(20X,8(X,E13.6))',err,error,*999)
          ENDIF
        ENDDO !matrixIdx
      ENDIF
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Node residual vector:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",nonlinearMatrices%updateResidual,err,error,*999)
      IF(nonlinearMatrices%updateResidual) THEN
        nodalVector=>nonlinearMatrices%nodalResidual
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",nodalVector%numberOfRows,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",nodalVector%maxNumberOfRows,err,error,*999)
        CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%rowDofs, &
          & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
        CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%vector, &
          & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
      ENDIF
      rhsVector=>vectorMatrices%rhsVector
      IF(ASSOCIATED(rhsVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Node RHS vector :",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",rhsVector%updateVector,err,error,*999)
        IF(rhsVector%updateVector) THEN
          nodalVector=>rhsVector%nodalVector
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",nodalVector%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",nodalVector%maxNumberOfRows,err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%rowDofs, &
            & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%vector, &
            & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
        ENDIF
      ENDIF
      sourceVector=>vectorMatrices%sourceVector
      IF(ASSOCIATED(sourceVector)) THEN
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Node source vector :",err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update vector = ",sourceVector%updateVector,err,error,*999)
        IF(sourceVector%updateVector) THEN
          nodalVector=>sourceVector%nodalVector
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",nodalVector%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",nodalVector%maxNumberOfRows,err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%rowDofs, &
            & '("  Row dofs         :",8(X,I13))','(20X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalVector%numberOfRows,8,8,nodalVector%vector, &
            & '("  Vector(:)        :",8(X,E13.6))','(20X,8(X,E13.6))',err,error,*999)
        ENDIF
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
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the nodes
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime2,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime2,err,error,*999)
      userElapsed=userTime2(1)-userTime1(1)
      systemElapsed=systemTime2(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Setup and initialisation",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    numberOfTimes=0
    !Loop over the internal nodes
    DO nodeIdx=nodalMapping%internalStart,nodalMapping%internalFinish
      nodeNumber=nodalMapping%domainList(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_JacobianNodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_JacobianNodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average node equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations Jacobian if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) &
      & CALL EquationsMatrices_JacobianOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
       
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
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
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

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatricesVector_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatricesVector_NodalInitialise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    DO nodeIdx=nodalMapping%internalStart,nodalMapping%internalFinish
      nodeNumber=nodalMapping%domainList(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatricesVector_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatricesVector_NodeAdd(vectorMatrices,err,error,*999)
    ENDDO !nodeIdx
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
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
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime5,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime5,err,error,*999)
      userElapsed=userTime5(1)-userTime3(1)
      systemElapsed=systemTime5(1)-systemTime3(1)
      nodeUserElapsed=nodeUserElapsed+userElapsed
      nodeSystemElapsed=nodeSystemElapsed+systemElapsed
      IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) &
        & CALL Profiling_TimingsOutput(0,"",userElapsed,systemElapsed,err,error,*999)
      CALL Profiling_TimingsOutput(1,"Boundary+ghost nodes equations assembly",userElapsed,systemElapsed,err,error,*999)
      IF(numberOfTimes>0) CALL Profiling_TimingsOutput(1,"Average node equations assembly", &
        & nodeUserElapsed/numberOfTimes,nodeSystemElapsed/numberOfTimes,err,error,*999)
    ENDIF
    !Finalise the nodal matrices
    CALL EquationsMatricesVector_NodalFinalise(vectorMatrices,err,error,*999)
    !Output timing information if required
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime6,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime6,err,error,*999)
      userElapsed=userTime6(1)-userTime1(1)
      systemElapsed=systemTime6(1)-systemTime1(1)
      CALL Profiling_TimingsOutput(1,"Total equations assembly",userElapsed,systemElapsed,err,error,*999)
    ENDIF
    !Output equations residual vector if required
    IF(equations%outputType>=EQUATIONS_MATRIX_OUTPUT) THEN
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
