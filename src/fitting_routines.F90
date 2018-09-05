!> \file
!> \author Sebastian Krittian
!> \brief This module handles all fitting routines.
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

!>This module handles all fitting routines.
MODULE FittingRoutines

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE DARCY_EQUATIONS_ROUTINES, ONLY: idebug1
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FLUID_MECHANICS_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE Maths
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
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

  PUBLIC Fitting_EquationsSetSetup
  PUBLIC Fitting_EquationsSetSpecificationSet
  PUBLIC Fitting_EquationsSetSolutionMethodSet

  PUBLIC Fitting_ProblemSetup
  PUBLIC Fitting_ProblemSpecificationSet

  PUBLIC Fitting_FiniteElementCalculate

  PUBLIC Fitting_PreSolve
  PUBLIC Fitting_PostSolve

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Galerkin projection finite element equations set.
   SUBROUTINE Fitting_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: ng,mh,mhs,ms,nh,nhs,ns,mi,ni
    INTEGER(INTG) :: dependentComponentColumnIdx,dependentComponentRowIdx,dependentElementParameterColumnIdx, &
      & dependentElementParameterRowIdx,dependentParameterColumnIdx,dependentParameterRowIdx,gaussPointIdx, &
      & meshComponentRow,meshComponentColumn,numberOfDataComponents
    REAL(DP) :: rwg,sum,jacobianGaussWeight
    REAL(DP) :: basisFunctionRow,basisFunctionColumn,PGM,PGN,PGMSI(3),PGNSI(3)
    REAL(DP) :: uValue(3)
    REAL(DP) :: phiM,phiN
    TYPE(DataProjectionType), POINTER :: dataProjection
    TYPE(DECOMPOSITION_TOPOLOGY_TYPE), POINTER :: decompositionTopology
    TYPE(DecompositionDataPointsType), POINTER :: dataPoints
    TYPE(BASIS_TYPE), POINTER :: dependentBasis,geometricBasis,sourceBasis,dependentBasisRow,dependentBasisColumn
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,independentField,materialsField,sourceField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dataVariable,dataWeightVariable,dependentVariable,fieldVariable
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: materialsInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: geometricInterpolatedPoint
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: referenceGeometricInterpolatedPoint
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: quadratureScheme,quadratureSchemeColumn,quadratureSchemeRow
    TYPE(VARYING_STRING) :: localError

    REAL(DP), POINTER :: independentVectorParameters(:),independentWeightParameters(:)
    REAL(DP) :: projectionXi(3)
    REAL(DP) :: porosity0, porosity, permOverVisParam0, permOverVisParam,tauParam,kappaParam
    REAL(DP) :: tension,curvature
    REAL(DP) :: materialFact
    REAL(DP) :: dXdY(3,3), dXdXi(3,3), dYdXi(3,3), dXidY(3,3), dXidX(3,3)
    REAL(DP) :: Jxy, Jyxi
    REAL(DP) :: dataPointWeight(99),dataPointVector(99)
    INTEGER(INTG) :: derivative_idx, component_idx, xi_idx, numberOfDimensions,smoothingType
    INTEGER(INTG) :: dataPointIdx,dataPointUserNumber,dataPointLocalNumber,dataPointGlobalNumber
    INTEGER(INTG) :: numberOfXi
    INTEGER(INTG) :: componentIdx
    INTEGER(INTG) :: dependentVariableType,variableType,localDof

    INTEGER(INTG) numberDofs
    INTEGER(INTG) meshComponent1,meshComponent2

    ENTERS("Fitting_FiniteElementCalculate",err,error,*999)

    NULLIFY(dependentBasis,geometricBasis)
    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(linearMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(linearMatrices)
    NULLIFY(rhsVector)
    NULLIFY(equationsMatrix)
    NULLIFY(dependentField,geometricField,materialsField)
    NULLIFY(dataPoints)
    NULLIFY(dataProjection)
    NULLIFY(decompositionTopology)
    NULLIFY(independentField)
    NULLIFY(independentVectorParameters)
    NULLIFY(independentWeightParameters)
    NULLIFY(fieldVariable)
    NULLIFY(dependentVariable)
    NULLIFY(quadratureScheme)
    NULLIFY(geometricInterpolatedPoint,materialsInterpolatedPoint)

    dataPointVector = 0.0_DP
    dataPointWeight = 0.0_DP

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%equations
      IF(ASSOCIATED(equations)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        IF(.NOT.ALLOCATED(equationsSet%specification)) &
          & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        IF(SIZE(equationsSet%specification,1)/=4) &
          & CALL FlagError("Equations set specification must have four entries for a fitting type equations set.",err,error,*999)
        smoothingType=equationsSet%specification(4)
        SELECT CASE(equationsSet%specification(2))
        CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_DATA_POINT_FITTING_SUBTYPE)
            geometricField=>equations%interpolation%geometricField
            dependentField=>equations%interpolation%dependentField
            independentField=>equations%interpolation%independentField
            vectorMatrices=>vectorEquations%vectorMatrices
            linearMatrices=>vectorMatrices%linearMatrices
            equationsMatrix=>linearMatrices%matrices(1)%ptr
            rhsVector=>vectorMatrices%rhsVector
            vectorMapping=>vectorEquations%vectorMapping
            linearMapping=>vectorMapping%linearMapping
            dependentVariable=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
            dependentVariableType=dependentVariable%VARIABLE_TYPE
            dataVariable=>independentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
            dataWeightVariable=>independentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%ptr
            dependentBasis=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            geometricBasis=>geometricField%decomposition%domain(geometricField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            quadratureScheme=>dependentBasis%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & dependentInterpParameters(dependentVariableType)%ptr,err,error,*999)
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL Field_NumberOfComponentsGet(independentField,FIELD_U_VARIABLE_TYPE,numberOfDataComponents,err,error,*999)
            IF(numberOfDataComponents>99) CALL FlagError("Increase the size of the data point vectors.",err,error,*999)
            numberOfXi = dependentBasis%NUMBER_OF_XI

            !Get data point vector parameters
            CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentVectorParameters,err,error,*999)
            !Get data point weight parameters
            CALL Field_ParameterSetDataGet(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentWeightParameters,err,error,*999)

            !===============================
            ! D a t a   P o i n t   F i t
            !===============================
            dataProjection=>independentField%dataProjection
            IF(.NOT.ASSOCIATED(dataProjection)) &
              & CALL FlagError("Data projection is not associated on independent field.",err,error,*999)
            decompositionTopology=>independentField%decomposition%topology
            IF(ASSOCIATED(decompositionTopology)) THEN
              dataPoints=>decompositionTopology%dataPoints
              IF(.NOT.ASSOCIATED(dataPoints)) &
                & CALL FlagError("Data points are not associated on the decomposition topology of the independent field.", &
                & err,error,*999)
            ELSE
              CALL FlagError("Decomposition topology is not associated on the independent field.",err,error,*999)
            ENDIF
            !Loop over data points
            DO dataPointIdx=1,dataPoints%elementDataPoint(elementNumber)%numberOfProjectedData
              dataPointUserNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%userNumber
              dataPointLocalNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%localNumber
              dataPointGlobalNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%globalNumber
              ! Need to use global number to get the correct projection results
              projectionXi(1:numberOfXi) = dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi(1:numberOfXi)
              CALL Field_InterpolateXi(FIRST_PART_DERIV,projectionXi,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateXi(FIRST_PART_DERIV,projectionXi,equations%interpolation% &
                & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              !Get data point vector value and weight
              DO componentIdx=1,numberOfDataComponents
                localDof=dataVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                  & DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
                dataPointVector(componentIdx)=independentVectorParameters(localDof)
                localDof=dataWeightVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                  & DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
                dataPointWeight(componentIdx)=independentWeightParameters(localDof)
              ENDDO !componentIdx

              dependentParameterRowIdx=0
              !Loop over element rows
              DO dependentComponentRowIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                meshComponentRow=dependentVariable%components(dependentComponentRowIdx)%MESH_COMPONENT_NUMBER
                dependentBasisRow=>dependentField%decomposition%domain(meshComponentRow)%ptr%topology%elements% &
                  & elements(elementNumber)%basis
                DO dependentElementParameterRowIdx=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                  dependentParameterRowIdx=dependentParameterRowIdx+1
                  dependentParameterColumnIdx=0
                  basisFunctionRow=Basis_EvaluateXi(dependentBasisRow,dependentElementParameterRowIdx,NO_PART_DERIV, &
                    & projectionXi,err,error)
                  IF(equationsMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO dependentComponentColumnIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                      meshComponentColumn=dependentVariable%components(dependentComponentColumnIdx)%MESH_COMPONENT_NUMBER
                      dependentBasisColumn=>dependentField%decomposition%domain(meshComponentColumn)%ptr% &
                        & topology%elements%elements(elementNumber)%basis
                      DO dependentElementParameterColumnIdx=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                        dependentParameterColumnIdx=dependentParameterColumnIdx+1
                        !Treat each component as separate and independent so only calculate the diagonal blocks
                        IF(dependentComponentColumnIdx==dependentComponentRowIdx) THEN
                          basisFunctionColumn=Basis_EvaluateXi(dependentBasisColumn,dependentElementParameterColumnIdx, &
                            & NO_PART_DERIV,projectionXi,err,error)
                          sum = basisFunctionRow*basisFunctionColumn*dataPointWeight(dependentComponentRowIdx)
                          equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)= &
                            & equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)+sum
                        ENDIF
                      ENDDO !dependentElementParameterColumnIdx
                    ENDDO !dependentComponentColumnIdx
                  ENDIF
                  IF(rhsVector%updateVector) THEN
                    sum = basisFunctionRow*dataPointVector(dependentComponentRowIdx)*dataPointWeight(dependentComponentRowIdx)
                    rhsVector%elementVector%vector(dependentParameterRowIdx)= &
                      & rhsVector%elementVector%vector(dependentParameterRowIdx)+sum
                  ENDIF
                ENDDO !dependentElementParameterRowIdx
              ENDDO !dependentComponentRowIdx
            ENDDO !dataPointIdx

            !Restore data point vector parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentVectorParameters,err,error,*999)
            !Restore data point weight parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentWeightParameters,err,error,*999)

            SELECT CASE(smoothingType)
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
              IF(equationsMatrix%updateMatrix) THEN
                materialsField=>equations%interpolation%materialsField
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
                  & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                !Loop over Gauss points
                DO gaussPointIdx=1,quadratureScheme%NUMBER_OF_GAUSS
                  !Interpolate fields
                  CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                    & equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                  CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                    & equations%interpolation%dependentInterpPoint(dependentVariableType)%ptr,err,error,*999)
                  CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                    & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                  !Get Sobolev smoothing parameters from interpolated material field
                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                    & equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                  tauParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(1,NO_PART_DERIV)
                  kappaParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(2,NO_PART_DERIV)
                  jacobianGaussWeight=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                    & quadratureScheme%GAUSS_WEIGHTS(gaussPointIdx)

                  !Loop over field components
                  dependentParameterRowIdx=0
                  DO dependentComponentRowIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                    !Loop over element rows
                    meshComponentRow=dependentVariable%components(dependentComponentRowIdx)%MESH_COMPONENT_NUMBER
                    dependentBasisRow=>dependentField%decomposition%domain(meshComponentRow)%ptr% &
                      & topology%elements%elements(elementNumber)%basis
                    quadratureSchemeRow=>dependentBasisRow%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                    DO dependentElementParameterRowIdx=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                      dependentParameterRowIdx=dependentParameterRowIdx+1
                      dependentParameterColumnIdx=0
                      !Loop over element columns
                      DO dependentComponentColumnIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                        meshComponentColumn=dependentVariable%components(dependentComponentColumnIdx)%MESH_COMPONENT_NUMBER
                        dependentBasisColumn=>dependentField%decomposition%domain(meshComponentColumn)%ptr% &
                          & topology%elements%elements(elementNumber)%basis
                        quadratureSchemeColumn=>dependentBasisColumn%quadrature%QUADRATURE_SCHEME_MAP( &
                          & BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        DO dependentElementParameterColumnIdx=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                          dependentParameterColumnIdx=dependentParameterColumnIdx+1

                          !Calculate Sobolev surface tension and curvature smoothing terms
                          tension = tauParam*2.0_DP* ( &
                            & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1, &
                            & gaussPointIdx)* &
                            & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1, &
                            & gaussPointIdx))
                          curvature = kappaParam*2.0_DP* ( &
                            & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S1, &
                            & gaussPointIdx)* &
                            & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S1, &
                            & gaussPointIdx))
                          IF(numberOfXi > 1) THEN
                            tension = tension + tauParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2, &
                              & gaussPointIdx)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2, &
                              & gaussPointIdx))
                            curvature = curvature + kappaParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2_S2, &
                              & gaussPointIdx)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2_S2, &
                              & gaussPointIdx) + &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S2, &
                              & gaussPointIdx)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S2, &
                              & gaussPointIdx))
                            IF(numberOfXi > 2) THEN
                              tension = tension + tauParam*2.0_DP* ( &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S3, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S3, &
                                & gaussPointIdx))
                              curvature = curvature + kappaParam*2.0_DP* ( &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S3_S3, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S3_S3, &
                                & gaussPointIdx)+ &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S3, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S3, &
                                & gaussPointIdx)+ &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2_S3, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2_S3, &
                                & gaussPointIdx))
                            ENDIF ! 3D
                          ENDIF ! 2 or 3D
                          sum = (tension + curvature) * jacobianGaussWeight

                          equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)= &
                            equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)+sum

                        ENDDO !dependentElementParameterColumnIdx
                      ENDDO !dependentComponentColumnIdx
                    ENDDO !dependentElementParameterRowIdx
                  ENDDO !dependentComponentRowIdx
                ENDDO !gaussPointIdx
              ENDIF

            CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
            &  EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
            dependentField=>equations%interpolation%dependentField
            independentField=>equations%interpolation%independentField
            dataProjection=>independentField%dataProjection
            IF(.NOT.ASSOCIATED(dataProjection)) CALL FlagError("Data projection is not associated on independent field.", &
              & err,error,*999)
            decompositionTopology=>independentField%decomposition%topology
            IF(.NOT.ASSOCIATED(decompositionTopology))  &
              & CALL FlagError("Decomposition topology is not associated on the independent field.",err,error,*999)
            dataPoints=>decompositionTopology%dataPoints
            IF(.NOT.ASSOCIATED(dataPoints)) &
              & CALL FlagError("Data points are not associated on the decomposition topology of the independent field.", &
              & err,error,*999)
            geometricField=>equations%interpolation%geometricField
            materialsField=>equations%interpolation%materialsField
            vectorMatrices=>vectorEquations%vectorMatrices
            linearMatrices=>vectorMatrices%linearMatrices
            equationsMatrix=>linearMatrices%matrices(1)%ptr
            rhsVector=>vectorMatrices%rhsVector
            vectorMapping=>vectorEquations%vectorMapping
            linearMapping=>vectorMapping%linearMapping
            dependentVariable=>linearMapping%equationsMatrixToVarMaps(1)%variable
            dependentVariableType=dependentVariable%VARIABLE_TYPE
            dependentBasis=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            geometricBasis=>geometricField%decomposition%domain(geometricField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            quadratureScheme=>dependentBasis%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & dependentInterpParameters(dependentVariableType)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            numberOfXi = dependentBasis%NUMBER_OF_XI
            projectionXi=0.0_DP
            ! Get data point vector parameters
            CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentVectorParameters,err,error,*999)
            ! Get data point weight parameters
            CALL Field_ParameterSetDataGet(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentWeightParameters,err,error,*999)

            !===========================================
            ! D a t a   P o i n t   V e c t o r    F i t
            !===========================================
            !Loop over data points
            DO dataPointIdx=1,dataPoints%elementDataPoint(elementNumber)%numberOfProjectedData
              dataPointUserNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%userNumber
              dataPointLocalNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%localNumber
              dataPointGlobalNumber = dataPoints%elementDataPoint(elementNumber)%dataIndices(dataPointIdx)%globalNumber
              ! Need to use global number to get the correct projection results
              projectionXi(1:numberOfXi) = dataProjection%dataProjectionResults(dataPointGlobalNumber)%elementXi(1:numberOfXi)
              CALL Field_InterpolateXi(FIRST_PART_DERIV,projectionXi,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateXi(FIRST_PART_DERIV,projectionXi,equations%interpolation% &
                & dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

              ! Get data point vector value
              variableType=independentField%variables(1)%VARIABLE_TYPE
              fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%ptr
!!TODO: Shouldn't this be over the number of dimensions???
              DO componentIdx=1,numberOfXi
                localDof=fieldVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                  & DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
                dataPointVector(componentIdx)=independentVectorParameters(localDof)
              ENDDO

              variableType=independentField%variables(2)%VARIABLE_TYPE
              fieldVariable=>independentField%VARIABLE_TYPE_MAP(variableType)%ptr
              localDof=fieldVariable%components(1)%PARAM_TO_DOF_MAP%DATA_POINT_PARAM2DOF_MAP%DATA_POINTS(dataPointLocalNumber)
              dataPointWeight(1)=independentWeightParameters(localDof)

              mhs=0
              !Loop over element rows
              DO mh=1,dependentVariable%NUMBER_OF_COMPONENTS
                meshComponent1=dependentVariable%components(mh)%MESH_COMPONENT_NUMBER
                dependentBasisRow=>dependentField%decomposition%domain(meshComponent1)%ptr%topology%elements% &
                  & elements(elementNumber)%basis
                DO ms=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  PGM=Basis_EvaluateXi(dependentBasisRow,ms,NO_PART_DERIV,projectionXi,err,error)
                  IF(equationsMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO nh=1,dependentVariable%NUMBER_OF_COMPONENTS
                      meshComponent2=dependentVariable%components(nh)%MESH_COMPONENT_NUMBER
                      dependentBasisColumn=>dependentField%decomposition%domain(meshComponent2)%ptr% &
                        & topology%elements%elements(elementNumber)%basis
                      DO ns=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        PGN=Basis_EvaluateXi(dependentBasisColumn,ns,NO_PART_DERIV,projectionXi,err,error)
                        sum=0.0_DP
                        IF(mh==nh) THEN
                          sum = sum + PGM * PGN * dataPointWeight(1)
                        ENDIF
                        equationsMatrix%elementMatrix%matrix(mhs,nhs)=equationsMatrix%elementMatrix%matrix(mhs,nhs)+sum
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                  sum=0.0_DP
                  IF(rhsVector%updateVector) THEN
                    sum = sum + PGM*dataPointVector(mh)*dataPointWeight(1)
                    rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs) + sum
                  ENDIF
                ENDDO !ms
              ENDDO !mh
            ENDDO !dataPointIdx

            !Restore data point vector parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentVectorParameters,err,error,*999)
            !Restore data point weight parameters
            CALL Field_ParameterSetDataRestore(independentField,FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & independentWeightParameters,err,error,*999)

            SELECT CASE(smoothingType)
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
            !===========================================
            ! S o b o l e v   S m o o t h i n g
            !===========================================
            !Loop over gauss points
            DO ng=1,quadratureScheme%NUMBER_OF_GAUSS
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & dependentInterpPoint(dependentVariableType)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              tauParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(1,NO_PART_DERIV)
              kappaParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(2,NO_PART_DERIV)
              !Loop over field components
              jacobianGaussWeight=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                & quadratureScheme%GAUSS_WEIGHTS(ng)

              mhs=0
              DO mh=1,dependentVariable%NUMBER_OF_COMPONENTS
                !Loop over element rows
                meshComponent1=dependentVariable%components(mh)%MESH_COMPONENT_NUMBER
                dependentBasisRow=>dependentField%decomposition%domain(meshComponent1)%ptr% &
                  & topology%elements%elements(elementNumber)%basis
                quadratureSchemeRow=>dependentBasisRow%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                DO ms=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF(equationsMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO nh=1,dependentVariable%NUMBER_OF_COMPONENTS
                      meshComponent2=dependentVariable%components(nh)%MESH_COMPONENT_NUMBER
                      dependentBasisColumn=>dependentField%decomposition%domain(meshComponent2)%ptr% &
                        & topology%elements%elements(elementNumber)%basis
                      quadratureSchemeColumn=>dependentBasisColumn%quadrature%QUADRATURE_SCHEME_MAP( &
                        & BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO ns=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        sum = 0.0_DP

                        !Calculate sobolev surface tension and curvature smoothing terms
                        tension = tauParam*2.0_DP* ( &
                          & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                          & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng))
                        curvature = kappaParam*2.0_DP* ( &
                          & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                          & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng))

                        IF(dependentVariable%NUMBER_OF_COMPONENTS > 1) THEN
                          tension = tension + tauParam*2.0_DP* ( &
                            & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                            & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng))
                          curvature = curvature + kappaParam*2.0_DP* ( &
                            & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                            & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng) + &
                            & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                            & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng))

                          IF(dependentVariable%NUMBER_OF_COMPONENTS > 2) THEN
                            tension = tension + tauParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng))
                            curvature = curvature + kappaParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))
                          ENDIF ! 3D
                        ENDIF ! 2 or 3D

                        ! Add in smoothing terms to the element matrix
                        equationsMatrix%elementMatrix%matrix(mhs,nhs) = &
                          & equationsMatrix%elementMatrix%matrix(mhs,nhs) + (tension + curvature) * jacobianGaussWeight

                      ENDDO !ns
                    ENDDO !nh
                  ENDIF ! update matrix
                ENDDO !ms
              ENDDO !mh
            ENDDO !ng
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
            & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Galerkin projection is put in.
            !Store all these in equations matrices/somewhere else?????
            dependentField=>equations%interpolation%dependentField
            geometricField=>equations%interpolation%geometricField
            materialsField=>equations%interpolation%materialsField

            vectorMatrices=>vectorEquations%vectorMatrices
            linearMatrices=>vectorMatrices%linearMatrices
            equationsMatrix=>linearMatrices%matrices(1)%ptr
            rhsVector=>vectorMatrices%rhsVector
            vectorMapping=>vectorEquations%vectorMapping
            linearMapping=>vectorMapping%linearMapping
            fieldVariable=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
            dependentVariableType=fieldVariable%VARIABLE_TYPE

            dependentBasis=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            geometricBasis=>geometricField%decomposition%domain(geometricField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis

            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
              & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
              & equations%interpolation%materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

            quadratureScheme=>dependentBasis%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

            !--- Loop over gauss points
            DO ng=1,quadratureScheme%NUMBER_OF_GAUSS

              !--- Interpolation of Reference Geometry
              CALL Field_InterpolationParametersElementGet(FIELD_INITIAL_VALUES_SET_TYPE,elementNumber, &
                & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              referenceGeometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & referenceGeometricInterpolatedPoint,err,error,*999)
              !--- Retrieve local map dYdXi
              DO component_idx=1,dependentBasis%NUMBER_OF_XI
                DO xi_idx=1,dependentBasis%NUMBER_OF_XI
                  derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                  dYdXi(component_idx,xi_idx)=referenceGeometricInterpolatedPoint%values(component_idx,derivative_idx) !dy/dxi (y = referential)
                ENDDO
              ENDDO

              !--- Interpolation of (actual) Geometry and Metrics
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
                & equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              geometricInterpolatedPoint=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & geometricInterpolatedPoint,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI, &
                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              !--- Retrieve local map dXdXi
              DO component_idx=1,dependentBasis%NUMBER_OF_XI
                DO xi_idx=1,dependentBasis%NUMBER_OF_XI
                  derivative_idx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xi_idx) !2,4,7
                  dXdXi(component_idx,xi_idx)=geometricInterpolatedPoint%values(component_idx,derivative_idx) !dx/dxi
                ENDDO
              ENDDO

              !--- Compute deformation gradient tensor dXdY and its Jacobian Jxy
              CALL Invert(dYdXi,dXidY,Jyxi,err,error,*999) !dy/dxi -> dxi/dy
              CALL MatrixProduct(dXdXi,dXidY,dXdY,err,error,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
              CALL Determinant(dXdY,Jxy,err,error,*999)

              !--- Interpolation of Materials Field
              materialsInterpolatedPoint => equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                & materialsInterpolatedPoint,err,error,*999)

              !--- Retrieve reference material parameters:
              porosity0            = materialsInterpolatedPoint%values(1,NO_PART_DERIV)
              permOverVisParam0 = materialsInterpolatedPoint%values(2,NO_PART_DERIV)

              !--- Material dependence on structural deformation
              IF( ABS(Jxy) > ZERO_TOLERANCE ) THEN
                porosity = 1.0_DP - ( 1.0_DP - porosity0 ) / Jxy
              ELSE
                CALL FlagError("Jacobian Jxy is zero.",err,error,*999)
              END IF

              IF(equationsSet%specification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE) THEN
                permOverVisParam = permOverVisParam0
              ELSE
                materialFact = ( Jxy * porosity / porosity0 )**2.0_DP
                permOverVisParam = materialFact * permOverVisParam0
                !material modeling could use gradient information, or solve some PDE
              END IF

              IF(DIAGNOSTICS2) THEN
                IF(idebug1) THEN
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"geometricInterpPointMetrics%jacobian = ", &
                    & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"Jxy = ",Jxy,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"porosity = ",porosity,err,error,*999)
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"permOverVisParam = ",permOverVisParam,err,error,*999)
                  CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
                  idebug1 = .FALSE.
                ENDIF
              ENDIF

!!TODO: Think about symmetric problems.
              rwg=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                & quadratureScheme%GAUSS_WEIGHTS(ng)

              !Loop over field components
              mhs=0
              DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
                !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
                DO ms=1,dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF(equationsMatrix%updateMatrix) THEN

                    !Loop over element columns
                    DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                      DO ns=1,dependentBasis%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1

                        PGM=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        PGN=quadratureScheme%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                        sum = 0.0_DP
                        IF(mh==nh) THEN
                          sum = sum + PGM * PGN
                        ENDIF
                        equationsMatrix%elementMatrix%matrix(mhs,nhs) = &
                          & equationsMatrix%elementMatrix%matrix(mhs,nhs) + sum * rwg

                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                  IF(rhsVector%updateVector) THEN
                    PGM=quadratureScheme%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)

                    sum = 0.0_DP
                    IF(mh==1) THEN
                      sum = sum + PGM * porosity
                    ELSE IF(mh==2) THEN
                      sum = sum + PGM * permOverVisParam
                    END IF
                    rhsVector%elementVector%vector(mhs) = rhsVector%elementVector%vector(mhs) + sum * rwg
                  ENDIF

                ENDDO !ms
              ENDDO !mh
            ENDDO !ng


            !-----------------------------------------------------------------------------------------------------------------------------------
            ! CHECK STIFFNESS MATRIX AND RHS VECTOR WITH CMHEART
            IF(DIAGNOSTICS5) THEN
              IF( elementNumber == 1 ) THEN
                numberDofs = 0
                DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
                  meshComponent1 = fieldVariable%components(mh)%MESH_COMPONENT_NUMBER
                  dependentBasisRow => dependentField%decomposition%domain(meshComponent1)%ptr% &
                    & topology%elements%elements(elementNumber)%basis
                  numberDofs = numberDofs + dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                END DO

                CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE,"Element Matrix for element number 1 (Galerkin Projection):",err,error,*999)
                DO mhs=1,numberDofs
                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"row number = ",mhs,err,error,*999)
                  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberDofs,numberDofs,numberDofs,&
                    & equationsMatrix%elementMatrix%matrix(mhs,:), &
                    & '("",4(X,E13.6))','4(4(X,E13.6))',err,error,*999)
                  CALL WriteString(DIAGNOSTIC_OUTPUT_TYPE," ",err,error,*999)
                END DO
              END IF
            END IF

          CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
            & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
            dependentField=>equations%interpolation%dependentField
            geometricField=>equations%interpolation%geometricField
            materialsField=>equations%interpolation%materialsField
            sourceField=>equations%interpolation%sourceField
            vectorMatrices=>vectorEquations%vectorMatrices
            linearMatrices=>vectorMatrices%linearMatrices
            equationsMatrix=>linearMatrices%matrices(1)%ptr
            rhsVector=>vectorMatrices%rhsVector
            sourceVector=>vectorMatrices%sourceVector
            vectorMapping=>vectorEquations%vectorMapping
            linearMapping=>vectorMapping%linearMapping
            fieldVariable=>linearMapping%equationsMatrixToVarMaps(1)%VARIABLE
            dependentVariableType=fieldVariable%VARIABLE_TYPE
            dependentBasis=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            geometricBasis=>geometricField%decomposition%domain(geometricField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            sourceBasis=>sourceField%decomposition%domain(sourceField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            quadratureScheme=>dependentBasis%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & dependentInterpParameters(dependentVariableType)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            !Loop over gauss points
            DO ng=1,quadratureScheme%NUMBER_OF_GAUSS
              !             CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & dependentInterpPoint(dependentVariableType)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              tauParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(1,NO_PART_DERIV)
              kappaParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(2,NO_PART_DERIV)
              ! WRITE(*,*)'tauParam ',tauParam
              uValue=0.0_DP
              IF(sourceVector%updateVector) THEN
                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
                  & equations%interpolation%sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                  & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                uValue(1)=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(1,NO_PART_DERIV)
                uValue(2)=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(2,NO_PART_DERIV)
                IF(dependentBasis%NUMBER_OF_XI==3) THEN
                  uValue(3)=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(3,NO_PART_DERIV)
                ENDIF
              ENDIF
              !Calculate rwg.
              rwg=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                & quadratureScheme%GAUSS_WEIGHTS(ng)
              !Loop over field components
              mhs=0
              DO mh=1,fieldVariable%NUMBER_OF_COMPONENTS
                !Loop over element rows
                meshComponent1=fieldVariable%components(mh)%MESH_COMPONENT_NUMBER
                dependentBasisRow=>dependentField%decomposition%domain(meshComponent1)%ptr% &
                  & topology%elements%elements(elementNumber)%basis
                quadratureSchemeRow=>dependentBasisRow%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                DO ms=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  IF(equationsMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO nh=1,fieldVariable%NUMBER_OF_COMPONENTS
                      meshComponent2=fieldVariable%components(nh)%MESH_COMPONENT_NUMBER
                      dependentBasisColumn=>dependentField%decomposition%domain(meshComponent2)%ptr% &
                        & topology%elements%elements(elementNumber)%basis
                      quadratureSchemeColumn=>dependentBasisColumn%quadrature%QUADRATURE_SCHEME_MAP &
                        & (BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO ns=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        PGM=quadratureSchemeRow%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        PGN=quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                        DO ni=1,dependentBasisColumn%NUMBER_OF_XI
                          DO mi=1,dependentBasisRow%NUMBER_OF_XI
                            dXidX(mi,ni)=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & dXi_dX(mi,ni)
                          END DO
                          PGMSI(ni)=quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          PGNSI(ni)=quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        END DO !ni
                        sum = 0.0_DP
                        !Calculate sum
                        IF(equationsSet%specification(3)==equations_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                          & equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                          IF(mh==nh) THEN
                            !This stiffness matrix contribution is without "integration" means ng=nd in fact = least square!
                            sum = sum + PGM * PGN
                          ENDIF
                          !
                          !                         IF(mh==nh) THEN
                          !                           !This stiffness matrix happens with "integration" so the integral error is reduced!
                          !                           sum = sum + PGM * PGN * rwg
                          !                         ENDIF
                          !REDUCED SOBOLEV SMOOTHING
                          !This stiffness matrix contribution is with "integration" means ng=ng in fact!
                          SELECT CASE(smoothingType)
                          CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                            !Do nothing
                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
                            sum = sum +    ( &
                              & tauParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng)) +&
                              & kappaParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))) !&
                            ! no weighting either?
                            !                             & * rwg

                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE DEFAULT
                            localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                              & " is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT

                          equationsMatrix%elementMatrix%matrix(mhs,nhs) = &
                            & equationsMatrix%elementMatrix%matrix(mhs,nhs) + sum

                        ELSE IF(equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                          & equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                          IF(mh==nh.AND.mh<=numberOfDimensions) sum = sum + PGM * PGN

                          SELECT CASE(smoothingType)
                          CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                            !Do nothing
                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
                            !REDUCED SOBOLEV SMOOTHING
                            !This stiffness matrix contribution is with "integration" means ng=ng in fact!
                            sum = sum +    ( &
                              & tauParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3,ng)) +&
                              & kappaParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S1,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S1,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S3_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S3_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S2,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S2,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S1_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S1_S3,ng)+ &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(ms,PART_DERIV_S2_S3,ng)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(ns,PART_DERIV_S2_S3,ng))) !&
                            ! no weighting either?
                            !                             & * rwg


                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE DEFAULT
                            localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                              & " is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT

                          equationsMatrix%elementMatrix%matrix(mhs,nhs) = &
                            & equationsMatrix%elementMatrix%matrix(mhs,nhs) + sum

                          IF(nh==fieldVariable%NUMBER_OF_COMPONENTS.AND.mh<=numberOfDimensions) THEN
                            sum=0.0_DP
                            !Calculate sum
                            DO ni=1,dependentBasisRow%NUMBER_OF_XI
                              sum=sum+PGN*PGMSI(ni)*dXidX(ni,mh)
                            ENDDO !ni
                            equationsMatrix%elementMatrix%matrix(mhs,nhs) = &
                              & equationsMatrix%elementMatrix%matrix(mhs,nhs) + sum * rwg
                            equationsMatrix%elementMatrix%matrix(nhs,mhs) = &
                              & equationsMatrix%elementMatrix%matrix(nhs,mhs) + sum * rwg
                          ENDIF
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh

                  ENDIF
                  IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
                  IF(sourceVector%updateVector) THEN
                    IF(equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                      & equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                      sum=0.0_DP
                      PGM=quadratureSchemeRow%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      sum=uValue(mh)*PGM
                      !                     sum=42.0_DP*PGM
                      sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+sum
                    ELSEIF(equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                      & equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                      IF(mh<=numberOfDimensions) THEN
                        sum=0.0_DP
                        PGM=quadratureSchemeRow%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                        sum=uValue(mh)*PGM
                        !                       sum=42.0_DP*PGM
                        sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+sum
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO !ms
              ENDDO !mh
            ENDDO !ng
          CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999)

          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is not valid for a data fitting equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
          SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_GAUSS_POINT_FITTING_SUBTYPE)
            dependentField=>equations%interpolation%dependentField
            independentField=>equations%interpolation%independentField
            decompositionTopology=>independentField%decomposition%topology
            geometricField=>equations%interpolation%geometricField
            vectorMatrices=>vectorEquations%vectorMatrices
            linearMatrices=>vectorMatrices%linearMatrices
            equationsMatrix=>linearMatrices%matrices(1)%ptr
            rhsVector=>vectorMatrices%rhsVector
            vectorMapping=>vectorEquations%vectorMapping
            linearMapping=>vectorMapping%linearMapping
            dependentVariable=>linearMapping%equationsMatrixToVarMaps(1)%variable
            dependentVariableType=dependentVariable%VARIABLE_TYPE
            dependentBasis=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            geometricBasis=>geometricField%decomposition%domain(geometricField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            quadratureScheme=>dependentBasis%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & dependentInterpParameters(dependentVariableType)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & independentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
              & independentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL Field_NumberOfComponentsGet(independentField,FIELD_U_VARIABLE_TYPE,numberOfDataComponents,err,error,*999)
            IF(numberOfDataComponents>99) CALL FlagError("Increase the size of the data point vector.",err,error,*999)
            numberOfXi = dependentBasis%NUMBER_OF_XI

            SELECT CASE(smoothingType)
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
              materialsField=>equations%interpolation%materialsField
              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,equations%interpolation% &
                & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

            !===============================
            ! G a u s s   P o i n t   F i t
            !===============================
            !Loop over Gauss points
            DO gaussPointIdx=1,quadratureScheme%NUMBER_OF_GAUSS
              ! Interpolate fields
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & equations%interpolation%dependentInterpPoint(dependentVariableType)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                & equations%interpolation%independentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(geometricBasis%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

              !Get fitting data from interpolated fields
              DO componentIdx=1,numberOfDataComponents
                dataPointVector(componentIdx)=equations%interpolation%independentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                  & values(componentIdx,NO_PART_DERIV)
                dataPointWeight(componentIdx)=equations%interpolation%independentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                  & values(componentIdx,NO_PART_DERIV)
              ENDDO !componentIdx

              SELECT CASE(smoothingType)
              CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                !Do nothing
              CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
                !Get Sobolev smoothing data from interpolated fields
                CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
                  & equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                tauParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(1,NO_PART_DERIV)
                kappaParam=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%values(2,NO_PART_DERIV)
              CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE DEFAULT
                localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT

              jacobianGaussWeight=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                & quadratureScheme%GAUSS_WEIGHTS(gaussPointIdx)

              dependentParameterRowIdx=0
              !Loop over element rows
              DO dependentComponentRowIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                meshComponentRow=dependentVariable%components(dependentComponentRowIdx)%MESH_COMPONENT_NUMBER
                dependentBasisRow=>dependentField%decomposition%domain(meshComponentRow)%ptr% &
                  & topology%elements%elements(elementNumber)%basis
                quadratureSchemeRow=>dependentBasisRow%quadrature%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                DO dependentElementParameterRowIdx=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
                  dependentParameterRowIdx=dependentParameterRowIdx+1
                  dependentParameterColumnIdx=0
                  phiM=quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,NO_PART_DERIV,gaussPointIdx)
                  IF(equationsMatrix%updateMatrix) THEN
                    !Loop over element columns
                    DO dependentComponentColumnIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                      meshComponentColumn=dependentVariable%components(dependentComponentColumnIdx)%MESH_COMPONENT_NUMBER
                      dependentBasisColumn=>dependentField%decomposition%domain(meshComponentColumn)%ptr% &
                        & topology%elements%elements(elementNumber)%basis
                      quadratureSchemeColumn=>dependentBasisColumn%quadrature%QUADRATURE_SCHEME_MAP( &
                        & BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO dependentElementParameterColumnIdx=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                        dependentParameterColumnIdx=dependentParameterColumnIdx+1
                        !Treat each component as separate and independent so only calculate the diagonal blocks
                        IF(dependentComponentColumnIdx==dependentComponentRowIdx) THEN
                          phiN=quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,NO_PART_DERIV, &
                            & gaussPointIdx)
                          sum=phiM*phiN*dataPointWeight(dependentComponentRowIdx)
                          SELECT CASE(smoothingType)
                          CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
                            !Do nothing
                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING)
                            !Calculate Sobolev surface tension and curvature smoothing terms
                            tension = tauParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1, &
                              & gaussPointIdx)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1, &
                              & gaussPointIdx))
                            curvature = kappaParam*2.0_DP* ( &
                              & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S1, &
                              & gaussPointIdx)* &
                              & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S1, &
                              & gaussPointIdx))
                            IF(numberOfXi > 1) THEN
                              tension = tension + tauParam*2.0_DP* ( &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2, &
                                & gaussPointIdx))
                              curvature = curvature + kappaParam*2.0_DP* ( &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2_S2, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2_S2, &
                                & gaussPointIdx) + &
                                & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S2, &
                                & gaussPointIdx)* &
                                & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S2, &
                                & gaussPointIdx))
                              IF(numberOfXi > 2) THEN
                                tension = tension + tauParam*2.0_DP* ( &
                                  & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S3, &
                                  & gaussPointIdx)* &
                                  & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S3, &
                                  & gaussPointIdx))
                                curvature = curvature + kappaParam*2.0_DP* ( &
                                  & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S3_S3, &
                                  & gaussPointIdx)* &
                                  & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S3_S3, &
                                  & gaussPointIdx)+ &
                                  & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S1_S3, &
                                  & gaussPointIdx)* &
                                  & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S1_S3, &
                                  & gaussPointIdx)+ &
                                  & quadratureSchemeRow%GAUSS_BASIS_FNS(dependentElementParameterRowIdx,PART_DERIV_S2_S3, &
                                  & gaussPointIdx)* &
                                  & quadratureSchemeColumn%GAUSS_BASIS_FNS(dependentElementParameterColumnIdx,PART_DERIV_S2_S3, &
                                  & gaussPointIdx))
                              ENDIF ! 3D
                            ENDIF ! 2 or 3D
                            sum = sum + (tension + curvature) * jacobianGaussWeight
                          CASE(EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
                            CALL FlagError("Not implemented.",err,error,*999)
                          CASE DEFAULT
                            localError="The fitting smoothing type of "//TRIM(NumberToVString(smoothingType,"*",err,error))// &
                              & " is invalid."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT

                          equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)= &
                            equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)+sum
                        ENDIF
                      ENDDO !dependentElementParameterColumnIdx
                    ENDDO !dependentComponentColumnIdx
                  ENDIF
                  IF(rhsVector%updateVector) THEN
                    rhsVector%elementVector%vector(dependentParameterRowIdx)= &
                      & rhsVector%elementVector%vector(dependentParameterRowIdx) + &
                      & phiM*dataPointVector(dependentComponentRowIdx)*dataPointWeight(dependentComponentRowIdx)
                  ENDIF
                ENDDO !dependentElementParameterRowIdx
              ENDDO !dependentComponentRowIdx
            ENDDO !gaussPointIdx

          CASE DEFAULT
            localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is not valid for a Gauss fitting equations set class."
            CALL FlagError(localError,err,error,*999)
          END SELECT


        CASE DEFAULT
          localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
            & " is not valid for a fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT

        !Scale factor adjustment
        IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
          CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,equations%interpolation% &
            & dependentInterpParameters(dependentVariableType)%ptr,err,error,*999)
          dependentParameterRowIdx=0
          !Loop over element rows
          DO dependentComponentRowIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
            meshComponentRow=dependentVariable%components(dependentComponentRowIdx)%MESH_COMPONENT_NUMBER
            dependentBasisRow=>dependentField%decomposition%domain(meshComponentRow)%ptr% &
              & topology%elements%elements(elementNumber)%basis
            DO dependentElementParameterRowIdx=1,dependentBasisRow%NUMBER_OF_ELEMENT_PARAMETERS
              dependentParameterRowIdx=dependentParameterRowIdx+1
              dependentParameterColumnIdx=0
              IF(equationsMatrix%updateMatrix) THEN
                !Loop over element columns
                DO dependentComponentColumnIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
                  meshComponentColumn=dependentVariable%components(dependentComponentColumnIdx)%MESH_COMPONENT_NUMBER
                  dependentBasisColumn=>dependentField%decomposition%domain(meshComponentColumn)%ptr% &
                    & topology%elements%elements(elementNumber)%basis
                  DO dependentElementParameterColumnIdx=1,dependentBasisColumn%NUMBER_OF_ELEMENT_PARAMETERS
                    dependentParameterColumnIdx=dependentParameterColumnIdx+1
                    equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)= &
                      & equationsMatrix%elementMatrix%matrix(dependentParameterRowIdx,dependentParameterColumnIdx)* &
                      & equations%interpolation%dependentInterpParameters(dependentVariableType)%ptr% &
                      & SCALE_FACTORS(dependentElementParameterRowIdx,dependentComponentRowIdx)* &
                      & equations%interpolation%dependentInterpParameters(dependentVariableType)%ptr% &
                      & SCALE_FACTORS(dependentElementParameterColumnIdx,dependentComponentColumnIdx)
                  ENDDO !dependentElementParameterColumnIdx
                ENDDO !dependentComponentColumnIdx
              ENDIF
              IF(rhsVector%updateVector) THEN
                rhsVector%elementVector%vector(dependentParameterRowIdx)= &
                  & rhsVector%elementVector%vector(dependentParameterRowIdx)* &
                  & equations%interpolation%dependentInterpParameters(dependentVariableType)%ptr% &
                  & SCALE_FACTORS(dependentElementParameterRowIdx,dependentComponentRowIdx)
              ENDIF
            ENDDO !dependentElementParameterRowIdx
          ENDDO !dependentComponentRowIdx
        ENDIF

      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Fitting_FiniteElementCalculate",err,error)
    RETURN 1

  END SUBROUTINE Fitting_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the update-materials Galerkin projection.
  SUBROUTINE FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,GEOMETRIC_COMPONENT_NUMBER,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,numberOfDimensions,I,MATERIAL_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
!     TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,geometricField
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: localError

    ENTERS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)/=4) THEN
        CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
          & err,error,*999)
      END IF
      IF(equationsSet%specification(3)==EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE.OR. &
        & equationsSet%specification(3)==EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE) THEN
        SELECT CASE(equationsSetSetup%SETUP_TYPE)

        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL Field_LabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL Field_MeshDecompositionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)


              CALL Field_NumberOfVariablesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & err,error,*999)

              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              !component 1: dependent porosity variable, component 2: dependent permeability variable
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)

              DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              END DO


              SELECT CASE(equationsSet%SOLUTION_METHOD)

              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO I=1,DEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
                !Other solutions not defined yet
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
            !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)

              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)

              !component 1: dependent porosity variable, component 2: dependent permeability variable
              DEPENDENT_FIELD_NUMBER_OF_COMPONENTS=2
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & DEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created independent field
              !start field creation with name 'INDEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !define new created field to be independent
              CALL Field_DependentTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_INDEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_MeshDecompositionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% &
                & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
              !set number of variables to 1 (1 for U)
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & [FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=numberOfDimensions
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              !Default to the geometric interpolation setup
              DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              END DO
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=numberOfDimensions
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
            ENDIF
          !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
            ENDIF
            CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
            MATERIAL_FIELD_NUMBER_OF_VARIABLES=1
            MATERIAL_FIELD_NUMBER_OF_COMPONENTS=2

            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                !start field creation with name 'MATERIAL_FIELD'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet% &
                  & MATERIALS%MATERIALS_FIELD,err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%MATERIALS%MATERIALS_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,equationsSet%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL Field_NumberOfVariablesSet(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                  & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,"Fitting Materials", &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                DO I = 1, MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL Field_CreateFinish(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
! ! !         CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
! ! !           SELECT CASE(equationsSetSetup%ACTION_TYPE)
! ! !           CASE(EQUATIONS_SET_SETUP_START_ACTION)
! ! !             IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               DEPENDENT_FIELD=>equationsSet%DEPENDENT%DEPENDENT_FIELD
! ! !               IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !                 geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
! ! !                 IF(ASSOCIATED(geometricField)) THEN
! ! !                   CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
! ! !                   SELECT CASE(equationsSetSetup%ANALYTIC_FUNCTION_TYPE)
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                     !Check that we are in 2D
! ! !                     IF(numberOfDimensions/=2) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                     !Check that we are in 2D
! ! !                     IF(numberOfDimensions/=2) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_2
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                     !Check that we are in 3D
! ! !                     IF(numberOfDimensions/=4) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                     !Check that we are in 3D
! ! !                     IF(numberOfDimensions/=4) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_2
! ! !                   CASE DEFAULT
! ! !                     localError="The specified analytic function type of "// &
! ! !                       & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                       & " is invalid for a moving mesh Galerkin projection."
! ! !                     CALL FlagError(localError,err,error,*999)
! ! !                   END SELECT
! ! !                 ELSE
! ! !                   CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
! ! !                 ENDIF
! ! !              ELSE
! ! !                 CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
! ! !             IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
! ! !               ANALYTIC_FIELD=>equationsSet%ANALYTIC%ANALYTIC_FIELD
! ! !               IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
! ! !                 IF(equationsSet%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
! ! !                   CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
! ! !                 ENDIF
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FlagError("Equations set analytic is not associated.",err,error,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
! ! !             IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
! ! !                 IF(equationsSet%ANALYTIC%ANALYTIC_FINISHED) THEN
! ! !                   CALL FITTING_ANALYTIC_CALCULATE(equationsSet,err,error,*999)
! ! !                 ELSE
! ! !                   CALL FlagError("Equations set analtyic has not been finished.",err,error,*999)
! ! !                 ENDIF
! ! !               ELSE
! ! !                 CALL FlagError("Equations set analtyic is not associated.",err,error,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FlagError("Equations set dependent has not been finished.",err,error,*999)
! ! !             ENDIF
! ! !           CASE DEFAULT
! ! !             localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
! ! !               & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
! ! !               & " is invalid for an update-materials Galerkin projection."
! ! !             CALL FlagError(localError,err,error,*999)
! ! !           END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(equationsSet%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   c a s e   d e f a u l t
        !-----------------------------------------------------------------
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for an update-materials Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " does not equal an update-materials Galerkin projection subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP")
    RETURN
999 ERRORSEXITS("FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Galerkin projection type of a data fitting equations set class.
  SUBROUTINE Fitting_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup a Galerkin projection on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSetup",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE
        IF(SIZE(equationsSet%specification,1)/=4) THEN
          CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
            & err,error,*999)
        ENDIF
      ENDIF
      SELECT CASE(equationsSet%specification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_DATA_POINT_FITTING_SUBTYPE)
          CALL Fitting_EquationsSetDataSetup(equationsSet,equationsSetSetup,err,error,*999)
        CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE)
          CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(equationsSet,equationsSetSetup,err,error,*999)
        CASE(EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
          CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(equationsSet,equationsSetSetup,err,error,*999)
        CASE(EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
          CALL FITTING_EQUATIONS_SET_VECTORDATA_SETUP(equationsSet,equationsSetSetup,err,error,*999)
        CASE(EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE)
          CALL FITTING_EQUATIONS_SET_MAT_PROPERTIES_SETUP(equationsSet,equationsSetSetup,err,error,*999)
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a data fitting equation set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_GAUSS_POINT_FITTING_SUBTYPE)
          CALL Fitting_EquationsSetGaussSetup(equationsSet,equationsSetSetup,err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a Gauss fitting equation set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set type "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fitting equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Galerkin projection type of an data fitting equations set class.
  SUBROUTINE Fitting_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE
        IF(SIZE(equationsSet%specification,1)/=4) &
          & CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
          & err,error,*999)
      ENDIF
      SELECT CASE(equationsSet%specification(2))
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_DATA_POINT_FITTING_SUBTYPE,EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE,EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE,EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
          localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a data fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_GAUSS_POINT_FITTING_SUBTYPE)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            equationsSet%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
          localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is not valid for a Gauss fitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " is not valid for a fitting equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetSolutionMethodSet",err,error)
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a data fitting equation set class.
  SUBROUTINE Fitting_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetType,equationsSetSubtype,equationsSetSmoothing
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=4) THEN
        CALL FlagError("Equations set specification must have four entries for a fitting class equations set.", &
          & err,error,*999)
      ENDIF
      equationsSetType=specification(2)
      SELECT CASE(equationsSetType)
      CASE(EQUATIONS_SET_DATA_FITTING_EQUATION_TYPE)
        equationsSetSubtype=specification(3)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_DATA_POINT_FITTING_SUBTYPE, &
          & EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_MAT_PROPERTIES_INRIA_MODEL_DATA_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
          & EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
        CASE(EQUATIONS_SET_GENERALISED_DATA_FITTING_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The third equations set specifiction of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a data fitting equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE)
        equationsSetSubtype=specification(3)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_GAUSS_POINT_FITTING_SUBTYPE)
          !OK
        CASE DEFAULT
          localError="The third equations set specifiction of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a Gauss fitting equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second equations set specification of "//TRIM(NumberToVstring(equationsSetType,"*",err,error))// &
          & " is not valid for a data fitting equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      equationsSetSmoothing=specification(4)
      SELECT CASE(equationsSetSmoothing)
      CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING,EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING, &
        & EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING,EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
        !ok
      CASE DEFAULT
        localError="The fourth equations set specification of "//TRIM(NumberToVstring(equationsSetSmoothing,"*",err,error))// &
          & " is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(4),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:4)=[EQUATIONS_SET_FITTING_CLASS,equationsSetType,equationsSetSubtype,equationsSetSmoothing]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Fitting_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Fitting_EquationsSetSpecificationSet",err,error)
    EXITS("Fitting_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Gauss fitting equations set.
  SUBROUTINE Fitting_EquationsSetGaussSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,geometricMeshComponent,geometricScalingType,numberOfComponents,numberOfComponents2, &
      & numberOfDependentComponents,numberOfIndependentComponents
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equationsMaterials
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetGaussSetup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE
        IF(SIZE(equationsSet%specification,1)/=4) THEN
          CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
            & err,error,*999)
        ENDIF
      ENDIF
      IF(equationsSet%specification(2)==EQUATIONS_SET_GAUSS_FITTING_EQUATION_TYPE) THEN
        SELECT CASE(equationsSet%specification(3))
        CASE(EQUATIONS_SET_GAUSS_POINT_FITTING_SUBTYPE)
          SELECT CASE(equationsSetSetup%SETUP_TYPE)
          CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
            !-----------------------------------------------------------------
            ! s o l u t i o n   m e t h o d
            !-----------------------------------------------------------------
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Gauss point fitting equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
            !-----------------------------------------------------------------
            ! g e o m e t r y   f i e l d
            !-----------------------------------------------------------------
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
            !-----------------------------------------------------------------
            ! S o u r c e   f i e l d
            !-----------------------------------------------------------------
            ! Do nothing
          CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
            !-----------------------------------------------------------------
            ! D e p e n d e n t   f i e l d
            !-----------------------------------------------------------------
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              !Set start action
              IF(equationsSet%dependent%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%region, &
                  & equationsSet%dependent%DEPENDENT_FIELD,err,error,*999)
                CALL Field_LabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                !start creation of a new field
                CALL Field_TypeSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !define new created field to be dependent
                CALL Field_DependentTypeSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL Field_MeshDecompositionGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,geometricDecomposition,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,equationsSet% &
                  & geometry%GEOMETRIC_FIELD,err,error,*999)
                !set number of variables to 2 (U)
!!TODO: We only really need 1 variable here as there is no real RHS variable. However, until the mapping system is redone so that we don't have a RHS we will create 2 variables.
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
                CALL Field_VariableLabelSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"delUdeln", &
                   err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                !Set the number of components.
                !If the independent field has been defined use that number of components
                IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                  IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                    CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfComponents,err,error,*999)
                  ELSE
                    numberOfComponents=1
                  ENDIF
                ELSE
                  numberOfComponents=1
                ENDIF
                CALL Field_NumberOfComponentsSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfComponents,err,error,*999)
                CALL Field_NumberOfComponentsSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & numberOfComponents,err,error,*999)
                !Default to the first geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & geometricMeshComponent,err,error,*999)
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentMeshComponentSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,numberOfComponents
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                  CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsSet%dependent%DEPENDENT_FIELD,geometricScalingType,err,error,*999)
                  !Other solutions not defined yet
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
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                !If the independent field has been defined checkk the number of components is the same
                IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                  IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                    CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfIndependentComponents,err,error,*999)
                    IF(numberOfComponents /= numberOfIndependentComponents) THEN
                      localError="The number of components for the specified dependent field of "// &
                        & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                        & " does not match the number of components for the independent field of "// &
                        & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF
                ENDIF
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents2, &
                  & err,error,*999)
                IF(numberOfComponents2/=numberOfComponents) THEN
                  localError="The number of components in the dependent field for variable type "// &
                     & TRIM(NumberToVstring(FIELD_DELUDELN_VARIABLE_TYPE,"*",err,error))//" of "// &
                     & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))// &
                     & " does not match the number of components for variable type "// &
                     & TRIM(NumberToVstring(FIELD_U_VARIABLE_TYPE,"*",err,error))//" of "// &
                     & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,numberOfComponents
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
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Check that we have the same number of components as the independent field
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents2, &
                  & err,error,*999)
                IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                  IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                    CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfIndependentComponents,err,error,*999)
                    IF(numberOfComponents /= numberOfIndependentComponents) THEN
                      localError="The number of components for the specified dependent field of "// &
                        & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                        & " does not match the number of components for the independentt field of "// &
                        & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF
                ENDIF
                IF(numberOfComponents2/=numberOfComponents) THEN
                  localError="The number of components in the dependent field for variable type "// &
                    & TRIM(NumberToVstring(FIELD_DELUDELN_VARIABLE_TYPE,"*",err,error))//" of "// &
                    & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))// &
                    & " does not match the number of components for variable type "// &
                    & TRIM(NumberToVstring(FIELD_U_VARIABLE_TYPE,"*",err,error))//" of "// &
                    & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Finish creating the field
                CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for an update-materials Galerkin projection"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
            !-----------------------------------------------------------------
            !   m a t e r i a l   f i e l d
            !-----------------------------------------------------------------
            SELECT CASE(equationsSet%specification(4))
            CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
              !Do nothing
            CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING,EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
              SELECT CASE(equationsSetSetup%ACTION_TYPE)
              CASE(EQUATIONS_SET_SETUP_START_ACTION)
                equationsMaterials=>equationsSet%materials
                IF(ASSOCIATED(equationsMaterials)) THEN
                  IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                    !Create the auto created materials field
                    CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%region,equationsSet% &
                      & materials%MATERIALS_FIELD,err,error,*999)
                    CALL Field_TypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                    CALL Field_DependentTypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                      & err,error,*999)
                    CALL Field_MeshDecompositionGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricDecomposition, &
                      & err,error,*999)
                    !apply decomposition rule found on new created field
                    CALL Field_MeshDecompositionSetAndLock(equationsSet%materials%MATERIALS_FIELD,geometricDecomposition, &
                      & err,error,*999)
                    !point new field to geometric field
                    CALL Field_GeometricFieldSetAndLock(equationsMaterials%MATERIALS_FIELD,equationsSet%geometry% &
                      & GEOMETRIC_FIELD,err,error,*999)
                    CALL Field_NumberOfVariablesSet(equationsMaterials%MATERIALS_FIELD,1,err,error,*999)
                    CALL Field_VariableTypesSetAndLock(equationsMaterials%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                      & err,error,*999)
                    CALL Field_DimensionSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                    CALL Field_DataTypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_DP_TYPE,err,error,*999)
                    !Sobolev smoothing material parameters- tau and kappa
                    CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 2,err,error,*999)
                    CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentMeshComponentSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 2,geometricMeshComponent,err,error,*999)
                    CALL Field_ComponentInterpolationSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    !Default the field scaling to that of the geometric field
                    CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                    CALL Field_ScalingTypeSet(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                  ELSE
                    !Check the user specified field
                    CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
                    CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                    CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                    CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                    CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                      & err,error,*999)
                    CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                    CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,2,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set materials is not associated.",err,error,*999)
                END IF
              CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                equationsMaterials=>equationsSet%materials
                IF(ASSOCIATED(equationsMaterials)) THEN
                  IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                    !Finish creating the materials field
                    CALL Field_CreateFinish(equationsMaterials%MATERIALS_FIELD,err,error,*999)
                    !Set the default values for the materials field
                    CALL Field_ComponentValuesInitialise(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                    CALL Field_ComponentValuesInitialise(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set materials is not associated.",err,error,*999)
                ENDIF
              CASE DEFAULT
                localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for an update-materials Galerkin projection."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The fitting smoothing type of "//TRIM(NumberToVString(equationsSet%specification(4),"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT

          CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
            !-----------------------------------------------------------------
            ! I n d e p e n d e n t   t y p e
            ! (this field holds the data point based field of vectors to map to the dependent field)
            !-----------------------------------------------------------------
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
              !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%region, &
                  & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL Field_LabelSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error, &
                  & *999)
                !define new created field to be independent
                CALL Field_DependentTypeSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL Field_MeshDecompositionGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,geometricDecomposition, &
                  & err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% &
                  & geometry%GEOMETRIC_FIELD,err,error,*999)
                !Create two variables: U for data and V for weights
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & 2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricMeshComponent,err,error,*999)
                !If the dependent field has been created then use that number of components
                IF(ASSOCIATED(equationsSet%DEPENDENT%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfComponents,err,error,*999)
                ELSE
                  numberOfComponents=1
                ENDIF
                ! U Variable: data points
                CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                !Default to the geometric interpolation setup
                ! V Variable: data point weights
                CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                 CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_V_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                  !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx = 1,numberOfComponents
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                  CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsSet%independent%INDEPENDENT_FIELD,geometricScalingType,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                ! U (vector) variable
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                  & numberOfComponents,err,error,*999)
                !If the dependent field has been defined use that number of components
                IF(ASSOCIATED(equationsSet%dependent%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDependentComponents,err,error,*999)
                  IF(numberOfComponents /= numberOfDependentComponents) THEN
                    localError="The number of components for the specified independent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the dependent field of "// &
                      & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDIF
                ! V (weight) variable
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,numberOfComponents
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Check that we have the same number of components as the dependent field
                IF(ASSOCIATED(equationsSet%dependent%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDependentComponents,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                  IF(numberOfComponents /= numberOfDependentComponents) THEN
                    localError="The number of components for the specified independent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the dependentt field of "// &
                      & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDIF
                !Specify finish action
                CALL Field_CreateFinish(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a fitting equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
            !-----------------------------------------------------------------
            !   e q u a t i o n s   t y p e
            !-----------------------------------------------------------------
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%dependent%DEPENDENT_FINISHED) THEN
                IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                  IF(equationsSet%INDEPENDENT%INDEPENDENT_FINISHED) THEN
                    CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
                    CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
                    CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
                  ELSE
                    CALL FlagError("Equations set independent field has not been finished.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set independent field has not been finished.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations creation
                CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a fitting equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a fitting equations set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set type of "//TRIM(NumberToVString(equationsSet%specification(2),"*",err,error))// &
          & " does not equal a Gauss fitting type."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_EquationsSetGaussSetup")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetGaussSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_EquationsSetGaussSetup

  !
  !================================================================================================================================
  !

  !>Sets up the standard data point fitting equations set.
  SUBROUTINE Fitting_EquationsSetDataSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,geometricMeshComponent,geometricScalingType,numberOfComponents,numberOfComponents2, &
      & numberOfDependentComponents,numberOfIndependentComponents
    TYPE(DECOMPOSITION_TYPE), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: equationsMaterials
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_EquationsSetDataSetup",err,error,*999)

    NULLIFY(equations)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(equationsMaterials)
    NULLIFY(geometricDecomposition)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      IF(SIZE(equationsSet%specification,1)/=4) &
        & CALL FlagError("Equations set specification must have four entries for a fitting type equations set.",err,error,*999)
      IF(equationsSet%specification(3)==EQUATIONS_SET_DATA_POINT_FITTING_SUBTYPE) THEN
        SELECT CASE(equationsSetSetup%SETUP_TYPE)

        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          !-----------------------------------------------------------------
          ! s o l u t i o n   m e t h o d
          !-----------------------------------------------------------------
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !-----------------------------------------------------------------
          ! g e o m e t r y   f i e l d
          !-----------------------------------------------------------------
          !Do nothing

        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          !-----------------------------------------------------------------
          ! d e p e n d e n t   f i e l d
          !-----------------------------------------------------------------
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL Field_LabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
              CALL Field_MeshDecompositionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricDecomposition, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              !Set the number of components.
              !If the independent field has been defined use that number of components
              IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfComponents,err,error,*999)
                ELSE
                  numberOfComponents=1
                ENDIF
              ELSE
                numberOfComponents=1
              ENDIF
              CALL Field_NumberOfComponentsSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfComponents,err,error,*999)
              CALL Field_NumberOfComponentsSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfComponents,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & geometricMeshComponent,err,error,*999)
              DO componentIdx=1,numberOfComponents
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & componentIdx,geometricMeshComponent,err,error,*999)
              ENDDO !componentIdx
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !componentIdx
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%dependent%DEPENDENT_FIELD,geometricScalingType,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
              !If the independent field has been defined check the number of components is the same
              IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfIndependentComponents,err,error,*999)
                  IF(numberOfComponents /= numberOfIndependentComponents) THEN
                    localError="The number of components for the specified dependent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the independent field of "// &
                        & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDIF
              ENDIF
              CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents2, &
                & err,error,*999)
              IF(numberOfComponents2/=numberOfComponents) THEN
                localError="The number of components in the independent field for variable type "// &
                  & TRIM(NumberToVstring(FIELD_DELUDELN_VARIABLE_TYPE,"*",err,error))//" of "// &
                  & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))// &
                  & " does not match the number of components for variable type "// &
                  & TRIM(NumberToVstring(FIELD_U_VARIABLE_TYPE,"*",err,error))//" of "// &
                     & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,componentIdx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !numberOfComponents
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Check that we have the same number of components as the independent field
              IF(ASSOCIATED(equationsSet%dependent%DEPENDENT_FIELD)) THEN
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfComponents2, &
                  & err,error,*999)
                IF(ASSOCIATED(equationsSet%INDEPENDENT)) THEN
                  IF(ASSOCIATED(equationsSet%INDEPENDENT%INDEPENDENT_FIELD)) THEN
                    CALL Field_NumberOfComponentsGet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfIndependentComponents,err,error,*999)
                    IF(numberOfComponents /= numberOfIndependentComponents) THEN
                      localError="The number of components for the specified dependent field of "// &
                        & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                        & " does not match the number of components for the independentt field of "// &
                        & TRIM(NumberToVString(numberOfIndependentComponents,"*",err,error))//"."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF
                ENDIF
                IF(numberOfComponents2/=numberOfComponents) THEN
                  localError="The number of components in the dependent field for variable type "// &
                    & TRIM(NumberToVstring(FIELD_DELUDELN_VARIABLE_TYPE,"*",err,error))//" of "// &
                    & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))// &
                    & " does not match the number of components for variable type "// &
                    & TRIM(NumberToVstring(FIELD_U_VARIABLE_TYPE,"*",err,error))//" of "// &
                    & TRIM(NumberToVstring(numberOfComponents2,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Finish creating the field
                CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              ELSE
                CALL FlagError("The equations set dependent field is not associated.",err,error,*999)
              ENDIF
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          !-----------------------------------------------------------------
          !   m a t e r i a l   f i e l d
          !-----------------------------------------------------------------
          SELECT CASE(equationsSet%specification(4))
          CASE(EQUATIONS_SET_FITTING_NO_SMOOTHING)
            !Do nothing
          CASE(EQUATIONS_SET_FITTING_SOBOLEV_VALUE_SMOOTHING,EQUATIONS_SET_FITTING_SOBOLEV_DIFFERENCE_SMOOTHING)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              equationsMaterials=>equationsSet%materials
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Create the auto created materials field
                  CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%region,equationsSet% &
                    & materials%MATERIALS_FIELD,err,error,*999)
                  CALL Field_TypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL Field_DependentTypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                    & err,error,*999)
                  CALL Field_MeshDecompositionGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricDecomposition, &
                    & err,error,*999)
                  !apply decomposition rule found on new created field
                  CALL Field_MeshDecompositionSetAndLock(equationsSet%materials%MATERIALS_FIELD,geometricDecomposition, &
                    & err,error,*999)
                  !point new field to geometric field
                  CALL Field_GeometricFieldSetAndLock(equationsMaterials%MATERIALS_FIELD,equationsSet%geometry% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL Field_NumberOfVariablesSet(equationsMaterials%MATERIALS_FIELD,1,err,error,*999)
                  CALL Field_VariableTypesSetAndLock(equationsMaterials%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL Field_DimensionSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL Field_DataTypeSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  !Sobolev smoothing material parameters- tau and kappa
                  CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,err,error,*999)
                  CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,geometricMeshComponent,err,error,*999)
                  CALL Field_ComponentInterpolationSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSet(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  !Default the field scaling to that of the geometric field
                  CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsMaterials%MATERIALS_FIELD,geometricScalingType,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
                  CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,2,err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              END IF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              equationsMaterials=>equationsSet%materials
              IF(ASSOCIATED(equationsMaterials)) THEN
                IF(equationsMaterials%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL Field_CreateFinish(equationsMaterials%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials field
                  CALL Field_ComponentValuesInitialise(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                  CALL Field_ComponentValuesInitialise(equationsMaterials%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set materials is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for an update-materials Galerkin projection."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_FITTING_STRAIN_ENERGY_SMOOTHING)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The fitting smoothing type of "//TRIM(NumberToVString(equationsSet%specification(4),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT

          CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
            !-----------------------------------------------------------------
            ! I n d e p e n d e n t   t y p e
            !
            ! (this field holds the data point based field of vectors to map to the dependent field)
            !-----------------------------------------------------------------
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
              !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created independent field
                !start field creation with name 'INDEPENDENT_FIELD'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%region, &
                  & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !start creation of a new field
                CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL Field_LabelSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error, &
                  & *999)
                !define new created field to be independent
                CALL Field_DependentTypeSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL Field_MeshDecompositionGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricDecomposition,err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,geometricDecomposition, &
                  & err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% &
                  & geometry%GEOMETRIC_FIELD,err,error,*999)
                !Create two variables: U for data and V for weights
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & 2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,geometricMeshComponent,err,error,*999)
                !If the dependent field has been created then use that number of components
                IF(ASSOCIATED(equationsSet%DEPENDENT%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfComponents,err,error,*999)
                ELSE
                  numberOfComponents=1
                ENDIF
                ! U Variable: data points
                CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                !Default to the geometric interpolation setup
                ! V Variable: data point weights
                CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                !Default to the geometric interpolation setup
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                 CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                   & FIELD_V_VARIABLE_TYPE,componentIdx,geometricMeshComponent,err,error,*999)
                ENDDO !componentIdx
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                  !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx = 1,numberOfComponents
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                      & FIELD_V_VARIABLE_TYPE,componentIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                  CALL Field_ScalingTypeGet(equationsSet%geometry%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                  CALL Field_ScalingTypeSet(equationsSet%independent%INDEPENDENT_FIELD,geometricScalingType,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
                ! U (vector) variable
                CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE, &
                  & numberOfComponents,err,error,*999)
                !If the dependent field has been defined use that number of components
                IF(ASSOCIATED(equationsSet%dependent%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDependentComponents,err,error,*999)
                  IF(numberOfComponents /= numberOfDependentComponents) THEN
                    localError="The number of components for the specified independent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the dependent field of "// &
                      & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDIF
                ! V (weight) variable
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  DO componentIdx=1,numberOfComponents
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,componentIdx, &
                      & FIELD_GAUSS_POINT_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !componentIdx
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%independent%INDEPENDENT_FIELD_AUTO_CREATED) THEN
                !Check that we have the same number of components as the dependent field
                IF(ASSOCIATED(equationsSet%dependent%DEPENDENT_FIELD)) THEN
                  CALL Field_NumberOfComponentsGet(equationsSet%dependent%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & numberOfDependentComponents,err,error,*999)
                  CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfComponents,err,error,*999)
                  IF(numberOfComponents /= numberOfDependentComponents) THEN
                    localError="The number of components for the specified independent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the dependentt field of "// &
                      & TRIM(NumberToVString(numberOfDependentComponents,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                  CALL Field_NumberOfComponentsGet(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,numberOfComponents2,err,error,*999)
                  IF(numberOfComponents /= numberOfComponents2) THEN
                    localError="The number of components for the U variable of the specified independent field of "// &
                      & TRIM(NumberToVString(numberOfComponents,"*",err,error))// &
                      & " does not match the number of components for the V variable of the field of "// &
                      & TRIM(NumberToVString(numberOfComponents2,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*999)
                  ENDIF
                ENDIF
                !Specify finish action
                CALL Field_CreateFinish(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
                !Initialise the weights to 1.0
                DO componentIdx=1,numberOfComponents
                  CALL Field_ComponentValuesInitialise(equationsSet%independent%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
                ENDDO !componentIdx
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a fitting equations set."
              CALL FlagError(localError,err,error,*999)
            END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          !-----------------------------------------------------------------
          !   s o u r c e   t y p e
          !-----------------------------------------------------------------
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          !-----------------------------------------------------------------
          !   e q u a t i o n s   t y p e
          !-----------------------------------------------------------------
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(equationsSet%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(equations%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        CASE DEFAULT
          !-----------------------------------------------------------------
          !   c a s e   d e f a u l t
          !-----------------------------------------------------------------
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " does not equal a standard Galerkin projection subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_EquationsSetDataSetup")
    RETURN
999 ERRORSEXITS("Fitting_EquationsSetDataSetup",err,error)
    RETURN 1
  END SUBROUTINE Fitting_EquationsSetDataSetup

  !
  !================================================================================================================================
  !

  !>Sets up the vector data Galerkin projection.
  SUBROUTINE FITTING_EQUATIONS_SET_VECTORDATA_SETUP(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,geometricScalingType,GEOMETRIC_COMPONENT_NUMBER
    INTEGER(INTG) :: numberOfDimensions,I !,MATERIAL_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,INDEPENDENT_FIELD_NUMBER_OF_VARIABLES
    INTEGER(INTG) :: dependentFieldNumberOfVariables
    INTEGER(INTG) :: dimensionIdx
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
!     TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: localError

    ENTERS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP",err,error,*999)

    NULLIFY(BOUNDARY_CONDITIONS)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)/=4) THEN
        CALL FlagError("Equations set specification must have four entries for a fitting type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(equationsSet%specification(3))
      CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing

        !-----------------------------------------------------------------
        ! d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL Field_LabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL Field_MeshDecompositionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension and one for pressure
              CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              IF(equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                & equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
! !                 CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
! !                   & 1,err,error,*999)

                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
! !                 CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
! !                   & 1,err,error,*999)

!                 DO I=1,1
                DO I=1,numberOfDimensions
                  !Default to the geometric interpolation setup
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO
              ELSE IF(equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE.OR. &
                equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE) THEN
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions+1,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & numberOfDimensions+1,err,error,*999)
                DO I=1,numberOfDimensions+1
                  !Default to the geometric interpolation setup
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                END DO
              ENDIF
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                IF(equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE) THEN
!                   DO I=1,numberOfDimensions
                  DO I=1,1
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                ELSE IF(equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,numberOfDimensions+1
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                ENDIF
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType, &
                  & err,error,*999)
                !Other solutions not defined yet
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
            !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,2,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                IF(equationsSet%specification(3)==EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,numberOfDimensions
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                ELSEIF(equationsSet%specification(3)==EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
                  DO I=1,numberOfDimensions+1
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_U_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                      & FIELD_DELUDELN_VARIABLE_TYPE,I,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  END DO
                ENDIF
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType, &
                  & err,error,*999)
                !Other solutions not defined yet
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        ! I N d e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created independent field
              !start field creation with name 'INDEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !define new created field to be independent
              CALL Field_DependentTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_INDEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_MeshDecompositionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% &
                & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
              !set number of variables to 1 (1 for U)
              INDEPENDENT_FIELD_NUMBER_OF_VARIABLES=1
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & INDEPENDENT_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & [FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=numberOfDimensions
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              !Default to the geometric interpolation setup
              DO I=1,INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS
                CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,I,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              END DO
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
              CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS=numberOfDimensions
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & INDEPENDENT_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
            ENDIF
          !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
            ENDIF
            CALL Field_ParameterSetCreate(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                !start field creation with name 'MATERIAL_FIELD'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet% &
                  & MATERIALS%MATERIALS_FIELD,err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%MATERIALS%MATERIALS_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,equationsSet%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL Field_NumberOfVariablesSet(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                  & 1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL Field_CreateFinish(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT

        !-----------------------------------------------------------------
        !   s o u r c e   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(equationsSet%specification(3))
          CASE(EQUATIONS_SET_VECTOR_DATA_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_FITTING_SUBTYPE, &
            & EQUATIONS_SET_VECTOR_DATA_PRE_FITTING_SUBTYPE,EQUATIONS_SET_DIVFREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
            SELECT CASE(equationsSetSetup%ACTION_TYPE)
              !Set start action
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(equationsSet%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                !start field creation with name 'sourceField'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                  & equationsSet%SOURCE%SOURCE_FIELD,err,error,*999)
                !start creation of a new field
                CALL Field_TypeSetAndLock(equationsSet%SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !label the field
                CALL Field_LabelSetAndLock(equationsSet%SOURCE%SOURCE_FIELD,"Source Field",err,error, &
                  & *999)
                !define new created field to be source
                CALL Field_DependentTypeSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                      & FIELD_INDEPENDENT_TYPE,err,error,*999)
                !look for decomposition rule already defined
                CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(equationsSet%SOURCE%SOURCE_FIELD,equationsSet% &
                  & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                  & 1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionSetAndLock(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentSet(equationsSet%SOURCE%SOURCE_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                  !Specify fem solution method
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%SOURCE%SOURCE_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                    & err,error,*999)
                  CALL Field_ScalingTypeSet(equationsSet%SOURCE%SOURCE_FIELD,geometricScalingType, &
                    & err,error,*999)
                  !Other solutions not defined yet
                CASE DEFAULT
                  localError="The solution method of " &
                    & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                !calculate number of components with one component for each dimension and one for pressure
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                SELECT CASE(equationsSet%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CASE DEFAULT
                  localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                    &"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
              !Specify finish action
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              IF(equationsSet%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                CALL Field_CreateFinish(equationsSet%SOURCE%SOURCE_FIELD,err,error,*999)
                !These 2 parameter sets will contain the fitted hermite/lagrange velocity field
!                 CALL Field_ParameterSetCreate(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
!                 CALL Field_ParameterSetCreate(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)

!                 CALL Field_ParameterSetCreate(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_INPUT_DATA3_SET_TYPE,err,error,*999)
!                 CALL Field_ParameterSetCreate(equationsSet%SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
!                   & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard PEE problem"
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " for a setup sub type of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
              & " is invalid for a PPE equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
! ! !         CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
! ! !           SELECT CASE(equationsSetSetup%ACTION_TYPE)
! ! !           CASE(EQUATIONS_SET_SETUP_START_ACTION)
! ! !             IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               DEPENDENT_FIELD=>equationsSet%DEPENDENT%DEPENDENT_FIELD
! ! !               IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
! ! !                 GEOMETRIC_FIELD=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
! ! !                 IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
! ! !                   CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
! ! !                   SELECT CASE(equationsSetSetup%ANALYTIC_FUNCTION_TYPE)
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_1)
! ! !                     !Check that we are in 2D
! ! !                     IF(numberOfDimensions/=2) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_TWO_DIM_2)
! ! !                     !Check that we are in 2D
! ! !                     IF(numberOfDimensions/=2) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 2 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_TWO_DIM_2
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_1)
! ! !                     !Check that we are in 3D
! ! !                     IF(numberOfDimensions/=4) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_1
! ! !                   CASE(EQUATIONS_SET_FITTING_THREE_DIM_2)
! ! !                     !Check that we are in 3D
! ! !                     IF(numberOfDimensions/=4) THEN
! ! !                       localError="The number of geometric dimensions of "// &
! ! !                         & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
! ! !                         & " is invalid. The analytic function type of "// &
! ! !                         & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                         & " requires that there be 3 geometric dimensions."
! ! !                       CALL FlagError(localError,err,error,*999)
! ! !                     ENDIF
! ! !                     !Create analytic field if required
! ! !                     !Set analtyic function type
! ! !                     equationsSet%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_FITTING_THREE_DIM_2
! ! !                   CASE DEFAULT
! ! !                     localError="The specified analytic function type of "// &
! ! !                       & TRIM(NumberToVString(equationsSetSetup%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
! ! !                       & " is invalid for a standard Galerkin projection."
! ! !                     CALL FlagError(localError,err,error,*999)
! ! !                   END SELECT
! ! !                 ELSE
! ! !                   CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
! ! !                 ENDIF
! ! !              ELSE
! ! !                 CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
! ! !             IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
! ! !               ANALYTIC_FIELD=>equationsSet%ANALYTIC%ANALYTIC_FIELD
! ! !               IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
! ! !                 IF(equationsSet%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
! ! !                   CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
! ! !                 ENDIF
! ! !               ENDIFstandard
! ! !             ELSE
! ! !               CALL FlagError("Equations set analytic is not associated.",err,error,*999)
! ! !             ENDIF
! ! !           CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
! ! !             IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
! ! !               IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
! ! !                 IF(equationsSet%ANALYTIC%ANALYTIC_FINISHED) THEN
! ! !                   CALL FITTING_ANALYTIC_CALCULATE(equationsSet,err,error,*999)
! ! !                 ELSE
! ! !                   CALL FlagError("Equations set analtyic has not been finished.",err,error,*999)
! ! !                 ENDIF
! ! !               ELSE
! ! !                 CALL FlagError("Equations set analtyic is not associated.",err,error,*999)
! ! !               ENDIF
! ! !             ELSE
! ! !               CALL FlagError("Equations set dependent has not been finished.",err,error,*999)
! ! !             ENDIF
! ! !           CASE DEFAULT
! ! !             localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
! ! !               & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
! ! !               & " is invalid for a standard Galerkin projection."
! ! !             CALL FlagError(localError,err,error,*999)
! ! !           END SELECT

        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(equationsSet%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT

      CASE(EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
        & EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
        SELECT CASE(equationsSetSetup%SETUP_TYPE)
        !-----------------------------------------------------------------
        ! s o l u t i o n   m e t h o d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Fitting_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! g e o m e t r y   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          ! Do nothing
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        ! (this field will hold the mesh fitted data from the data points field)
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              !start field creation with name 'DEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                & equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                !define new created field to be dependent
              CALL Field_DependentTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_DEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_MeshDecompositionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,equationsSet% &
                & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
              !set number of variables to 2 (U, delUdelN)
              dependentFieldNumberOfVariables=2
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                & dependentFieldNumberOfVariables,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                & [FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
              CALL Field_VariableLabelSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              DO I=1,numberOfDimensions
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,I, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              END DO
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%DEPENDENT%DEPENDENT_FIELD,geometricScalingType, &
                  & err,error,*999)
                !Other solutions not defined yet
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
            !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        !   m a t e r i a l   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                !start field creation with name 'MATERIAL_FIELD'
                CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION,equationsSet% &
                  & MATERIALS%MATERIALS_FIELD,err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                !apply decomposition rule found on new created field
                CALL Field_MeshDecompositionSetAndLock(equationsSet%MATERIALS%MATERIALS_FIELD, &
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                !point new field to geometric field
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,equationsSet%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL Field_NumberOfVariablesSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                ! Sobolev smoothing material parameters- tau and kappa
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,err,error,*999)
                CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%MATERIALS_FIELD,geometricScalingType,err,error,*999)
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(equationsSetSetup%FIELD,1,err,error,*999)
                CALL Field_VariableTypesCheck(equationsSetSetup%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,err,error,*999)
                CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            END IF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>equationsSet%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL Field_CreateFinish(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for an update-materials Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! I n d e p e n d e n t   t y p e
        ! (this field holds the data point based field of vectors to map to the dependent field)
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
            !Set start action
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created independent field
              !start field creation with name 'INDEPENDENT_FIELD'
              CALL Field_CreateStart(equationsSetSetup%FIELD_USER_NUMBER,equationsSet%REGION, &
                & equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
              !start creation of a new field
              CALL Field_TypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              !label the field
              CALL Field_LabelSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,"Independent Field",err,error, &
                & *999)
              !define new created field to be independent
              CALL Field_DependentTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_INDEPENDENT_TYPE,err,error,*999)
              !look for decomposition rule already defined
              CALL Field_MeshDecompositionGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              !apply decomposition rule found on new created field
              CALL Field_MeshDecompositionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & GEOMETRIC_DECOMPOSITION,err,error,*999)
              !point new field to geometric field
              CALL Field_GeometricFieldSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,equationsSet% &
                & GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
              CALL Field_NumberOfVariablesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & 2,err,error,*999)
              CALL Field_VariableTypesSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & [FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
              ! U Variable: data point vectors
              CALL Field_DimensionSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
              CALL Field_ComponentMeshComponentGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ! V Variable: data point weights
              CALL Field_DimensionSetAndLock(equationsSet%independent%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(equationsSet%independent%INDEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(equationsSet%independent%INDEPENDENT_FIELD, &
                & FIELD_V_VARIABLE_TYPE,1,err,error,*999)
              CALL Field_ComponentMeshComponentGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentSet(equationsSet%independent%INDEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
                !Specify fem solution method
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO dimensionIdx = 1,numberOfDimensions
                  CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,dimensionIdx,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                ENDDO
                CALL Field_ComponentInterpolationSetAndLock(equationsSet%INDEPENDENT%INDEPENDENT_FIELD, &
                  & FIELD_V_VARIABLE_TYPE,1,FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ScalingTypeGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,geometricScalingType, &
                  & err,error,*999)
                CALL Field_ScalingTypeSet(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,geometricScalingType, &
                  & err,error,*999)
                !Other solutions not defined yet
              CASE DEFAULT
                localError="The solution method of " &
                  & //TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(equationsSetSetup%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(equationsSetSetup%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              ! U (vector) variable
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(equationsSet%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              !calculate number of components with one component for each dimension
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              ! V (weight) variable
              CALL Field_DimensionCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              SELECT CASE(equationsSet%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                  & FIELD_DATA_POINT_BASED_INTERPOLATION,err,error,*999)
              CASE DEFAULT
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD, &
                  &"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            !Specify finish action
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(equationsSet%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL Field_CreateFinish(equationsSet%INDEPENDENT%INDEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard PEE problem"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        !   e q u a t i o n s   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(equationsSetSetup%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              IF (equationsSet%specification(3)==EQUATIONS_SET_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE) THEN
                CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
              ELSE IF (equationsSet%specification(3)==EQUATIONS_SET_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE) THEN
                CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(equationsSet%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
              CALL Equations_CreateFinish(equations,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
                & err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// &
                  & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
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
                localError="The solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NumberToVString(equationsSet%specification(3),"*",err,error))// &
          & " does not equal a vector data Galerkin projection subtype."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP")
    RETURN
999 ERRORSEXITS("FITTING_EQUATIONS_SET_VECTORDATA_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FITTING_EQUATIONS_SET_VECTORDATA_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the Fitting problem.
  SUBROUTINE Fitting_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem set to setup a fitting problem on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_ProblemSetup",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
      END IF
      SELECT CASE(problem%specification(3))
      CASE(PROBLEM_STATIC_FITTING_SUBTYPE)
        CALL Fitting_ProblemStaticSetup(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_STANDARD_SETUP(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
        CALL FITTING_PROBLEM_VECTORDATA_SETUP(problem,problemSetup,err,error,*999)
      CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " is not valid for a Galerkin projection type of a data fitting problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_ProblemSetup")
    RETURN
999 ERRORSEXITS("Fitting_ProblemSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Sets up the standard Galerkin projections problem.
  SUBROUTINE FITTING_PROBLEM_STANDARD_SETUP(PROBLEM,problemSetup,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

    ENTERS("FITTING_PROBLEM_STANDARD_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
      END IF
      IF(problem%specification(3)==PROBLEM_STANDARD_DATA_FITTING_SUBTYPE) THEN
        SELECT CASE(problemSetup%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL ControlLoop_CreateStart(PROBLEM,CONTROL_LOOP,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL ControlLoop_CreateFinish(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
          CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL Solvers_CreateStart(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_NumberSet(SOLVERS,1,err,error,*999)
            !Set the solver to be a linear solver
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            CALL Solver_TypeSet(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL Solvers_CreateFinish(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SolverEquations_CreateStart(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SolverEquations_TimeDependenceTypeSet(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SolverEquations_SparsityTypeSet(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            CALL Solver_SolverEquationsGet(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SolverEquations_CreateFinish(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
       CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " does not equal a standard Galerkin projection subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FITTING_PROBLEM_STANDARD_SETUP")
    RETURN
999 ERRORSEXITS("FITTING_PROBLEM_STANDARD_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !

  SUBROUTINE Fitting_ProblemStaticSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop, controlLoopRoot
    TYPE(SOLVER_TYPE), POINTER :: solver
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_ProblemStaticSetup",err,error,*999)

    NULLIFY(controlLoop)
    NULLIFY(solver)
    NULLIFY(solverEquations)
    NULLIFY(solvers)
    IF(ASSOCIATED(problem)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
      END IF
      IF(problem%specification(3)==PROBLEM_STATIC_FITTING_SUBTYPE) THEN
        SELECT CASE(problemSetup%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static fitting problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            controlLoopRoot=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
            CALL ControlLoop_CreateFinish(controlLoop,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static fitting problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          controlLoopRoot=>problem%CONTROL_LOOP
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
            CALL Solvers_NumberSet(solvers,1,err,error,*999)
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
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static fitting problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            controlLoopRoot=>problem%CONTROL_LOOP
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
            controlLoopRoot=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
            !Get the solver equations
            CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
            CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
            CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
            !Finish the solver equations creation
            CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static fitting problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
       CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a static fitting problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " does not equal a static fitting problem."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_ProblemStaticSetup")
    RETURN
999 ERRORSEXITS("Fitting_ProblemStaticSetup",err,error)
    RETURN 1

  END SUBROUTINE Fitting_ProblemStaticSetup

  !
  !================================================================================================================================
  !

  !>Sets up the vector data Galerkin projections problem.
  SUBROUTINE FITTING_PROBLEM_VECTORDATA_SETUP(PROBLEM,problemSetup,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

    ENTERS("FITTING_PROBLEM_VECTORDATA_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
      END IF
      IF(problem%specification(3)==PROBLEM_VECTOR_DATA_FITTING_SUBTYPE.OR. &
        & problem%specification(3)==PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE .OR. &
        & problem%specification(3)==PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE .OR. &
        & problem%specification(3)==PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE) THEN
        SELECT CASE(problemSetup%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL ControlLoop_CreateStart(PROBLEM,CONTROL_LOOP,err,error,*999)
            IF(problem%specification(3)==PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE) THEN
              CALL ControlLoop_TypeSet(CONTROL_LOOP,PROBLEM_CONTROL_SIMPLE_TYPE,err,error,*999)
            ELSE
              CALL ControlLoop_TypeSet(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
            ENDIF
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL ControlLoop_CreateFinish(CONTROL_LOOP,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
          CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL Solvers_CreateStart(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_NumberSet(SOLVERS,1,err,error,*999)
            !Set the solver to be a linear solver
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            CALL Solver_TypeSet(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL Solver_LibraryTypeSet(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL Solvers_CreateFinish(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(problemSetup%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SolverEquations_CreateStart(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SolverEquations_LinearityTypeSet(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            IF(problem%specification(3)==PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE) THEN
              CALL SolverEquations_TimeDependenceTypeSet(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            ELSE
              CALL SolverEquations_TimeDependenceTypeSet(SOLVER_EQUATIONS,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
            ENDIF
            CALL SolverEquations_SparsityTypeSet(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>problem%CONTROL_LOOP
            CALL ControlLoop_Get(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL ControlLoop_SolversGet(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL Solvers_SolverGet(SOLVERS,1,SOLVER,err,error,*999)
            CALL Solver_SolverEquationsGet(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SolverEquations_CreateFinish(SOLVER_EQUATIONS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(problemSetup%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a vector data Galerkin projection."
            CALL FlagError(localError,err,error,*999)
          END SELECT
       CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(problemSetup%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a vector data Galerkin projection."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
          & " does not equal a vector data Galerkin projection subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("FITTING_PROBLEM_VECTORDATA_SETUP")
    RETURN
999 ERRORSEXITS("FITTING_PROBLEM_VECTORDATA_SETUP",err,error)
    RETURN 1
  END SUBROUTINE FITTING_PROBLEM_VECTORDATA_SETUP

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a data fitting problem class.
  SUBROUTINE Fitting_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The proboem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType,problemSubtype

    ENTERS("Fitting_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemType=problemSpecification(2)
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemType)
        CASE(PROBLEM_DATA_FITTING_TYPE)
          SELECT CASE(problemSubtype)
          CASE(PROBLEM_STATIC_FITTING_SUBTYPE, &
            & PROBLEM_STANDARD_DATA_FITTING_SUBTYPE, &
            & PROBLEM_VECTOR_DATA_FITTING_SUBTYPE, &
            & PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE, &
            & PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE, &
            & PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
            !ok
          CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
            CALL FLAG_ERROR("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
              & " is not valid for a Galerkin projection type of a data fitting problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a data fitting problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FITTING_CLASS,problemType,problemSubtype]
      ELSE
        CALL FlagError("Fitting problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated",err,error,*999)
    END IF

    EXITS("Fitting_ProblemSpecificationSet")
    RETURN
999 ERRORS("Fitting_ProblemSpecificationSet",err,error)
    EXITS("Fitting_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE Fitting_ProblemSpecificationSet
  !
  !================================================================================================================================
  !

  !>Evaluates the deformation gradient tensor at a given Gauss point
  SUBROUTINE Fitting_GaussDeformationGradientTensor(referenceGeometricInterpolatedPoint,geometricInterpolatedPoint, &
    & dXdY,Jxy,err,error,*)

    !Argument variables
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: referenceGeometricInterpolatedPoint, geometricInterpolatedPoint
    REAL(DP) :: dXdY(3,3)  !dXdY - Deformation Gradient Tensor
    REAL(DP) :: Jxy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: derivativeIdx,componentIdx,xiIdx
    REAL(DP) :: dXdXi(3,3),dYdXi(3,3),dXidY(3,3)
    REAL(DP) :: Jyxi

    ENTERS("Fitting_GaussDeformationGradientTensor",err,error,*999)

    !--- ToDo: Needs to be generalized such that it also works for 2D
    DO componentIdx=1,3 !Always 3 components - 3D
      DO xiIdx=1,3 !Thus 3 element coordinates
        derivativeIdx=PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx) !2,4,7
        dXdXi(componentIdx,xiIdx)=geometricInterpolatedPoint%values(componentIdx,derivativeIdx) !dx/dxi
        dYdXi(componentIdx,xiIdx)=referenceGeometricInterpolatedPoint%values(componentIdx,derivativeIdx) !dy/dxi (y = referential)
      ENDDO
    ENDDO

    CALL Invert(dYdXi,dXidY,Jyxi,err,error,*999) !dy/dxi -> dxi/dy
    CALL MatrixProduct(dXdXi,dXidY,dXdY,err,error,*999) !dx/dxi * dxi/dy = dx/dy (deformation gradient tensor, F)
    CALL Determinant(dXdY,Jxy,err,error,*999)

    EXITS("Fitting_GaussDeformationGradientTensor")
    RETURN
999 ERRORSEXITS("Fitting_GaussDeformationGradientTensor",err,error)
    RETURN 1

  END SUBROUTINE Fitting_GaussDeformationGradientTensor

  !
  !================================================================================================================================
  !

 !>Sets up the output type for a data fitting problem class.
  SUBROUTINE Fitting_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_PreSolve",err,error,*999)

    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%solvers
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop)) THEN
          problem=>controlLoop%problem
          IF(ASSOCIATED(problem)) THEN
            IF(.NOT.ALLOCATED(problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE
              IF(SIZE(problem%specification,1)<3) THEN
                CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
              END IF
            ENDIF
            SELECT CASE(problem%specification(3))
            CASE(PROBLEM_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
              !IF(controlLoop%WHILE_LOOP%ITERATION_NUMBER==1)THEN
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
              !Update indpendent data fields
              CALL Fitting_PreSolveUpdateInputData(solver,err,error,*999)
              !  CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
              !ELSE
              !  CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
              !ENDIF
            CASE DEFAULT
              localError="The third problem specification of "// &
                & TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a data fitting problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Control loop problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solvers control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_PreSolve")
    RETURN
999 ERRORSEXITS("Fitting_PreSolve",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PreSolve

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a data fitting problem class.
  SUBROUTINE Fitting_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("Fitting_PostSolve",err,error,*999)

    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%solvers
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop)) THEN
          problem=>controlLoop%problem
          IF(ASSOCIATED(problem)) THEN
            IF(.NOT.ALLOCATED(problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE
              IF(SIZE(problem%specification,1)<3) THEN
                CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
              ENDIF
            ENDIF
            SELECT CASE(problem%specification(3))
            CASE(PROBLEM_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE,PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE, &
              & PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE, &
              & PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE)
              CALL Fitting_PostSolveOutputData(solver,err,error,*999)
            CASE(PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_PRE_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
              !Do nothing
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a fitting problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Control loop problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solvers control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_PostSolve")
    RETURN
999 ERRORSEXITS("Fitting_PostSolve",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PostSolve

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE Fitting_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(PROBLEM_TYPE), POINTER :: problem  !<A pointer to the solver equations
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError
    REAL(DP) :: currentTime,timeIncrement
    INTEGER(INTG) :: equationsSetIdx,currentLoopIteration,outputIterationNumber
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlTimeLoop !<A pointer to the control loop to solve.
    LOGICAL :: exportField
    CHARACTER(7) :: outputFile

    ENTERS("Fitting_PostSolveOutputData",err,error,*999)

    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%solvers
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop)) THEN
          problem=>controlLoop%problem
          IF(ASSOCIATED(problem)) THEN
            IF(.NOT.ALLOCATED(problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE
              IF(SIZE(problem%specification,1)<3) THEN
                CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
              ENDIF
            ENDIF
            SELECT CASE(problem%specification(3))
            CASE(PROBLEM_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE,PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE, &
              & PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE, &
              & PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_PRE_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_PRE_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
              controlTimeLoop=>controlLoop
              CALL ControlLoop_CurrentTimesGet(controlTimeLoop,currentTime,timeIncrement,err,error,*999)
              solverEquations=>solver%SOLVER_EQUATIONS
              IF(ASSOCIATED(solverEquations)) THEN
                solverMapping=>solverEquations%SOLVER_MAPPING
                IF(ASSOCIATED(solverMapping)) THEN
                  !Make sure the equations sets are up to date
                  DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
                    equationsSet=>solverMapping%EQUATIONS_SETS(equationsSetIdx)%ptr
                    currentLoopIteration=controlTimeLoop%TIME_LOOP%ITERATION_NUMBER
                    outputIterationNumber=controlTimeLoop%TIME_LOOP%OUTPUT_NUMBER
                    IF(outputIterationNumber/=0) THEN
                      IF(controlTimeLoop%TIME_LOOP%CURRENT_TIME<=controlTimeLoop%TIME_LOOP%STOP_TIME) THEN
                        IF(currentLoopIteration<10) THEN
                          WRITE(outputFile,'("DATA_0",I0)') currentLoopIteration
                        ELSE IF(currentLoopIteration<100) THEN
                          WRITE(outputFile,'("DATA_",I0)') currentLoopIteration
                        ENDIF
                        exportField=.TRUE.
                        IF(exportField) THEN
                          IF(MOD(currentLoopIteration,outputIterationNumber)==0)  THEN
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                            CALL FLUID_MECHANICS_IO_WRITE_FITTED_FIELD(equationsSet%region,equationsSet%GLOBAL_NUMBER, &
                              & outputFile,err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                            CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          ENDIF
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a fitting equation of a classical field problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Control loop problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solvers control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Fitting_PostSolveOutputData",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Update input data conditions for field fitting
  SUBROUTINE Fitting_PreSolveUpdateInputData(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop,controlTimeLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(SOLVERS_TYPE), POINTER :: solvers
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: numberOfDimensions,currentLoopIteration
    INTEGER(INTG) :: inputType,inputOption
    REAL(DP), POINTER :: inputVelNewData(:)
    REAL(DP) :: currentTime,timeIncrement
    LOGICAL :: boundaryUpdate

    boundaryUpdate=.FALSE.

    ENTERS("Fitting_PreSolveUpdateInputData",err,error,*999)

    NULLIFY(inputVelNewData)

    IF(ASSOCIATED(solver)) THEN
      solvers=>solver%solvers
      IF(ASSOCIATED(solvers)) THEN
        controlLoop=>solvers%CONTROL_LOOP
        IF(ASSOCIATED(controlLoop)) THEN
          problem=>controlLoop%problem
          IF(ASSOCIATED(problem)) THEN
            IF(.NOT.ALLOCATED(problem%specification)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE
              IF(SIZE(problem%specification,1)<3) THEN
                CALL FlagError("Problem specification must have three entries for a fitting problem.",err,error,*999)
              ENDIF
            ENDIF
            SELECT CASE(problem%specification(3))
            CASE(PROBLEM_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_STANDARD_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_GENERALISED_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_MAT_PROPERTIES_DATA_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_DATA_POINT_VECTOR_STATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_DATA_POINT_VECTOR_QUASISTATIC_FITTING_SUBTYPE)
              !Do nothing
            CASE(PROBLEM_VECTOR_DATA_FITTING_SUBTYPE,PROBLEM_DIV_FREE_VECTOR_DATA_FITTING_SUBTYPE)
              !Do nothing
              controlTimeLoop=>controlLoop
              CALL ControlLoop_CurrentTimesGet(controlTimeLoop,currentTime,TimeIncrement,err,error,*999)
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
              solverEquations=>solver%SOLVER_EQUATIONS
              IF(ASSOCIATED(solverEquations)) THEN
                solverMapping=>solverEquations%SOLVER_MAPPING
                equations=>solverMapping%EQUATIONS_SET_TO_SOLVER_MAP(1)%equations
                IF(ASSOCIATED(equations)) THEN
                  equationsSet=>equations%equationsSet
                  IF(ASSOCIATED(equationsSet)) THEN
                    CALL Field_NumberOfComponentsGet(equationsSet%geometry%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & numberOfDimensions,err,error,*999)
                    currentLoopIteration=controlTimeLoop%TIME_LOOP%ITERATION_NUMBER
                    !this is the current time step
                    !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
                    inputType=1
                    inputOption=1
                    CALL Field_ParameterSetDataGet(equationsSet%source%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,inputVelNewData,err,error,*999)
                    CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelNewData, &
                      & numberOfDimensions,inputType,inputOption,currentLoopIteration,1.0_DP,err,error,*999)
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("Solver equations are not associated.",err,error,*999)
              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
                & " is not valid for a vector data type of a fitting field problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL Field_ParameterSetUpdateStart(equationsSet%source%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(equationsSet%source%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VALUES_SET_TYPE,err,error,*999)
          ELSE
            CALL FlagError("Control loop problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solvers control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("Fitting_PreSolveUpdateInputData")
    RETURN
999 ERRORSEXITS("Fitting_PreSolveUpdateInputData",err,error)
    RETURN 1

  END SUBROUTINE Fitting_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !


END MODULE FittingRoutines
