!> \file
!> \author Chris Bradley
!> \brief This module handles all Poisson equations routines.
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

!>This module handles all Poisson equations routines.
MODULE PoissonEquationsRoutines

  USE BaseRoutines 
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines  
  USE DistributedMatrixVector
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
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
! temporary input for vector-source
  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Poisson_BoundaryConditionsAnalyticCalculate
  
  PUBLIC Poisson_EquationsSetSetup

  PUBLIC Poisson_EquationsSetSolutionMethodSet

  PUBLIC Poisson_EquationsSetSpecificationSet

  PUBLIC Poisson_FiniteElementCalculate

  PUBLIC Poisson_FiniteElementJacobianEvaluate,Poisson_FiniteElementResidualEvaluate

  PUBLIC Poisson_PreSolve,Poisson_PostSolve
    
  PUBLIC Poisson_ProblemSpecificationSet

  PUBLIC Poisson_ProblemSetup

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Poisson_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,derivativeGlobalIndex,dimensionIdx,localDOFIdx,nodeIdx, &
      & numberOfComponents,numberOfDimensions,numberOfNodeDerivatives,numberOfNodes,variableIdx,variableType
    REAL(DP) :: dependentValue,x(3)
    REAL(DP), POINTER :: geometricParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
   TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: dependentField,geometricField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(VARYING_STRING) :: localError    
   
    ENTERS("Poisson_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
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
    DO variableIdx=1,numberOfVariables !U and deludeln
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents !u,v,w
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
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
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          !Loop over the derivatives
          DO derivativeIdx=1,numberOfNodesDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,derivativeGlobalIndex,err,error,*999)
            SELECT CASE(analyticFunctionType)
            CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
              !u=ln(4/(x+y+1^2))
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=LOG(4.0_DP/((x(1)+x(2)+1.0_DP)**2))
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  !This is believed to be incorrect, should be: dependentValue=-2.0_DP/(x(1)+x(2)+1.0_DP)
                  dependentValue=-2.0_DP*(x(1)+x(2))/(x(1)+x(2)+1.0_DP)
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              END SELECT
            CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_2)
              CALL FlagError("The analytic function type is not implemented yet.",err,error,*999)
            CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_3)
              CALL FlagError("The analytic function type is not implemented yet.",err,error,*999)
            CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
              !u=ln(6/(x+y+z+1^2))
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=LOG(6.0_DP/((x(1)+x(2)+x(3)+1.0_DP)**2))
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=-3.0_DP/(x(1)+x(2)+x(3)+1.0_DP)
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE DEFAULT
                localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_2)
              CALL FlagError("The analytic function type is not implemented yet.",err,error,*999)
            CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_3)
              CALL FlagError("The analytic function type is not implemented yet.",err,error,*999)
            CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1,EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2)
              !\todo: This test case has been set up for the Pressure Poisson equation -
              !needs to be formulated in a more general way
              !u=ln(4/(x+y+1^2))
              SELECT CASE(variableType)
              CASE(FIELD_U_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  dependentValue=SIN(2.0_DP*PI*x(1)/10.0_DP)*SIN(2.0_DP*PI*x(2)/10.0_DP)*SIN(2.0_DP*PI*x(3)/10.0_DP)
                  !  dependentValue=2*x(1)*x(1)+2*x(2)
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                SELECT CASE(derivativeGlobalIndex)
                CASE(NO_GLOBAL_DERIV)
                  !do nothing, tbd
                CASE(GLOBAL_DERIV_S1)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE(GLOBAL_DERIV_S1_S2)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  localError="The global derivative index of "//TRIM(NumberToVString(derivativeGlobalIndex,"*",err,error))// &
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
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
              & dependentValue,err,error,*999)
            CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            IF(variableType==FIELD_U_VARIABLE_TYPE.AND.boundaryNode) THEN !For Dirichlet
              CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentField,variableType, &
                & localDOFIdx,BOUNDARY_CONDITION_FIXED,dependentValue,err,error,*999)
            ENDIF
            IF(variableType==FIELD_DELUDELN_VARIABLE_TYPE.AND.nodeIdx/=1) THEN
              IF(boundaryNode) THEN
                !If we are a boundary node then set the analytic value on the boundary
                !CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentField,variableType,localDOFIdx, &
                !  & BOUNDARY_CONDITION_FIXED,dependentValue,err,error,*999)
                !Do nothing at present
              ENDIF
            ENDIF

          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ELSE
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)

    ENDDO !variableIdx

    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
   
    EXITS("Poisson_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Poisson_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Sets up the Poisson equation.
  SUBROUTINE Poisson_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricComponentNumber,geometricScalingType,numberOfDimensions, &
      & numberOfMaterialsComponents,geometricMeshComponent,numberOfSourceComponents,numberOfSourceVariables,solutionMethod
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Poisson_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE,EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid."
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
        CALL Poisson_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Poisson equation."
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
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE) THEN
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
              & err,error,*999)    
          ELSE IF(esSpecification(3)==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"P",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del P/del n", &
              & err,error,*999)    
          ELSE
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
              & err,error,*999)    
          ENDIF
          !Default to the geometric interpolation setup
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricComponentNumber,err,error,*999)
          CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
            & geometricComponentNumber,err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
          IF(esSpecification(3)==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
            !These 3 parameter sets will contain the original velocity field with x,y,z components
            CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_VEL1_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_VEL2_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_VEL3_SET_TYPE,err,error,*999)
            !This parameter set will contain the inside/outside labelling information for the PPE approach
            CALL Field_ParameterSetCreate(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_LABEL_SET_TYPE,err,error,*999)
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Poisson equation"
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
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_LabelSet(equationsSet%MATERIALS%materialsField,"Material Field",err,error,*999)                   
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,equationsSet%geometry%geometricField, &
            & err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
            & err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
            !Constant source. Materials field components are for the conductivity tensor \sigma
            !i.e., div(\sigma(x).grad(u(x)))+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err, &
              & error,*999)
          CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
            !Linear source. Materials field components are for the conductivity tensor \sigma plus a
            !i.e., div(\sigma(x).grad(u(x)))+a(x)u(x)+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+1
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err, &
              & error,*999)
          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
            !Quadratic source. Materials field components are for the conductivity tensor \sigma plus a and b
            !i.e., div(\sigma(x).grad(u(x)))+b(x)u(x)^2+a(x)u(x)+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+2
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err, &
              & error,*999)
          CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            !Exponential source. Matierals field components are for the conductivity tensor \sigma plus a and b
            !i.e., div(\sigma(x).grad(u(x)))+a(x)e^[b(x)u(x)]+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+2
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err, &
              & error,*999)
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err, &
              & error,*999)
            !\mu and \rho for the source
            numberOfMaterialsComponents=2
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            !Extracellular bidomain. Material field components are the intra and extra-cellular conductivity tensors
            !i.e.. div([\sigma_i(x)+\simga_e(x)].grad(u(x))+div(\simga_i(x).grad(V_m(x)))=0
            CALL Field_VariableLabelSet(equationsSet%materials%materialsField,FIELD_U_VARIABLE_TYPE,"Conductivity",err, &
              & error,*999)
            numberOfMaterialsComponents=2*NUMBER_OF_VOIGT(numberOfDimensions)
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
          !Set the number of materials components
          CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          !Default the  materials components to the first component geometric interpolation with constant interpolation
          CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
            & 1,geometricComponentNumber,err,error,*999)
          DO componentIdx=1,numberOfMaterialsComponents
            CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricComponentNumber,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
          ENDDO !components_idx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDimensions,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
            !Constant source. Materials field components are for the conductivity tensor \sigma
            !i.e., div(\sigma(x).grad(u(x)))+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)
          CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
            !Linear source. Materials field components are for the conductivity tensor \sigma plus a
            !i.e., div(\sigma(x).grad(u(x)))+a(x)u(x)+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+1
          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
            !Quadratic source. Materials field components are for the conductivity tensor \sigma plus a and b
            !i.e., div(\sigma(x).grad(u(x)))+b(x)u(x)^2+a(x)u(x)+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+2
          CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            !Exponential source. Matierals field components are for the conductivity tensor \sigma plus a and b
            !i.e., div(\sigma(x).grad(u(x)))+a(x)e^[b(x)u(x)]+s(x)=0
            numberOfMaterialsComponents=NUMBER_OF_VOIGT(numberOfDimensions)+2
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            !\rho and \mu for the source
            numberOfMaterialsComponents=2
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            !Extracellular bidomain. Material field components are the intra and extra-cellular conductivity tensors
            !i.e.. div([\sigma_i(x)+\simga_e(x)].grad(u(x))+div(\simga_i(x).grad(V_m(x)))=0
            numberOfMaterialsComponents=2*NUMBER_OF_VOIGT(numberOfDimensions)
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
         CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
            & err,error,*999)
        ENDIF        
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values for the materials field
          CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDimensions,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
            & numberOfMaterialsComponents,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            !Default the intra and extra cellular conductivity tensors to a diagonal identity tensor
            DO componentIdx=1,numberOfDimensions
              componentNumber=TENSOR_TO_VOIGT(componentIdx,componentIdx,numberOfDimensions)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentNumber,1.0_DP,err,error,*999)
            ENDDO !componentIdx
            !Default the remainder of the constants to 1.0
            DO componentIdx=NUMBER_OF_VOIGT(numberOfDimensions)+1,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            !\mu and \rho for the source
            !Default the values to 1.0
            DO componentIdx=1,numberOfMaterialsComponents
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx            
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            !Default the intra and extra cellular conductivity tensors to a diagonal identity tensor
            DO componentIdx=1,numberOfDimensions
              componentNumber=TENSOR_TO_VOIGT(componentIdx,componentIdx,numberOfDimensions)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentNumber,1.0_DP,err,error,*999)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentNumber+NUMBER_OF_VOIGT(numberOfDimensions),1.0_DP,err,error,*999)
            ENDDO !componentIdx
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a velocity source Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d 
      !-----------------------------------------------------------------
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        !Set start action
        IF(equationsSet%source%sourceFieldAutoCreated) THEN
          !Create the auto created source field
          !Start field creation with name 'Source Field'
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%region,equationsSet%source%sourceField, &
            & err,error,*999)
          !start creation of a new field
          CALL Field_TypeSetAndLock(equationsSet%source%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
          !define new created field to be source
          CALL Field_DependentTypeSetAndLock(equationsSet%source%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          !look for decomposition rule already defined
          CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
          !apply decomposition rule found on new created field
          CALL Field_DecompositionSetAndLock(equationsSet%source%sourceField,geometricDecomposition,err,error,*999)
          !point new field to geometric field
          CALL Field_GeometricFieldSetAndLock(equationsSet%source%sourceField,equationsSet%geometry%geometricField, &
            & err,error,*999)
          !set number of variables to 1 (1 for U)
          numberOfSourceVariables=1
          CALL Field_NumberOfVariablesSetAndLock(equationsSet%source%sourceField,numberOfSourceVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSet%source%sourceField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
            & err,error,*999)
          !calculate number of components with one component for each dimension
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            !label the field
            CALL Field_LabelSet(equationsSet%source%sourceField,"Source Field",err,error,*999)
            numberOfSourceComponents=1
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            !label the field
            CALL Field_LabelSet(equationsSet%source%sourceField,"Source Field",err,error,*999)
            numberOfSourceComponents=numberOfDimensions
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            !label the field
            CALL Field_LabelSet(equationsSet%source%sourceField,"Vm",err,error,*999)
            numberOfSourceComponents=1
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE, &
            & numberOfSourceComponents,err,error,*999)
          CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,1, &
            & geometricMeshComponent,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfSourceComponents
            CALL Field_ComponentMeshComponentSet(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          !Specify fem solution method
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfSourceComponents
              CALL Field_ComponentInterpolationSetAndLock(equationsSet%source%sourceField, &
                & FIELD_U_VARIABLE_TYPE,componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsSet%source%sourceField,geometricScalingType,err,error,*999)
            !Other solutions not defined yet
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDimensions,err,error,*999)
          !calculate number of components with one component for each dimension and one for pressure
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &          
            & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            numberOfSourceComponents=1
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            numberOfSourceComponents=numberOfDimensions
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            numberOfSourceComponents=1
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfSourceComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          CASE DEFAULT
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ENDIF
        !Specify finish action
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%source%sourceFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%source%sourceField,err,error,*999)
          CALL Field_NumberOfComponentsGet(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE,numberOfSourceComponents, &
            & err,error,*999)
          !Initialise field
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &          
            & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
            !Default source field to 1.0
            DO componentIdx=1,numberOfSourceComponents
              CALL Field_ComponentValuesInitialise(equationsSet%source%surceField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
          CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
            & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
            !Do nothing, default source field to 0.0
            !Create parameter sets
            !These 2 parameter sets will contain the fitted hermite/lagrange velocity field
            CALL Field_ParameterSetCreate(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
            !CALL Field_ParameterSetCreate(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE, &
            !  & FIELD_INPUT_DATA3_SET_TYPE,err,error,*999)
            !CALL Field_ParameterSetCreate(equationsSet%source%sourceField,FIELD_U_VARIABLE_TYPE, &
            !  & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            !Do nothing, default Vm field to 0.0
          CASE(DEFAULT)
            localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)            
          END SELECT
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! I n d e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      !define an independent field for ALE information
      IF(esSpecification(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
        & esSpecification(3)==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
        NULLIFY(equationsIndependent)
        CALL EquationsSet_IndependentGet(equationsSet,equationsIndependent,err,error,*999)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            !Create the auto created independent field
            !start field creation with name 'independentField'
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,equationsSet%region, &
              & equationsIndependent%independentField,err,error,*999)
            !start creation of a new field
            CALL Field_TypeSetAndLock(equationsIndependent%independentField,FIELD_GENERAL_TYPE,err,error,*999)
            !label the field
            CALL Field_LabelSetAndLock(equationsIndependent%independentField,"Independent Field",err,error,*999)
            !define new created field to be independent
            CALL Field_DependentTypeSetAndLock(equationsIndependent%independentField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            !look for decomposition rule already defined
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(equationsSet%geometry%geometricField,geometricDecomposition,err,error,*999)
            !apply decomposition rule found on new created field
            CALL Field_DecompositionSetAndLock(equationsIndependent%independentField,geometricDecomposition,err,error,*999)
            !point new field to geometric field
            CALL Field_GeometricFieldSetAndLock(equationsIndependent%independentField,equationsSet%geometry%geometricField, &
              & err,error,*999)
            !set number of variables to 1 (1 for U)
            CALL Field_NumberOfVariablesSetAndLock(equationsIndependent%independentField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsIndependent%independentField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
              & err,error,*999)
            !calculate number of components with one component for each dimension
            CALL Field_NumberOfComponentsSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfDimensions,err,error,*999)
            CALL Field_ComponentMeshComponentGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,1, &
              & geometricMeshComponent,err,error,*999)
            !Default to the geometric interpolation setup
            DO componentIdx=1,numberOfDimensions
              CALL Field_ComponentMeshComponentSet(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Specify fem solution method
              DO componentIdx=1,numberOfDimensions
                CALL Field_ComponentInterpolationSetAndLock(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
                  & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              ENDDO !componentIdx
              CALL Field_ScalingTypeGet(equationsSet%geometry%geometricField,geometricScalingType,err,error,*999)
              CALL Field_ScalingTypeSet(equationsIndependent%independentField,geometricScalingType,& err,error,*999)
            CASE DEFAULT
              !Other solutions not defined yet
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsGet(equationsSet%geometry%geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
              & err,error,*999)
            !calculate number of components with one component for each dimension and one for pressure
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
            CASE DEFAULT
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
          !Specify finish action
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(equationsIndependent%independentFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsIndependent%independentField,err,error,*999)
            CALL Field_ParameterSetCreate(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_MESH_VELOCITY_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetCreate(equationsIndependent%independentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_BOUNDARY_SET_TYPE,err,error,*999)
          ENDIF
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a standard PPE fluid"
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
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
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
          IF(numberOfDimensions==2) THEN
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
              CALL EquationsSet_AnalyticFunctionTypeSet(equationsSet,EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "// &
                & /TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a two-dimensional constant source Poisson equation."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE IF(numberOfDimensions==3) THEN
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
              CALL EquationsSet_AnalyticFunctionTypeSet(equationsSet,EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "// &
                & /TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a three-dimensional constant source Poisson equation."
              CALL FlagError(localError,err,error,*999)
            ENDIF
         ELSE
            localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid for a constant source Poisson equation."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
          IF(numberOfDimensions==2) THEN
!!TODO: Why do we have the same analytic function type here as for the constant source???
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_POISSON_EQUATION_TWO_DIM_1)
              CALL EquationsSet_AnalyticFunctionTypeSet(equationsSet,EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "// &
                & /TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a two-dimensional exponential source Poisson equation."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE IF(numberOfDimensions==3) THEN
!!TODO: Why do we have the same analytic function type here as for the constant source???
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1)
              CALL EquationsSet_AnalyticFunctionTypeSet(equationsSet,EQUATIONS_SET_POISSON_EQUATION_THREE_DIM_1,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "// &
                & /TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a three-dimensional exponential source Poisson equation."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid for an exponential source Poisson equation."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
          & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
          IF(numberOfDimensions==3) THEN
            SELECT CASE(equationsSetSetup%analyticFunctionType)
            CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1)
              CALL EquationsSet_AnalyticFunctionTypeSet(equationsSet,EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1,err,error,*999)
            CASE DEFAULT
              localError="The analytic function type of "// &
                & /TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " is invalid for a three-dimensional pressure Poisson equation."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ELSE
            localError="The number of dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
              & " is invalid for a pressure Poisson equation."
            CALL FlagError(localError,err,error,*999)
          ENDIF
        CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(DEFAULT)
          localError="The third equations set specification of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)            
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(.NOT.ASSOCIATED(equationsSet%analytic)) CALL FlagError("Equations set analytic is not associated.",err,error,*999)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsSet%analytic%analyticFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
          ENDIF          
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a velocity source Poisson equation."
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
        IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        ELSE
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
        ENDIF
        IF(esSpecification(3)==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE) THEN
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_QUASISTATIC,err,error,*999)
        ELSE
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)          
        SELECT CASE(solutionMethod)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          !Finish the creation of the equations
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          CALL Equations_CreateFinish(equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          !Create the equations mapping.
          CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
          CALL EquationsMapping_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
          CALL EquationsMappingVector_SourceVariableTypeSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
          SELECT CASE(sparsityType)
          CASE(EQUATIONS_MATRICES_FULL_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
            ENDIF
          CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
            CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE.OR. &
            & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE) THEN
            !Use the analytic Jacobian calculation
            CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
              & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
          ENDIF
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
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a Poisson equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Poisson_EquationsSetSetup")
    RETURN
999 ERRORS("Poisson_EquationsSetSetup",err,error)
    EXITS("Poisson_EquationsSetSetup")
    RETURN 1
    
  END SUBROUTINE Poisson_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Poisson equation type of an classical field equations set class.
  SUBROUTINE Poisson_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Poisson_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
   
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
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
      localError="Equations set subtype of "//TRIM(NumberToVString(equationsSet%SPECIFICATION(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Poisson_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Poisson_EquationsSetSolutionMethodSet",err,error)
    EXITS("Poisson_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Poisson_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Poisson equation type of a classical field equations set class.
  SUBROUTINE Poisson_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poisson_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    subtype=specification(3)
    
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE, &
      & EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
        & " is not valid for a Poisson equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_POISSON_EQUATION_TYPE,subtype]

    EXITS("Poisson_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Poisson_EquationsSetSpecificationSet",err,error)
    EXITS("Poisson_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Poisson_EquationsSetSpecificationSet
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE Poisson_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnXiIdx,columnElementParameterIdx,dofIdx,gaussPointIdx, &
      & nodeIdx,numberOfElementParameters,numberOfDimensions,numberOfXi,rowComponentIdx,rowElementDOFIdx,rowBasisIdx, &
      & rowElementParameterIdx
    REAL(DP) :: aParam,b(3),colsdPhidXi(3),colsPhi,conductivity(3,3),deltaT,diffCoeff1,diffCoeff2,dXidX(3,3),d2XidX2(3,3), &
      & extraConductivity(3,3),intraConductivity(3,3),jacobianGaussWeight,muParam,pDeriv(3),rhoParam,rowsdPhidXi(3),rowsPhi, &
      & sourceValue,sum,sum2,uDeriv(3,3),uOld(3),uSecond(3,3,3),uValue(3),Vm(64),wValue(3),x(3)
    LOGICAL :: between,inside,update,updateMatrix,updateRHS,updateSource
    TYPE(BasisType), POINTER :: columnBasis,geometricBasis,independentBasis,rowBasis,sourceBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,sourceField,materialsField,independentField,fibreField
    TYPE(FieldVariableType), POINTER :: columnVariable,fieldVariable,rowVariable,sourceVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme
    TYPE(FieldInterpolatedPointType), POINTER :: geometricInterpolatedPoint,fibreInterpolatedPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Poisson_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    !Get the fields and check the equations
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
    NULLIFY(fibreField)
    NULLIFY(independentField)
    NULLIFY(sourceField)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
      !div(\sigma.grad u) + s = 0
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.u + s = 0
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + b.u^2 + a.u + s = 0
      CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.e^(b.u) + s = 0
      CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
    CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
      !div(grad p) -  div( \mu div(grad u) - \rho (\del u/\del t - u . grad u ) ) = 0
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE) THEN
        CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      ENDIF
    CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      !div((\sigma_i+\sigma_e).grad \phi_e) + div(\sigma_i. grad V_m) = 0
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      Vm=0.0_DP
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*9999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(rowsVariable)
    CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
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
    updateMatrix=equationsMatrix%updateMatrix
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      updateSource=sourceVector%updateVector
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      updateRHS=rhsVector%updateVector
    ENDIF

    update=(updateMatrix.OR.updateSource.OR.updateRHS)

    IF(update) THEN
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      
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
      
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(colsDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(colsDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
      NULLIFY(colsBasis)
      CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(colsQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(colsQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(materialsInterpParamters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
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
      
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      IF(ASSOCIATED(sourceField)) THEN
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
      ENDIF
    
      NULLIFY(independentInterpParameters)
      NULLIFY(independentInterpPoint)
      IF(ASSOCIATED(independentField)) THEN
        CALL EquationsInterpolation_IndependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpParameters,err,error,*999)
        CALL EquationsInterpolation_IndependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & independentInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,independentInterpParameters, &
          & err,error,*999)
      ENDIF
      
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
        !Do nothing
      CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
        !Do nothing
      CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
        & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
        CALL EquationsSet_TimeIncrementGet(equationsSet,deltaT,err,error,*999)
        diffCoeff1=1.0_DP
        diffCoeff2=1.0_DP
        !Determine inside outside nodes/elements
        CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_LABEL_SET_TYPE,inputLabel,err,error,*999)
        CALL Basis_NumberOfLocalNodesGet(colsBasis,numberOfLocalNodes,err,error,*999)
        DO localNodeIdx=1,numberOfLocalNodes
          CALL DomainElements_ElementNodeGet(colsDomainElements,localNodeIdx,elementNumber,elementNode,err,error,*999)
          !this needs to be changed even at the right inputLabel
          IF(inputLabel(elementNode)<0.5_DP) THEN
            !labelling if !!!ELEMENT!!! is inside or outside
            !so even if only one node is outside node, the whole element is outside element.
            inside=.FALSE.
          ENDIF
        ENDDO !localNodeIdx
        CALL Basis_NumberOfLocalNodesGet(geometricBasis,numberOfLocalNodes,err,error,*999)
        DO localNodeIdx=1,numberOfLocalNodes
          CALL DomainElements_ElementNodeGet(geometricDomainElements,localNodeIdx,elementNumber,elementNode,err,error,*999)
          !this needs to be changed even at the right inputLabel
          IF(.NOT.inside) THEN
            IF(inputLabel(elementNode)>0.5_DP) THEN
              !labelling if !!!ELEMENT!!! is inside or outside
              !so even if only one node is outside node, the whole element is outside element.
              between=.TRUE.
            ENDIF
          ENDIF
        ENDDO !localNodeIdx
        
        IF(inside) THEN
          diffCoeff1=1.0_DP
          diffCoeff2=1.0_DP
        ELSE IF(between) THEN
          diffCoeff1=1.0_DP/1000.0_DP
          diffCoeff2=0.0_DP
        ELSE
          diffCoeff1=0.0_DP
          diffCoeff2=0.0_DP
        ENDIF
        NULLIFY(sourceVariable)
        CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
        CALL FieldVariable_InterpolationParametersInitialise(sourceVariable,oldSourceInterpParameters,err,error,*999)
        NULLIFY(oldSourceInterpPoint)
        CALL Field_InterpolatedPointInitialise(oldSourceInterpParameters,oldSourceInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_INPUT_DATA2_SET_TYPE,elementNumber,oldSourceInterpParameters, &
          & err,error,*999)
        IF(ASSOCIATED(equationsSet%analytic)) &
          & CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
      CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
        Vm=0.0_DP
        !get correct Vm-vector entries (rows columnElementDOFIdx)      
        CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
        DO columnElementDOFIdx=1,numberOfColsElementParameters
          !get correct global node index for this element
          nodeIdx=sourceVector%elementVector%rowDOFS(columnElementDOFIdx)
          !get corresponding dofIdx (in this case nodeIdx=dofIdx as just 1 component, but leave for generality)
          CALL FieldVariable_LocalNodeDOFGet(sourceVariable,1,1,nodeIdx,1,dofIdx,err,error,*999)
          !get the correct value from global source field
          CALL FieldVariable_ParameterSetGetLocalDOF(sourceVariable,FIELD_VALUES_SET_TYPE, &
            & dofIdx,Vm(columnElementDOFIdx),err,error,*999)
        ENDDO !columnElementDOFIdx
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a Poisson equation type of a classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      
      !Loop over the Gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        !Interpolate fields
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
          IF(esSpecification(3)==EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE) &
            & aParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+1,NO_PART_DERIV)
          IF(updateSource) THEN
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
              & err,error,*999)
            sourceValue=sourceInterpPoint%values(1,NO_PART_DERIV)
          ENDIF
          !Calculate conductivity tensor
          IF(ASSOCIATED(fibreField)) &
            & CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)  
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & materialsInterpPoint%values(1:NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)        
        CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
          wValue=0.0_DP
          IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE) THEN
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,independentInterpPoint, &
              & err,error,*999)
            wValue(1:numberOfDimensions)=independentInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
          ENDIF
          uValue=0.0_DP
          uOld=0.0_DP
          dXidX=0.0_DP
          d2XidX2=0.0_DP       
          uDeriv=0.0_DP
          uSecond=0.0_DP
          IF(updateSource) THEN
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
              & err,error,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,oldSourceInterpPoint, &
              & err,error,*999)
            DO componentIdx=1,numberOfDimensions
              uValue(componentIdx)=sourceInterpPoint%values(componentIdx,NO_PART_DERIV)
              uOld(componentIdx)=oldSourceInterpPoint%values(componentIdx,NO_PART_DERIV)
              DO xiIdx1=1,numberOfXi
                uDeriv(componentIdx,xiIdx1)=sourceInterpPoint%values(componentIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx1))
                DO xiIdx2=1,numberOfXi
                  uSecond(componentIdx,xiIdx1,xiIdx2)=sourceInterpPoint%values(componentIdx, &
                    & PARTIAL_DERIVATIVE_SECOND_DERIVATIVES_MAP(xiIdx1,xiIdx2))              
                ENDDO !xiIdx2            
              ENDDO !xiIdx1
            ENDDO !componentIdx
            DO xiIdx1=1,numberOfXi
              pDeriv(xiIdx1)=dependentInterpPoint%values(1,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx1))
            ENDDO !xiIdx
          ENDIF
          muParam=materialsInterpPoint%values(1,NO_PART_DERIV)
          rhoParam=materialsInterpPoint%values(2,NO_PART_DERIV)
        CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          !Calculate conductivity tensor
          IF(ASSOCIATED(fibreField)) &
            & CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & materialsInterpPoint%values(1:NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV)),intraConductivity,err,error,*999)
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+1: &
            & 2*NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),extraConductivity,err,error,*999)
          conductivity(1:numberOfXi,1:numberOfXi)=intraConductivity(1:numberOfXi,1:numberOfXi)+ &
            & extraConductivity(1:numberOfXi,1:numberOfXi)
        CASE DEFAULT
          !Do nothing
        END SELECT
        
        !Calculate Jacobian and Gauss weight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
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
          CALL Basis_QuadratureSchemeGet(rowsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowsQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowsPhi,err,error,*999)
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowsdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(equationsMatrix%updateMatrix) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & colsPhi,err,error,*999)
                  DO xiIdx1=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx1),gaussPointIdx,colsdPhidXi(xiIdx1),err,error,*999)
                  ENDDO !xiIdx1
                  sum=0.0_DP
                  SELECT CASE(esSpecification(3))
                  CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx                    
                  CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    sum=sum+aParam*rowsPhi*colsPhi
                  CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
                    & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        sum=sum+rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)*geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    sum=sum*diffCoeff1
                  CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  CASE DEFAULT
                    localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                      & " is not valid for a Laplace equation type of a classical field equations set class."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  
                  equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !updateMatrix
            IF(updateSource) THEN
              SELECT CASE(esSpecification(3))
              CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
                sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+  &
                  & sourceValue*rowsPhi*jacobianGaussWeight           
              CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
                & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
                sum=0.0_DP
                sum2=0.0_DP
                b=0.0_DP
                IF(equationsSet%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE) THEN
                  IF(inside) THEN
                    DO componentIdx1=1,numberOfDimensions
                      !Dynamic term (should be zero for static problems)
                      b(componentIdx1)=-rhoParam*((uValue(componentIdx1)-uOld(componentIdx1))/deltaT)
                      !Second viscous term
                      DO xiIdx1=1,numberOfXi
                        DO componentIdx2=1,numberOfDimensions
                          !ALE term
                          b(componentIdx1)=b(componentIdx1)+rhoParam*wValue(componentIdx2)* &
                            & uDeriv(componentIdx1,xiIdx1)*dXidX(xiIdx1,componentIdx2)
                          DO xiIdx2=1,numberOfXi
                            b(componentIdx1)=b(componentIdx1)+muParam*uSecond(componentIdx1,xiIdx1,xiIdx2)* &
                              & dXidX(xiIdx1,componentIdx2)*dXidX(xiIdx2,componentIdx2)
                          ENDDO !xiIdx2
                        ENDDO !componentIdx2
                      ENDDO !xiIdx1
                      DO xiIdx1=1,numberOfXi
                        sum=sum+b(componentIdx1)*rowsdPhidXi(xiIdx1)*dXidX(xiIdx1,componentIdx1)
                      ENDDO !xiIdx1
                    ENDDO !componentIdx1
                  ENDIF
                ELSE
                  IF(inside) THEN
                    DO componentIdx1=1,numberOfDimensions
                      !Dynamic term (should be zero for static problems)
                      b(componentIdx1)=-rhoParam*((uValue(componentIdx1)-uOld(componentIdx1))/deltaT)
                      !Nonlinear term (should be zero for Stokes flow)
                      DO xiIdx1=1,numberOfXi
                        DO componentIdx2=1,numberOfDimensions
                          !Nonlinear/ALE term
                          b(componentIdx1)=b(componentIdx1)-rhoParam*(uValue(componentIdx2)-wValue(componentIdx2))* &
                            & uDeriv(componentIdx1,xiIdx1)*dXidX(xiIdx1,componentIdx2)
                          DO xiIdx2=1,numberOfXi
                            !Viscous term
                            b(componentIdx1)=b(componentIdx1)+muParam*uSecond(componentIdx1,xiIdx1,xiIdx2)* &
                              & dXidX(xiIdx1,componentIdx2)*dXidX(xiIdx2,componentIdx2)
                          ENDDO !xiIdx2
                        ENDDO !componentIdx2
                      ENDDO !xiIdx1
                      DO xiIdx1=1,numberOfXi
                        sum=sum+b(componentIdx1)*rowsdPhidXi(xiIdx1)*dXidX(xiIdx1,componentIdx1)
                      ENDDO !xiIdx1
                    ENDDO !componentIdx1
                  ENDIF
                  IF(ASSOCIATED(equationsSet%analytic)) THEN
                    x=0.0_DP
                    x(1:numberOfDimensions)=geometricInterpPoint%values(1:numberOfDimensions,NO_PART_DERIV)
                    SELECT CASE(analyticFunctionType
                    CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_1)
                      sum2=-4.0_DP*PI*PI/100.0_DP*(3.0_DP*SIN(2.0_DP*PI*x(1)/10.0_DP)*SIN(2.0_DP*PI*x(2)/10.0_DP)* &
                        & SIN(2.0_DP*PI*x(3)/10.0_DP)-6.0_DP*rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)**2+ &
                        & 8.0_DP*rhoParam*COS(2.0_DP*PI*x(1)/10.0_DP)**2*COS(2.0_DP*PI*x(3)/10.0_DP)**2- &
                        & 2.0_DP*rhoParam*COS(2.0_DP*PI*x(3)/10.0_DP)**2+2.0_DP*rhoParam*COS(2.0_DP*PI*X(1)/10.0_DP)**2* &
                        & COS(2.0_DP*PI*x(2)/10.0_DP)**2+4.0_DP*rhoParam*COS(2.0_DP*PI*x(2)/10.0_DP)**2- &
                        & 2.0_DP*rhoParam*COS(2.0_DP*PI*x(2)/10.0_DP)**2*COS(2.0_DP*PI*x(3)/10.0_DP)**2)
                      sum=sum-sum2*rowsPhi
                    CASE(EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2)
                      sum2=-12.0_DP*SIN(2.0_DP*PI*x(1)/10.0_DP)*PI*PI/100.0_DP*SIN(2.0_DP*PI*x(2)/10.0_DP)* &
                        & SIN(2.0_DP*PI*x(3)/10.0_DP)
                      sum=-sum2*rowsPhi
                    ELSE 
                      CALL FlagError("Not implemented.",err,error,*999)
                    ENDIF
                  ENDIF
                ENDIF
              CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
!!TODO eventually handle this as an additional linear matrix that is not mapped to the sovler matrix. 
                !Loop over field components
                columnElementDOFIdx=0
                !Loop over element columns
                DO columnComponentIdx=1,numberOfColsComponents
                  NULLIFY(colsDomain)
                  CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                  NULLIFY(colsDomainTopology)
                  CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                  NULLIFY(colsDomainElements)
                  CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                  NULLIFY(colsBasis)
                  CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                  CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
                  CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                  DO columnElementParameterIdx=1,numberOfColsElementParameters
                    columnElementDOFIdx=columnElementDOFIdx+1
                    DO xiIdx1=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx1),gaussPointIdx,colsdPhidXi(xiIdx1),err,error,*999)
                    ENDDO !xiIdx1
                    sum2=0.0_DP
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum2=sum2+intraConductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    
                    sum=sum+sum2*jacobianGaussWeight*Vm(columnElementDOFIdx)
                    
                  ENDDO !columnElementParameterIdx
                ENDDO !columnComponentIdx  
                
              CASE DEFAULT
                sum=0.0_DP
              END SELECT
              
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                & sum*jacobianGaussWeight
              
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
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                    & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                    & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update
    
    EXITS("Poisson_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Poisson_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE Poisson_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)
    
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element Jacobian evaluation on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,esSpecification(3),gaussPointIdx, &
      & rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,
    REAL(DP) :: aParam,bParam,jacobianGaussWeight,jacobianValue,uValue
    TYPE(BasisType), POINTER :: colsBasis,geometricBasis,rowsBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(QuadratureSchemeType), POINTER :: colsQuadratureScheme,geometricQuadratureScheme,rowsQuadratureScheme
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Poisson_FiniteElementJacobianEvaluate",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
      !div(\sigma.grad u) + s = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.u + s = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + b.u^2 + a.u + s = 0
      !OK
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.e^(b.u) + s = 0
      !OK
    CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
      !div(grad p) -  div( \mu div(grad u) - \rho (\del u/\del t - u . grad u ) ) = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      !div((\sigma_i+\sigma_e).grad \phi_e) + div(\sigma_i. grad V_m) = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    NULLIFY(equations)
    CALL EquationSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    NULLIFY(jacobianMatrix)
    CALL EquationsMatricesResidual_JacobianMatrixGet(residualVector,1,jacobianMatrix,err,error,*999)
    updateJacobian=jacobianMatrix%updateJacobian
    
    IF(updateJacobian) THEN
      
      NULLIFY(lhsMapping)
      CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      
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

      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(colsDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(colsDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
      NULLIFY(colsBasis)
      CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(colsQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(colsQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(dependentInterpParamters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
      
      NULLIFY(materialsInterpParamters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)

        !Find material parameters and u value at this Gauss point
        uValue=dependentInterpPoint%values(1,NO_PART_DERIV)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
          bParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+2,NO_PART_DERIV)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
          aParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+1,NO_PART_DERIV)
          bParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+2,NO_PART_DERIV)
        END SELECT
              
        !Calculate jacobianGaussWeight.
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
              
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
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & gaussPointIdx,rowsPhi,err,error,*999)
            columnElementDOFIdx=0
            IF(equationsMatrix%updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,colElementParameterIdx,NO_PART_DERIV, &
                    & colsPhi,err,error,*999)
                  SELECT CASE(esSpecification(3))
                  CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
                    jacobianValue=-2.0_DP*bParam*rowsPhi*colsPhi*uValue
                  CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
                    jacobianValue=-aParam*bParam*rowsPhi*colsPhi*EXP(bParam*uValue)
                  END SELECT
                  jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                    & jacobianValue*jacobianGaussWeight
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
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
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            columnElementDOFIdx=0
            !Loop over element columns
            DO columnComponentIdx=1,numberOfColsComponents
              NULLIFY(colsDomain)
              CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
              NULLIFY(colsDomainTopology)
              CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
              NULLIFY(colsDomainElements)
              CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
              NULLIFY(colsBasis)
              CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
              CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
              DO columnElementParameterIdx=1,numberOfColsElementParameters
                columnElementDOFIdx=columnElementDOFIdx+1
                jacobianMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  & jacobianMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                  & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                  & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
              ENDDO !columnElementParameterIdx
            ENDDO !columnComponentIdx
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
      
    ENDIF !update Jacobian
       
    EXITS("Poisson_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("Poisson_FiniteElementJacobianEvaluate",err,error)
    EXITS("Poisson_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE Poisson_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Poisson equation finite element equations set.
  SUBROUTINE Poisson_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnXiIdx,columnElementParameterIdx,componentIdx,gaussPointIdx, &
      & rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx
    REAL(DP) :: aParam,bParam,colsPhi(3),gaussWeight,jacobianGaussWeight,rowsPhi(3),sourceValue,sum1,sum2,uValue
    LOGICAL :: update
    TYPE(BasisType), POINTER :: colsBasis,geometricBasis,rowsBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField,sourceField
    TYPE(FieldVariableType), POINTER :: colsVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: colsQuadratureScheme,quadratureScheme,rowsQuadratureScheme
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Poisson_FiniteElementResidualEvaluate",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_POISSON_SUBTYPE)
      !div(\sigma.grad u) + s = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.u + s = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + b.u^2 + a.u + s = 0
      !OK
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
      !div(\sigma.grad u) + a.e^(b.u) + a = 0
      !OK
    CASE(EQUATIONS_SET_NONLINEAR_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & EQUATIONS_SET_ALE_PRESSURE_POISSON_SUBTYPE,EQUATIONS_SET_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE(EQUATIONS_SET_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      !div((\sigma_i+\sigma_e).grad \phi_e) + div(\sigma_i. grad V_m) = 0
      CALL FlagError("Can not evaluate a Jacobian for a linear Poisson equation.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(equations)
    CALL EquationSet_EquationsGet(equationsSet,equations,err,error,*999)
    NULLIFY(vectorEquations)
    CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
    NULLIFY(vectorMapping)
    CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
    NULLIFY(lhsMapping)
    CALL EquationsMappingVector_LHSMappingGet(vectorMapping,lhsMapping,err,error,*999)
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    NULLIFY(linearMapping)
    CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
    NULLIFY(residualMapping)
    CALL EquationsMappingNonlinear_ResidualMappingGet(nonlinearMapping,1,residualMapping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    NULLIFY(sourceMapping)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
    ENDIF
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(nonlinearMatrices)
    CALL EquationsMatricesVector_NonlinearMatricesGet(vectorMatrices,nonlinearMatrices,err,error,*999)
    NULLIFY(residualVector)
    CALL EquationsMatricesNonlinear_ResidualVectorGet(nonlinearMatrices,1,residualVector,err,error,*999)
    updateResidual=residualVector%updateVector
    NULLIFY(linearMatrices)
    CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
    NULLIFY(equationsMatrix)
    CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,equationsMatrix,err,error,*999)
    updateMatrix=equationsMatrix%updateMatrix
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      updateSource=sourceVector%updateVector
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,1,rhsVector,err,error,*999)
      updateRHS=rhsVector%updateVector
    ENDIF

    update=(updateMatrix.OR.updateSource.OR.updateResidual.OR.updateRHS)

    IF(update) THEN
      
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
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
      NULLIFY(dependentDecomposition)
      CALL Field_DecompositionGet(dependentField,dependentDecomposition,err,error,*999)
      NULLIFY(colsDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,colsDomain,err,error,*999)
      NULLIFY(colsDomainTopology)
      CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
      NULLIFY(colsDomainElements)
      CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
      NULLIFY(colsBasis)
      CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
    
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(colsQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
      CALL BasisQuadrature_NumberOfGaussGet(colsQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(dependentInterpParamters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,colsVariableType,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
    
      NULLIFY(materialsInterpParamters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(sourceInterpParamters)
      NULLIFY(sourceInterpPoint)
      IF(ASSOCIATED(sourceField)) THEN
        CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & sourceInterpParameters,err,error,*999)
        CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
      ENDIF
        
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      IF(ASSOCIATED(fibreField)) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF
            
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint,err,error,*999)
        IF(ASSOCIATED(fibreField)) &
          & CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint,err,error,*999)
        IF(ASSOCIATED(sourceField)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint,err,error,*999)
          sourceValue=sourceInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        
        uValue=dependentInterpPoint%values(1,NO_PART_DERIV)
        aParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+1,NO_PART_DERIV)
        bParam=materialsInterpPoint%values(NUMBER_OF_VOIGT(numberOfDimensions)+2,NO_PART_DERIV)
        
        CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
          & materialsInterpPoint%values(1:NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)        
        
        !Calculate jacobianGaussWeight.      
!!TODO: Think about symmetric problems.
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
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
          CALL Basis_QuadratureSchemeGet(rowsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,rowsQuadratureScheme,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
              & rowsPhi,err,error,*999)
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowsdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            columnElementDOFIdx=0
            IF(updateMatrix) THEN
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_QuadratureSchemeGet(colsBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,colsQuadratureScheme,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                    & colsPhi,err,error,*999)
                  DO xiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                      & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,colsdPhidXi(xiIdx),err,error,*999)
                  ENDDO !xiIdx
                  sum=0.0_DP
                  SELECT CASE(esSpecification(3))
                  CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    sum=sum+aParam*rowsPhi*colsPhi
                  CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
                    DO rowIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowsdPhidXi(rowXiIdx)*colsdPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                  CASE DEFAULT
                    localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                      & " is not valid for a Laplace equation type of a classical field equations set class."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  
                  equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+sum*jacobianGaussWeight
                  
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !updateMatrix
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                & sourceValue*rowsPhi*jacobianGaussWeight
              & quadratureScheme%gaussBasisFunctions(rowElementParameterIdx,NO_PART_DERIV,gaussPointIdx)*aParam*
            ENDIF
            IF(updateResidual) THEN
              SELECT CASE(esSpecification(3))
              CASE(EQUATIONS_SET_QUADRATIC_SOURCE_POISSON_SUBTYPE)
                residualVector%elementVector%vector(rowElementDOFIdx)=residualVector%elementVector%vector(rowElementDOFIdx)+ &
                  & bParam*uValue**2*rowsPhi*jacobianGaussWeight
              CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_POISSON_SUBTYPE)
!!TODO: Handle floating point exceptions better
                IF((bParam*uValue)>20000.0_DP) THEN
                  localError="The value of "//TRIM(NumberToVString(bParam*uValue,"*",err,error))// &
                    & " is out of range for an exponential function."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                residualVector%elementVector%vector(rowElementDOFIdx)=residualVector%elementVector%vector(rowElementDOFIdx)+ &
                  & aParam*EXP(bParam*uValue)*rowsPhi*jacobianGaussWeight
              CASE DEFAULT
                localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                  & " is not valid for a Laplace equation type of a classical field equations set class."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
            ENDIF
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
        DO rowComponentIdx=1,numberOfRowsComponents
          NULLIFY(rowsDomain)
          CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowsDomain,err,error,*999)
          NULLIFY(rowsDomainTopology)
          CALL Domain_DomainTopologyGet(rowsDomain,rowsDomainTopology,err,error,*999)
          NULLIFY(rowsDomainElements)
          CALL DomainTopology_DomainElementsGet(rowsDomainTopology,rowsDomainElements,err,error,*999)
          NULLIFY(rowsBasis)
          CALL DomainElements_ElementBasisGet(rowsDomainElements,elementNumber,rowsBasis,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(rowsBasis,numberOfRowsElementParameters,err,erorr,*999)
          DO rowElementParameterIdx=1,numberOfRowsElementParameters
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrix) THEN
              columnElementDOFIdx=0
              !Loop over element columns
              DO columnComponentIdx=1,numberOfColsComponents
                NULLIFY(colsDomain)
                CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,colsDomain,err,error,*999)
                NULLIFY(colsDomainTopology)
                CALL Domain_DomainTopologyGet(colsDomain,colsDomainTopology,err,error,*999)
                NULLIFY(colsDomainElements)
                CALL DomainTopology_DomainElementsGet(colsDomainTopology,colsDomainElements,err,error,*999)
                NULLIFY(colsBasis)
                CALL DomainElements_ElementBasisGet(colsDomainElements,elementNumber,colsBasis,err,error,*999)
                CALL Basis_NumberOfElementParametersGet(colsBasis,numberOfColsElementParameters,err,error,*999)
                DO columnElementParameterIdx=1,numberOfColsElementParameters
                  columnElementDOFIdx=columnElementDOFIdx+1
                  equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                    & equationsMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                    & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                    & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateResidual) THEN
              residualVector%elementVector%vector(rowElementDOFIdx)=residualVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update
    
    EXITS("Poisson_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Poisson_FiniteElementResidualEvaluate",err,error)
    EXITS("Poisson_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Poisson_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Performs pre-solver operations for the Poisson problem.
  SUBROUTINE Poisson_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,currentIteration,maximumNumberOfIterations,outputType,pSpecification(3), &
      & solverGlobalNumber
    REAL(DP) :: absoluteTolerance,relativeTolerance
    LOGICAL :: continueLoop
    TYPE(ControlLoopType), POINTER :: controlLoop 
    TYPE(EquationsSetType), POINTER :: equationsSet 
    TYPE(ProblemType), POINTER :: problem 
    TYPE(SolverEquationsType), POINTER :: solverEquations 
    TYPE(SolverMappingType), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Poisson_PreSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
  
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL ControlLoop_WhileInformationGet(controlLoop,currentIteration,maximumNumberOfIterations,absoluteTolerance, &
        & relativeTolerance,continueLoop,err,error,*999)
      IF(currentIteration==1)THEN
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquation_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        IF(ASSOCIATED(equationsSet%analytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2) THEN
            !do nothing
          ELSE
            IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
            ENDIF
            !Update indpendent data fields
            CALL Poisson_PreSolveUpdateInputData(CONTROLLOOP,solver,err,error,*999)
            IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
              CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN                  
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
          ENDIF
          !Update indpendent data fields
          CALL Poisson_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
          IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
            CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
          ENDIF
        ENDIF
      ELSE
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN       
          CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
        ENDIF
      ENDIF
    CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGloblaNumber,err,error,*999)
      !Pre solve for the linear solver: fitting problem
      IF(solverGlobalNumber==1) THEN
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
        ENDIF
        !Update independent data fields
        CALL Poisson_PreSolveUpdateInputData(controlLoop,solver,err,error,*999)
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Solving fitting problem... ",err,error,*999)
        ENDIF
      ELSE IF(solverGlobalNumber==2) THEN
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Update source field ... ",err,error,*999)
        ENDIF
        !First update source field for PPE calculation
        CALL Poisson_PreSolveUpdatePPESource(solver,err,error,*999)
        !Read in the lable information
        !CALL Poisson_PreSolveUpdateInputData(solver,err,error,*999)
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Solving PPE problem... ",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver global number not associated for PPE problem.",err,error,*999)
      ENDIF
    CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL ControlLoop_WhileInformationGet(controlLoop,currentIteration,maximumNumberOfIterations,absoluteTolerance, &
        & relativeTolerance,continueLoop,err,error,*999)
      IF(currentIteration==1)THEN
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquation_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        IF(ASSOCIATED(equationsSet%analytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          IF(analyticFunctionType==EQUATIONS_SET_PRESSURE_POISSON_THREE_DIM_2) THEN
            !do nothing
          ELSE
            IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
            ENDIF
            !Update indpendent data fields
            CALL Poisson_PreSolveUpdateInputData(solver,err,error,*999)
            IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
              CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
            ENDIF
          ENDIF
        ELSE
          IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Read in vector data... ",err,error,*999)
          ENDIF
          !Update indpendent data fields
          CALL Poisson_PreSolveUpdateInputData(solver,err,error,*999)
          IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
            CALL WriteString(GENERAL_OUTPUT_TYPE,"While loop... ",err,error,*999)
          ENDIF
        ENDIF
        !Update mesh
        CALL Poisson_PreSolveUpdatePPEMesh(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Poisson_PreSolve")
    RETURN
999 ERRORSEXITS("Poisson_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PreSolve
  
  !
  !================================================================================================================================
  !
 

  !>Update boundary conditions for Poisson pre solve
  SUBROUTINE Poisson_PreSolveUpdateInputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDimensions,currentIteration,inputIteration,inputType,inputOption,outputIteration,outputType, &
      & solverGlobalNumber
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    REAL(DP), POINTER :: inputVelocityNewData(:),inputVelocityOldData(:),inputVelocityLabelData(:),inputVelocityUData(:), &
      & inputVelocityVData(:),inputVelocityWData(:)
    LOGICAL :: boundaryUpdate
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsType), POINTER :: equations
    TYPE(FieldType), POINTER :: dependentField,geometricField,sourceField
    TYPE(SolverEquationsType), POINTER :: solverEquations 
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    boundaryUpdate=.FALSE.

    ENTERS("Poisson_PreSolveUpdateInputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    NULLIFY(sourceField)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==1) THEN
        !only read in data if it is the fitting problem
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
          & outputIteration,inputIteration,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(sourceField)
        CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
        inputType=1
        inputOption=1
        CALL Field_ParameterSetDataGet(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,inputVelocityNewData,err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityNewData,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)
      ELSE IF(solverGlobalNumber==2) THEN
        !this is the interior flag
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
          & outputIteration,inputIteration,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(dependentField)
        CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        !Set pressure field to zero     
        !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
        inputType=1
        inputOption=3
        CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_LABEL_SET_TYPE,inputVelocityLabelData, &
          & err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityLabelData,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)
        !this is the reference U velocity
        inputType=1
        inputOption=4
        CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL1_SET_TYPE,inputVelocityUData, &
          & err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityUData,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)
        !this is the reference V velocity
        inputType=1
        inputOption=5
        CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL2_SET_TYPE,inputVelocityVData, &
          & err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityVData,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)
        !this is the reference W velocity
        inputType=1
        inputOption=6
        CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL3_SET_TYPE,inputVelocityWData, &
          & err,error,*999)
        CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityWData,numberOfDimensions,inputType,inputOption, &
          & currentIteration,1.0_DP,err,error,*999)
      ENDIF
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
        CALL WriteString(GENERAL_OUTPUT_TYPE,"Read input data... ",err,error,*999)
      ENDIF
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(sourceField)
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
      inputType=1
      inputOption=1
      CALL Field_ParameterSetDataGet(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,inputVelocityNewData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityNewData,numberOfDimensions,inputType,inputOption, &
        & CONTROL_TIME_LOOP%timeLoop%iterationNumber,1.0_DP,err,error,*999)
      !this is the previous time step
      !\todo: Provide possibility for user to define input type and option (that's more or less an IO question)
      inputType=1
      inputOption=2
      CALL Field_ParameterSetDataGet(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,inputVelocityOldData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityOldData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      !this is the interior flag
      inputType=1
      inputOption=3
      CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_LABEL_SET_TYPE,inputVelocityLabelData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityLabelData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      !this is the reference U velocity
      inputType=1
      inputOption=4
      CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL1_SET_TYPE,inputVelocityUData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityUData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      !this is the reference V velocity
      inputType=1
      inputOption=5
      CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL2_SET_TYPE,inputVelocityVData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityVData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      !this is the reference W velocity
      inputType=1
      inputOption=6
      CALL Field_ParameterSetDataGet(dependentField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_VEL3_SET_TYPE,inputVelocityWData, &
        & err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputVelocityWData,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    IF(ASSOCIATED(sourceField)) THEN    
      CALL Field_ParameterSetUpdateStart(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateStart(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
    ENDIF
      
    EXITS("Poisson_PreSolveUpdateInputData")
    RETURN
999 ERRORSEXITS("Poisson_PreSolveUpdateInputData",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PreSolveUpdateInputData

  !
  !================================================================================================================================
  !
  
  !>Update mesh velocity and move mesh for ALE PPE problem
  SUBROUTINE Poisson_PreSolveUpdatePPEMesh(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfDimensions,geometricMeshComponent
    INTEGER(INTG) :: inputIteration,inputType,inputOption,componentIdx,derivativeIdx,localDOFIdx,nodeIdx,variableIdx,variable_type
    REAL(DP) :: currentTime,timeIncrement,alpha
    REAL(DP), POINTER :: meshDisplacementValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop to solve.
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: independentField
    TYPE(FieldVariableType), POINTER :: fieldVariable
    TYPE(SolverType), POINTER :: solverALEPPE
    TYPE(SolverEquationsType), POINTER :: solverEquations  !<A pointer to the solver equations
    TYPE(SolverMappingType), POINTER :: solverMapping !<A pointer to the solver mapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poisson_PreSolveUpdatePPEMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      !Update mesh within the dynamic solver
      IF(solver%solveType/=SOLVER_LINEAR_TYPE) &
        & CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
      !Get the independent field for the ALE PPE problem
      NULLIFY(solvers)
      CALL Solver_SolversGet(solver,solvers,err,error,*999)
      NULLIFY(solverALEPPE)
      CALL Solvers_SolverGet(solvers,1,solverALEPPE,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquations(solverALEPPE,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(independentField)
      CALL EquationsSet_IndependentFieldGet(equationsSet,independentField,err,error,*999)
      !Get the data
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      !\todo: Introduce flags set by the user (42/1 only for testings purpose)
      !Copy input to PPE' independent field
      inputType=42
      inputOption=1
      NULLIFY(meshDisplacementValues)
      CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & meshDisplacementValues,err,error,*999)
      CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,meshDisplacementValues,numberOfDimensions,inputType,inputOption, &
        & currentIteration,1.0_DP,err,error,*999)
      CALL Field_ParameterSetUpdateStart(independentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(independentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE,err,error,*999)
      !Use calculated values to update mesh
      CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
      NULLIFY(equations)
      CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      NULLIFY(vectorMapping)
      CALL VectorEquations_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
      CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
      DO variableIdx=1,numberOfVariables
        NULLIFY(fieldVariable)
        CALL Field_VariableIndexGet(dependentField,variableIdx,fieldVariable,variableType,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(fieldVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          NULLIFY(domain)
          CALL FieldVariable_ComponentDomainGet(fieldVariable,componentIdx,domain,err,error,*999)
          NULLIFY(domainTopology)
          CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
          NULLIFY(domainNodes)
          CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
          !Loop over the local nodes excluding the ghosts.
          CALL DomainNodes_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
          DO nodeIdx=1,numberOfNodes
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodesIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(fieldVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
              CALL Field_ParameterSetAddLocalDOF(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,localDOFIdx, & 
                & meshDisplacementValues(localDOFIdx),err,error,*999)
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !componentIdx
      ENDDO !variableIdx
      CALL Field_ParameterSetUpdateStart(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL Field_ParameterSetUpdateFinish(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
      !Now use displacement values to calculate velocity values
      alpha=1.0_DP/timeIncrement
      CALL Field_ParameterSetsCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
        & FIELD_MESH_VELOCITY_SET_TYPE,alpha,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a PPE equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Poisson_PreSolveUpdatePPEMesh")
    RETURN
999 ERRORSEXITS("Poisson_PreSolveUpdatePPEMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PreSolveUpdatePPEMesh

  !
  !================================================================================================================================
  !
  
  !>Update source for fitted PPE problem
  SUBROUTINE Poisson_PreSolveUpdatePPESource(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: I,numberOfDimensionsPPE,numberOfDimensionsFitted,geometricMeshComponent,inputIteration
    REAL(DP) :: currentTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(SolverType), POINTER :: solverFitted, solverPPE
    TYPE(FieldType), POINTER :: dependentFieldFitted,sourceFieldPPE
    TYPE(SolverEquationsType), POINTER :: solverEquationsFitted,solverEquationsPPE
    TYPE(SolverMappingType), POINTER :: solverMappingFitted,solverMappingPPE
    TYPE(EquationsSetType), POINTER :: equationsSetFitted,equationsSetPPE
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poisson_PreSolveUpdatePPESource",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_ProblemSpecificationGet(problem,pSpecification,err,error,*999)
    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
      & outputIteration,inputIteration,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverGlobalNumber,err,error,*999)
      IF(solverGlobalNumber==1) THEN
        !do nothing for the fitting process
      ELSE IF(solverGlobalNumber==2) THEN
        !Get the dependent field for the three component Laplace problem
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)        
        CALL Solvers_SolverGet(solvers,1,solverFitted,err,error,*999)
        NULLIFY(solverEquationsFitted)
        CALL Solver_SolverEquationsGet(solverFitted,solverEquationsFitted,err,error,*999)
        NULLIFY(solverMappingFitted)
        CALL SolverEquations_SolverMappingGet(solverEquationsFitted,solverMappingFitted,err,error,*999)
        NULLIFY(equationsSetFitted)
        CALL SolverMapping_EquationsSetGet(solverMappingFitted,1,equationsSetFitted,err,error,*999)
        NULLIFY(geometricFieldFitted)
        CALL EquationsSet_GeometricFieldGet(equationsSetFitted,geometricFieldFitted,err,error,*999)
        NULLIFY(dependentFieldFitted)
        CALL EquationsSet_DependentFieldGet(equationsSetFitted,dependentFieldFitted,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricFieldFitted,FIELD_U_VARIABLE_TYPE,numberOfDimensionsFitted,err,error,*999)
        !Get the source field for the PPE problem
        CALL Solvers_SolverGet(solvers,2,solverPPE,err,error,*999)
        NULLIFY(solverEquationsPPE)
        CALL Solver_SolverEquationsGet(solverPPE,solverEquationsPPE,err,error,*999)
        NULLIFY(solverMappingPPE)
        CALL SolverEquations_SolverMappingGet(solverEquationsPPE,solverMappingPPE,err,error,*999)
        NULLIFY(equationsSetPPE)
        CALL SolverMapping_EquationsSetGet(solverMappingPPE,1,equationsSetPPE,err,error,*999)
        NULLIFY(geometricFieldPPE)
        CALL EquationsSet_GeometricFieldGet(equationsSetPPE,geometricFieldPPE,err,error,*999)
        NULLIFY(sourceFieldPPE)
        CALL EquationsSet_SourceFieldGet(equationsSetPPE,sourceFieldPPE,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricFieldPPE,FIELD_U_VARIABLE_TYPE,numberOfDimensionsPPE,err,error,*999)
        !Copy result from FITTING to PPE's source field
        IF(numberOfDimensionsPPE/=numberOfDimensionsFitted) &
          & CALL FlagError("Mesh update is not defined for non-dynamic problems.",err,error,*999)
        IF(currentIteration/=1) THEN
          DO componentIdx=1,numberOfDimensionsPPE
            IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
              CALL WriteString(GENERAL_OUTPUT_TYPE,"Save old source field... ",err,error,*999)
            ENDIF
            CALL Field_ParametersToFieldParametersCopy(sourceFieldPPE,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
              & componentIdx,sourceFieldPPE,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,componentIdx,err,error,*999)
          ENDDO !componentIdx
        ENDIF
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Update new source field... ",err,error,*999)
        ENDIF
        DO componentIdx=1,numberOfDimensionsPPE
          CALL Field_ParametersToFieldParametersCopy(dependentFieldFitted, & 
            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,sourceFieldPPE, & 
            & FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE,componentIdx,err,error,*999)
        ENDDO !componentIdx
        IF(controlLoop%timeLoop%iterationNumber==1) THEN
          IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
            CALL WriteString(GENERAL_OUTPUT_TYPE,"Old source field is new source field!... ",err,error,*999)
          ENDIF
          DO componentIdx=1,numberOfDimensionsPPE
            CALL Field_ParametersToFieldParametersCopy(sourceFieldPPE,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA1_SET_TYPE, &
              & componentIdx,sourceFieldPPE,FIELD_U_VARIABLE_TYPE,FIELD_INPUT_DATA2_SET_TYPE,componentIdx,err,error,*999)
          ENDDO !componentIdx
        ENDIF
        !Use calculated values to update mesh
        CALL Field_ComponentMeshComponentGet(equationsSetPPE%geometry%geometricField, & 
          & FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
        CALL Field_ParameterSetUpdateStart(equationsSetPPE%source%sourceField,FIELD_U_VARIABLE_TYPE, & 
          & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(equationsSetPPE%source%sourceField,FIELD_U_VARIABLE_TYPE, & 
          & FIELD_INPUT_DATA1_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateStart(equationsSetPPE%source%sourceField,FIELD_U_VARIABLE_TYPE, & 
          & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(equationsSetPPE%source%sourceField,FIELD_U_VARIABLE_TYPE, & 
          & FIELD_INPUT_DATA2_SET_TYPE,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a PPE equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Poisson_PreSolveUpdatePPESource")
    RETURN
999 ERRORSEXITS("Poisson_PreSolveUpdatePPESource",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PreSolveUpdatePPESource

  !
  !================================================================================================================================
  !

  !>Performs post-solve operations for a Poisson problem
  SUBROUTINE Poisson_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: outputType,pSpecification(3),solverGlobalNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poisson_PostSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    NULLIFY(solver2)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      CALL Poisson_PostSolveOutputData(solver,err,error,*999)
    CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL Solver_GlobalNumberGet(solver,solverGloblaNumber,err,error,*999)
      !Post solve for the linear solver
      IF(solverGlobalNumber==1) THEN
        !do nothing for fitting problems
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN        
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Transferring fitting results to PPE... ",err,error,*999)
        ENDIF
      ELSE IF(solverGlobalNumber==2) THEN
        IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN                  
          CALL WriteString(GENERAL_OUTPUT_TYPE,"PPE post solve... ",err,error,*999)
        ENDIF
        CALL Poisson_PostSolveOutputData(solver,err,error,*999)
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
 
    EXITS("Poisson_PostSolve")
    RETURN
999 ERRORSEXITS("Poisson_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PostSolve

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE Poisson_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equationsSetIdx,currenIteration,inputIteration,numberOfEquationsSets,outputIteration,outputType, &
      & pSpecification(3)
    REAL(DP) :: absoluteTolerance,currentTime,relativeTolerance,startTime,stopTime,timeIncrement
    LOGICAL :: continueLoop,exportField
    CHARACTER(14) :: file,outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(VARYING_STRING) :: localError,method

    ENTERS("Poisson_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
      ! do nothing
    CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE,PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_WhileInformationGet(controlLoop,currentIteration,maximumNumberOfIterations,absoluteTolerance, &
        & relativeTolerance,continueLoop,err,error,*999)
      IF(currentIteration==maximumNumberOfIterations)THEN
        CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
          & outputIteration,inputIteration,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Make sure the equations sets are up to date
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSet,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          IF(outputIteration/=0) THEN
            IF(currentTime<stopTime) THEN
              IF(currentIteration<10) THEN
                WRITE(outputFile,'("TIME_STEP_000",I0)') currentIteration
              ELSE IF(currentIteration<100) THEN
                WRITE(outputFile,'("TIME_STEP_00",I0)') currentIteration
              ELSE IF(currentIteration<1000) THEN
                WRITE(outputFile,'("TIME_STEP_0",I0)') currentIteration
              ELSE IF(currentIteration<10000) THEN
                WRITE(outputFile,'("TIME_STEP_",I0)') currentIteration
              END IF
              file=outputFile
              method="FORTRAN"
              exportField=.TRUE.
              IF(exportField) THEN          
                IF(MOD(currentIteration,outputIteration)==0)  THEN
                  IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                  ENDIF
                  CALL FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK(equationsSet%region,equationsSet%globalNumber, &
                    & outputFile,err,error,*999)
                  IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                    CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                    CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx
      ENDIF
    CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
      CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSet,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(outputIteration/=0) THEN
          IF(currentTime<=stopTime) THEN
            IF(currentIteration<10) THEN
              WRITE(outputFile,'("TIME_STEP_000",I0)') currentIteration
            ELSE IF(currentIteration<100) THEN
              WRITE(outputFile,'("TIME_STEP_00",I0)') currentIteration
            ELSE IF(currentIteration<1000) THEN
              WRITE(outputFile,'("TIME_STEP_0",I0)') currentIteration
            ELSE IF(currentIteration<10000) THEN
              WRITE(outputFile,'("TIME_STEP_",I0)') currentIteration
            END IF
            file=outputFile
            method="FORTRAN"
            exportField=.TRUE.
            IF(exportField) THEN          
              IF(MOD(currentIteration,outputIteration)==0)  THEN   
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                ENDIF
                CALL FLUID_MECHANICS_IO_WRITE_ENCAS_BLOCK(equationsSet%region,equationsSet%globalNumber, &
                  & outputFile,err,error,*999)
                IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                  CALL WriteString(GENERAL_OUTPUT_TYPE,outputFile,err,error,*999)
                  CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Poisson equation fluid type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Poisson_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Poisson_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Poisson equation type.
  SUBROUTINE Poisson_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poisson_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE, &
      & PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
      & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Poisson problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_POISSON_EQUATION_TYPE,problemSubtype]

    EXITS("Poisson_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Poisson_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Poisson equations problem.
  SUBROUTINE Poisson_ProblemSetup(problem,problemSetup,err,error,*)

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
   
    ENTERS("Poisson_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
      & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
      & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid."
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
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          !Do nothing as the default is a simple loop
        CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
            !Set up a time control loop with a while subloop
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL ControlLoop_NumberOfSubLoopsSet(controlLoop,1,err,error,*999)
          CALL ControlLoop_SubLoopGet(controlLoop,1,subLoop,err,error,*999)
          CALL ControlLoop_TypeSet(subLoop,CONTROL_WHILE_LOOP_TYPE,err,error,*999)
          !Set the number of iterations for the while loop to 3 for now \todo make more general
          CALL ControlLoop_MaximumIterationsSet(subLoop,1,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
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
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          !Get the sub loop
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoopRoot,1,subLoop,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(subLoop,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Start the solvers creation
        NULLIFY(solvers)
        CALL Solvers_CreateStart(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification)
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &          
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE, &
          & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          !Create one linear solver
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be a linear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Start the linear solver creation
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
          !Create one nonlinear solver
          CALL Solvers_NumberOfSolversSet(solvers,1,err,error,*999)
          !Set the solver to be nonlinear solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
          !Create two solver
          CALL Solvers_NumberOfSolversSet(solvers,2,err,error,*999)
          !Set the first solver to be a linear solver for the velocity field fitting
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Start the linear solver creation
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
          !Set the second solver to be a linear solver for the pressure Poisson equation
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Start the linear solver creation
          CALL Solver_TypeSet(solver,SOLVER_LINEAR_TYPE,err,error,*999)
          !Set solver defaults
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          !Get the sub loop
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoopRoot,1,subLoop,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(subLoop,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Get the solvers
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          !Get the sub loop
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoopRoot,1,subLoop,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(subLoop,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Get the solver
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        SELECT CASE(pSpecification)
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &          
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE, &
          & PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          SELECT CASE(pSpecification(3))
          CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
            & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
            & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
            & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
            CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE)
          !Get the solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE(PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE)
          !Get the first solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEqutions,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
          !Get the second solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,2,solver,err,error,*999)
          !Create the solver equations
          NULLIFY(solverEquations)
          CALL SolverEquations_CreateStart(solver,solverEquations,err,error,*999)
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_QUASISTATIC,err,error,*999)
          CALL SolverEquations_SparsityTypeSet(solverEquations,SOLVER_SPARSE_MATRICES,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the control loop
        NULLIFY(controlLoopRoot)
        CALL Problem_ControlLoopRootGet(problem,controlLoopRoot,err,error,*999)
        SELECT CASE(pSpecification(3))
        CASE(PROBLEM_LINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_POISSON_SUBTYPE, &
          & PROBLEM_FITTED_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_EXTRACELLULAR_BIDOMAIN_POISSON_SUBTYPE)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(controlLoopRoot,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE(PROBLEM_LINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_NONLINEAR_PRESSURE_POISSON_SUBTYPE, &
          & PROBLEM_ALE_PRESSURE_POISSON_SUBTYPE)
          !Get the sub loop
          NULLIFY(subLoop)
          CALL ControlLoop_SubLoopGet(controlLoopRoot,1,subLoop,err,error,*999)
          NULLIFY(controlLoop)
          CALL ControlLoop_Get(subLoop,CONTROL_LOOP_NODE,controlLoop,err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Get the solver equations
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        CALL Solvers_NumberOfSolversGet(solvers,numberOfSolvers,err,error,*999)
        DO solverIdx=1,numberOfSolvers
          !Get the solver
          NULLIFY(solver)
          CALL Solvers_SolverGet(solvers,solverIdx,solver,err,error,*999)
          NULLIFY(solverEquations)
          CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
          !Finish the solver equations creation
          CALL SolverEquations_CreateFinish(solverEquations,err,error,*999)
        ENDDO !solverIdx
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a Poisson equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a Poisson equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Poisson_ProblemSetup")
    RETURN
999 ERRORSEXITS("Poisson_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Poisson_ProblemSetup

  !
  !================================================================================================================================
  !

END MODULE PoissonEquationsRoutines
