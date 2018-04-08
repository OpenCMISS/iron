!> \file
!> \author Chris Bradley
!> \brief This module handles all diffusion equation routines.
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

!>This module handles all diffusion equation routines.
MODULE DIFFUSION_EQUATION_ROUTINES

  USE ANALYTIC_ANALYSIS_ROUTINES
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE MeshAccessRoutines
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE SolverAccessRoutines
  USE SolverMappingAccessRoutines
  USE Timer
  USE Types

  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE
  
  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Diffusion_AnalyticFunctionsEvaluate,Diffusion_BoundaryConditionAnalyticCalculate

  PUBLIC DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP  

  PUBLIC DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC Diffusion_EquationsSetSolutionMethodSet
  
  PUBLIC Diffusion_EquationsSetSpecificationSet

  PUBLIC DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  PUBLIC Diffusion_FiniteElementJacobianEvaluate
  
  PUBLIC Diffusion_FiniteElementResidualEvaluate
  
  PUBLIC Diffusion_ProblemSpecificationSet

  PUBLIC DIFFUSION_EQUATION_PROBLEM_SETUP

  PUBLIC Diffusion_PreSolve,Diffusion_PostSolve

!!TODO: should the following two routines really be public???

  PUBLIC Diffusion_PreSolveGetSourceValue
  
  PUBLIC Diffusion_PreSolveStoreCurrentSolution


CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a two-dimensional unsteady heat equation solution (diffusion coefficient is 1)
  SUBROUTINE Diffusion_BoundaryConditionAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to calculate the boundary conditions for
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions !<A pointer to the boundary conditions to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,derivativeIdx,dimensionIdx,globalDerivativeIndex,imyMatrixIdx, &
      & localDofIdx,nodeIdx,numberOfDimensions,variableIdx,variableType,versionIdx
    INTEGER(INTG), POINTER :: equationsSetParameters(:)
    REAL(DP) :: initialValue,normal(3),tangents(3,3),time,VALUE,x(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(FIELD_TYPE), POINTER :: analyticField,dependentField,equationsSetField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dependentVariable,geometricVariable
 
    ENTERS("Diffusion_BoundaryConditionAnalyticCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a diffusion type equations set.", &
      & err,error,*999)
    CALL EquationsSet_AnalyticCreated(equationsSet,err,error,*999)
    
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
    NULLIFY(geometricParameters)
    CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(analyticField)
    NULLIFY(analyticParameters)    
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    IF(ASSOCIATED(analyticField)) &
      & CALL Field_ParameterSetDataGet(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    NULLIFY(materialsField)
    NULLIFY(materialsParameters)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    IF(ASSOCIATED(materialsField)) &
      CALL Field_ParameterSetDataGet(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)           
    analyticFunctionType=equationsSet%analytic%ANALYTIC_FUNCTION_TYPE
    CALL EquationsSet_AnalyticTimeGet(equationsSet,time,err,error,*999)
    IF(equationsSet%specification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)THEN
      !If a multi-comp model, we will use the equations set field information to assign only the appropriate field variable boundary conditions
      !Use predetermined mapping from equations set field compartment number to field variable type
      NULLIFY(equationsSetField)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      NULLIFY(equationsSetParameters)
      CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetParameters, &
        & err,error,*999)
      imyMatrixIdx = equationsSetParameters(1)
      DO variableIdx=0,1
        variableType=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(imyMatrixIdx-1))+variableIdx
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
        CALL Field_ParameterSetEnsureCreated(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        DO componentidx=1,dependentVariable%NUMBER_OF_COMPONENTS
          IF(dependentVariable%components(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
            NULLIFY(domain)
            CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
              DO dimensionIdx=1,numberOfDimensions
                !Default to version 1 of each node derivative
                localDofIdx=geometricVariable%components(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & nodes(nodeIdx)%derivatives(1)%versions(1)
                x(dimensionIdx)=geometricParameters(localDofIdx)
              ENDDO !dimensionIdx
              !Loop over the derivatives
              DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                globalDerivativeIndex=domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX
                CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
                  & globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,VALUE,err,error,*999)
                !Default to version 1 of each node derivative
                localDofIdx=dependentVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & nodes(nodeIdx)%derivatives(derivativeIdx)%versions(1)                
                CALL Field_ParameterSetUpdateLocalDof(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,localDofIdx, &
                  & value,err,error,*999)
                IF(MOD(variableType,FIELD_NUMBER_OF_VARIABLE_SUBTYPES)==FIELD_U_VARIABLE_TYPE) THEN
                  IF(domainNodes%nodes(nodeIdx)%BOUNDARY_NODE) THEN
                    !If we are a boundary node then set the analytic value on the boundary
                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType,localDofIdx, &
                      & BOUNDARY_CONDITION_FIXED,value,err,error,*999)
                  ELSE
                    CALL Field_ParameterSetUpdateLocalDof(dependentField,variableType,FIELD_VALUES_SET_TYPE,localDofIdx, &
                      & VALUE,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ELSE
            CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
          ENDIF
        ENDDO !component_idx
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
    ELSE
      !for single physics diffusion problems use standard analytic calculate           
      DO variableIdx=1,dependentField%NUMBER_OF_VARIABLES
        variableType=dependentField%variables(variableIdx)%VARIABLE_TYPE
        NULLIFY(dependentVariable)
        CALL Field_VariableGet(dependentField,variableType,dependentVariable,err,error,*999)
        CALL Field_ParameterSetEnsureCreated(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        DO componentIdx=1,dependentVariable%NUMBER_OF_COMPONENTS
          IF(dependentVariable%components(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
            NULLIFY(domain)
            CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
            NULLIFY(domainTopology)
            CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
            NULLIFY(domainNodes)
            CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
            !Loop over the local nodes excluding the ghosts.
            DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
              DO dimensionIdx=1,numberOfDimensions
                !Default to version 1 of each node derivative
                localDofIdx=geometricVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                  & nodes(nodeIdx)%derivatives(1)%versions(1)
                x(dimensionIdx)=geometricParameters(localDofIdx)
              ENDDO !dimensionIdx
              !Loop over the derivatives
              DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                globalDerivativeIndex=domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX
                CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,0.0_DP, &
                  & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,initialValue, &
                  & err,error,*999)
                CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time, &
                  & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,value, &
                  & err,error,*999)
                DO versionIdx=1,domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%numberOfVersions
                  localDofIdx=dependentVariable%components(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                    & nodes(nodeIdx)%derivatives(derivativeIdx)%versions(versionIdx)
                  CALL Field_ParameterSetUpdateLocalDof(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & localDofIdx,value,err,error,*999)
                  IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                    IF(domainNodes%nodes(nodeIdx)%BOUNDARY_NODE) THEN
                      !If we are a boundary node then set the analytic value on the boundary
                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(boundaryConditions,dependentField,variableType, &
                        & localDofIdx,BOUNDARY_CONDITION_FIXED,initialValue,err,error,*999)
                    ELSE
                      !Set the initial condition.
                      CALL Field_ParameterSetUpdateLocalDof(dependentField,variableType,FIELD_VALUES_SET_TYPE, &
                        & localDofIdx,initialValue,err,error,*999)
                    ENDIF
                  ENDIF
                ENDDO !versionIdx
              ENDDO !derivIdx
            ENDDO !nodeIdx
          ELSE
            CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
          ENDIF
        ENDDO !component_idx
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateStart(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL Field_ParameterSetUpdateFinish(dependentField,variableType,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
    ENDIF
    IF(ASSOCIATED(materialsField)) &
      & CALL Field_ParameterSetDataRestore(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & materialsParameters,err,error,*999)
    IF(ASSOCIATED(analyticField)) &
      & CALL Field_ParameterSetDataRestore(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & analyticParameters,err,error,*999)            
    CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
      & geometricParameters,err,error,*999)            

    EXITS("Diffusion_BoundaryConditionAnalyticCalculate")
    RETURN
999 ERRORS("Diffusion_BoundaryConditionAnalyticCalculate",err,error)
    EXITS("Diffusion_BoundaryConditionAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Diffusion_BoundaryConditionAnalyticCalculate

  !
  !================================================================================================================================
  !
  
  !>Evaluate the analytic solutions for a diffusion equation
  SUBROUTINE Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
    & globalDerivativeIndex,componentNumber,analyticParameters,materialsParameters,VALUE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<The equations set to evaluate
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: x(:) !<X(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<TANGENTS(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<NORMAL(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivativeIndex !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: value !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: equationsSubType
    REAL(DP) :: k,phi,A,B,C,D,A1,A2,A3,A4
    REAL(DP) :: aParam,bParam,cParam,kParam,lParam,constParam,betaParam,lambdaParam,muParam
    TYPE(VARYING_STRING) :: localError

    !These are parameters for the analytical solution
!     k = 1.0_DP !this is a time decay constant for the exponential term
!     phi = 0.785398163397_DP !pi/4 - this sets the orientation of the solution relative to the axes
    k = 10.0_DP !this is a time decay constant for the exponential term
    phi = 1.0_DP !pi/4 - this sets the orientation of the solution relative to the axes
 
    !Solution parameters for 
    A1 = 0.4_DP
    A2 = 0.3_DP
    A3 = 0.2_DP
    A4 = 0.1_DP

    ENTERS("Diffusion_AnalyticFunctionsEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.ALLOCATED(equationsSet%specification)) &
      & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
    IF(SIZE(equationsSet%specification,1)<3) &
      & CALL FlagError("Equations set specification must have at least three entries for a diffusion type equations set.", &
      & err,error,*999)
    
    equationsSubType=equationsSet%specification(3)
    SELECT CASE(equationsSubType)
    CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
        !For \del u/\del t = a \del^2 u/\del x^2
        !u(x,t)=A.exp(-4.\mu^2.t)cos(\mu.x+B)+C
        !see http://eqworld.ipmnet.ru/en/solutions/lpde/lpde101.pdf
        !OpenCMISS has \del u/\del t + k \del^2 u/\del x^2 = 0, thereform with \mu=2.\pi/L we have
        !u(x,t)=A.exp(4.\pi^2.k.t/L^2)cos(2.\pi.x/L+B)+C
        kParam=materialsParameters(1)
        aParam=analyticParameters(1)
        bParam=analyticParameters(2)
        cParam=analyticParameters(3)
        lParam=analyticParameters(4)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=aParam*EXP(4.0_DP*PI**2*kParam*time/lParam**2)*COS(2.0_DP*PI*x(1)/lParam+bParam)+cParam
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=0.0_DP
          CASE(GLOBAL_DERIV_S1)
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
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
        !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=EXP(-k*time)*SIN((SQRT(k))*(x(1)*COS(phi)+x(2)*SIN(phi)))!Need to specify time, k and phi!
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
        !u=A1*exp(-t)*(x^2+y^2+z^2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The specified analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
          & " is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
        !u=exp(-kt)*sin(sqrt(k)*(x*cos(phi)+y*sin(phi)))
        !These are parameters for the 3D analytical solution with a linear source
        A = -0.25_DP
        B = 0.5_DP   
        C = 0.5_DP
        D = 0.5_DP
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=EXP(A*time)*EXP(B*x(1))*EXP(C*x(2))*EXP(D*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
     SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
        !For del u/del t = del^2 u/del x^2 + a + bu + cu^m with a = 0 then
        !u(x,t) = [+/-\beta + C.exp(\lamba.t+/\mu.x)]^(2/1-m) where
        !\beta=\sqrt(-c/b); \lamba=b(1-m)(m+3)/(2(m+1)); \mu = \sqrt((b(1-m)^2)/(2.(m+1))
        !see http://eqworld.ipmnet.ru/en/solutions/npde/npde1104.pdf
        aParam=materialsParameters(1)
        bParam=materialsParameters(2)
        cParam=materialsParameters(3)
        betaParam=SQRT(-cParam/bParam)
        lambdaParam=-5.0_DP*bParam/6.0_DP
        muParam=SQRT(bParam/6.0_DP)
        constParam=1.0_DP
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=1.0_DP/(betaParam+constParam*EXP(lambdaParam*time+muParam*x(1)))**2
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            value=0.0_DP                                 
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
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
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
        !For del u/del t = del^2 u/del x^2 + a + b.e^c.u
        !u(x,t) = -2/c.ln[+/-\beta + C.exp(+/-\mu.x-a.c.t/2)] where
        !\beta=\sqrt(-b/a); \mu = \sqrt(a.c/2)
        !see http://eqworld.ipmnet.ru/en/solutions/npde/npde1105.pdf
        aParam=materialsParameters(1)
        bParam=materialsParameters(2)
        cParam=materialsParameters(3)
        constParam=1.0_DP
        betaParam=SQRT(-bParam/aParam)
        muParam=SQRT(aParam*cParam/2.0_DP)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=-2.0_DP/cParam*LOG(betaParam+constParam*EXP(muParam*X(1)-aParam*cParam*time/2.0_DP))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP 
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
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
    CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
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
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
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
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
        CASE(FIELD_V_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELVDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
        CASE(FIELD_U1_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A3*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELU1DELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
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
        CASE(FIELD_U2_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A4*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S2)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(GLOBAL_DERIV_S1_S2)
            CALL FlagError("Not implmented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivativeIndex,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELU2DELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivativeIndex)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding                                    
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
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(equationsSubType,"*",err,error))// &
        & " is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT
     
    EXITS("Diffusion_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("Diffusion_AnalyticFunctionsEvaluate",err,error)
    RETURN 1
  END SUBROUTINE Diffusion_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equation type of a classical field equations set class.
  SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a diffusion equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
        !Need to define the functions diffusion_equation_equations_set_linear_source_setup etc
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        CALL Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetNonlinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
        CALL Diffusion_EquationsSetNonlinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL Diffusion_EquationsSetNonlinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
         CALL FlagError("Not implemented.",err,error,*999)
!        CALL Diffusion_EquationsSetNonlinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a diffusion equation type of a classical field equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE Diffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Diffusion_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,& 
         & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE,&
         & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%SOLUTION_METHOD=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
          localError="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT    
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a diffusion equation type of an classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Diffusion_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Diffusion_EquationsSetSolutionMethodSet",err,error)
    EXITS("Diffusion_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Diffusion_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a diffusion equation type of a classical field equations set class.
  SUBROUTINE Diffusion_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Diffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
        & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a diffusion type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Diffusion_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Diffusion_EquationsSetSpecificationSet",err,error)
    EXITS("Diffusion_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Diffusion_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the linear diffusion equation.
  SUBROUTINE Diffusion_EquationsSetLinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS, NUMBER_OF_MATERIALS_COMPONENTS, NUMBER_OF_SOURCE_COMPONENTS,imy_matrix,Ncompartments, &
      & GEOMETRIC_COMPONENT_NUMBER
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SOURCE
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD,EQUATIONS_SET_FIELD_FIELD
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: num_var,num_var_count,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS    
    INTEGER(INTG) :: EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS    
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
    INTEGER(INTG), ALLOCATABLE :: VARIABLE_TYPES(:),VARIABLE_U_TYPES(:),COUPLING_MATRIX_STORAGE_TYPE(:), &
      & COUPLING_MATRIX_STRUCTURE_TYPE(:)

    ENTERS("Diffusion_EquationsSetLinearSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Diffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES = 1
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                !Create the auto created equations set field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                  & EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,"Equations Set Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_GENERAL_TYPE,&
                  & err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                  & EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & [FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Equations",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,EQUATIONS_SET_FIELD_NUMBER_OF_VARIABLES, &
                  & err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                  & EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
              ENDIF
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
              IF(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1, 1_INTG, ERR, ERROR, *999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 2, 1_INTG, ERR, ERROR, *999)
              ENDIF
            ENDIF
!!TODO: Check valid setup
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
            !do nothing 
          CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS = 2
              EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
              IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,&
                  & GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,& 
                  & EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)                
                DO component_idx = 1, EQUATIONS_SET_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                END DO
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD,GEOMETRIC_SCALING_TYPE, &
                  & err,error,*999)
              ELSE
                !Do nothing
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_DISPLACEMENT_SET_TYPE, ERR, ERROR, *999)
                CALL FIELD_PARAMETER_SET_CREATE(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, FIELD_U_VARIABLE_TYPE, &
                  & FIELD_MESH_VELOCITY_SET_TYPE, ERR, ERROR, *999)
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              !Do nothing
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a linear diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          END SELECT
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
            CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE, &
              & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
                !Create the auto created dependent field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & DEPENDENT_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "U",err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & "del U/del n",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                !Default to the geometric interpolation setup
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  !Default the scaling to the geometric field scaling
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                  localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                  !Check the field created by advection-diffusion routines for the coupled problem
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, & 
                    & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                    & err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    component_idx=1
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, & 
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT

                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                  !uses number of compartments to check that appropriate number and type of variables have been set on the
                  !dependent field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)                 
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2*Ncompartments,err,error,*999)
                  !Create & populate array storing all of the relevant variable types against which to check the field variables
                  ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                  DO num_var=1,Ncompartments
                    VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                    VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  ENDDO
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES,err,error,*999)

                  DO num_var=1,2*Ncompartments
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var), & 
                      & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                    CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),FIELD_DP_TYPE,err,error,*999)
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),1, &
                      & err,error,*999)
                  ENDDO
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    component_idx=1
                    DO num_var=1,2*Ncompartments
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,VARIABLE_TYPES(num_var),component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    ENDDO
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
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
                  & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                  CALL FlagError("Not implemented.",err,error,*999)
                CASE DEFAULT
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, & 
                    & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                    & err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
                  CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                    DO component_idx=1,NUMBER_OF_DIMENSIONS
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                        & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                        & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    ENDDO !component_idx
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
                    localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                      & " is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                END SELECT
              ENDIF
            END SELECT
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                    & err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "Materials",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                  ELSEIF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    !Linear source. Materials field components are 1 for each dimension and 1 for the linear source
                    !i.e., k and a in div(k.grad(u(x)))=a(x)u(x)+c(x)
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                  ELSE
                    NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2
                  ENDIF
                  !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                  !Default the k materials components to the geometric interpolation setup with constant interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                  !Default the source materials components to the first component geometric interpolation with constant
                  !interpolation
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                      CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                      CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    ENDDO !components_idx
                  ENDIF
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                  CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                    & MATERIALS_FIELD,err,error,*999)
                  CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                  CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                  CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                    & err,error,*999)
                  CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                    & GEOMETRIC_FIELD,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & "Materials",err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                  CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS
                  !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,err,error,*999)
                  !Default the k materials components to the geometric interpolation setup with constant interpolation
                  DO component_idx=1,NUMBER_OF_DIMENSIONS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO !component_idx
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  DO component_idx=1,NUMBER_OF_MATERIALS_COUPLING_COMPONENTS
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  ENDDO
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)

                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                    & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                      & err,error,*999)
                  ELSE
                    CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS+1, &
                      & err,error,*999)
                  ENDIF
                CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                    & err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                    & err,error,*999)
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  NUMBER_OF_MATERIALS_COUPLING_COMPONENTS=Ncompartments
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COUPLING_COMPONENTS,err,error,*999)
                END SELECT
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS             
                ELSE IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                  !Constant source. Materials field components are 1 for each dimension and 1 for the constant source
                  !i.e., k and c in div(k.grad(u(x)))=c(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+1
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                  EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                  CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                  Ncompartments=EQUATIONS_SET_FIELD_DATA(2)
                  DO component_idx=1,Ncompartments
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,0.0_DP,err,error,*999)
                  ENDDO !component_idx
                ENDIF
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                  !Now set the linear source values to 1.0
                  DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                  ENDDO !component_idx
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Create the auto created source field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SOURCE% &
                  & SOURCE_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,"Source Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Source",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                NUMBER_OF_SOURCE_COMPONENTS=1
                !Set the number of source components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOURCE_COMPONENTS,err,error,*999)
                !Default the source components to the geometric interpolation setup with constant interpolation
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                  DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                ENDIF
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SOURCE%SOURCE_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_SOURCE=>EQUATIONS_SET%SOURCE
            IF(ASSOCIATED(EQUATIONS_SOURCE)) THEN
              IF(EQUATIONS_SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN
                !Finish creating the source field
                CALL FIELD_CREATE_FINISH(EQUATIONS_SOURCE%SOURCE_FIELD,err,error,*999)
                !Set the default values for the source field
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                  NUMBER_OF_SOURCE_COMPONENTS=1
                ELSE
                  NUMBER_OF_SOURCE_COMPONENTS=0
                ENDIF
                !Now set the source values to 1.0
                DO component_idx=1,NUMBER_OF_SOURCE_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_SOURCE%SOURCE_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set source is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c  T y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                  IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                    EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                    IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                      IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                        CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                          & err,error,*999)
                        EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=0.0_DP
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a no source diffusion equation requires that there be 1 geometric dimension."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=4
                            !Set analytic function type
                            EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a no source diffusion equation requires that there be 2 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a no source diffusion equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a constant source diffusion equation requires that there be 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a constant source diffusion equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a linear source diffusion equation requires that there be 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a linear source diffusion equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
                            !Check that domain is 2D
                            IF(NUMBER_OF_DIMENSIONS/=2) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a multi-compartment diffusion equation requires that there be 2 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM
                          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM)
                            !Check that domain is 3D
                            IF(NUMBER_OF_DIMENSIONS/=3) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " for a multi-compartment diffusion requires that there be 3 geometric dimensions."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Set number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= & 
                              & EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a multi-compartment diffusion equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE DEFAULT
                          localError="The equation set subtype of "// &
                            & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is invalid for an analytical diffusion equation."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                        !Create analytic field if required
                        IF(NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                          IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                            !Create the auto created source field
                            CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                              & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                            CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",err,error,*999)
                            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                            CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                              & err,error,*999)
                            CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                              & err,error,*999)
                            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                              & GEOMETRIC_DECOMPOSITION,err,error,*999)
                            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                              & GEOMETRIC_FIELD,err,error,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,err,error,*999)
                            CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                              & err,error,*999)
                            CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & "Analytic",err,error,*999)
                            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_DP_TYPE,err,error,*999)
                            !Set the number of analytic components
                            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                            CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                            DO component_idx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            ENDDO !component_idx
                            !Default the field scaling to that of the geometric field
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                          ELSE
                            !Check the user specified field
                            CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                            CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                            CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                            IF(NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                            ELSE
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                            ENDIF
                            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                              & err,error,*999)
                            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                          ENDIF
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set materials has not been finished.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set materials is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  !Finish creating the analytic field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                  !Set the default values for the analytic field
                  CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
                    SELECT CASE(EQUATIONS_ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
                      !Set A
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                      !Set B
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
                      !Set C
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,3,1.0_DP,err,error,*999)
                      !Set L
                      CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                        & FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)                      
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                      !Do nothing
                    CASE DEFAULT
                      localError="The specified analytic function type of "// &
                        & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                        & " is invalid for a no source diffusion equation."
                      CALL FlagError(localError,err,error,*999)
                     END SELECT
                  CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE DEFAULT
                    localError="The equation set subtype of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                      & " is invalid for an analytical linear diffusion equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))      
              CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_V_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
                CALL FlagError("Not implemented.",err,error,*999)
              CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
                EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)
                imy_matrix = EQUATIONS_SET_FIELD_DATA(1)
                Ncompartments = EQUATIONS_SET_FIELD_DATA(2)    
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,Ncompartments-1,err,error,*999)

                ALLOCATE(VARIABLE_TYPES(2*Ncompartments))
                ALLOCATE(VARIABLE_U_TYPES(Ncompartments-1))
                DO num_var=1,Ncompartments
                  VARIABLE_TYPES(2*num_var-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                  VARIABLE_TYPES(2*num_var)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(num_var-1))
                ENDDO
                num_var_count=0
                DO num_var=1,Ncompartments
                  IF(num_var/=imy_matrix)THEN
                    num_var_count=num_var_count+1
                    VARIABLE_U_TYPES(num_var_count)=VARIABLE_TYPES(2*num_var-1)
                  ENDIF
                ENDDO
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix-1),err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,VARIABLE_U_TYPES,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,VARIABLE_TYPES(2*imy_matrix),err,error,*999)
                CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CASE DEFAULT
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              END SELECT
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN              
                CALL EquationsMapping_SourceVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              ENDIF
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                  [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
              ELSE
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)THEN
                    ALLOCATE(COUPLING_MATRIX_STORAGE_TYPE(Ncompartments-1))
                    ALLOCATE(COUPLING_MATRIX_STRUCTURE_TYPE(Ncompartments-1))
                    DO num_var=1,Ncompartments-1
                      COUPLING_MATRIX_STORAGE_TYPE(num_var)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                      COUPLING_MATRIX_STRUCTURE_TYPE(num_var)=EQUATIONS_MATRIX_FEM_STRUCTURE
                    ENDDO
                    CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                      & COUPLING_MATRIX_STORAGE_TYPE, &
                      & err,error,*999)      
                    CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                      COUPLING_MATRIX_STRUCTURE_TYPE,err,error,*999)
                  ENDIF
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
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
              localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a linear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not a linear diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Diffusion_EquationsSetLinearSetup")
    RETURN
999 ERRORSEXITS("Diffusion_EquationsSetLinearSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_EquationsSetLinearSetup

  !
  !================================================================================================================================
  !

  !>Sets up the non-linear diffusion equation.
  SUBROUTINE Diffusion_EquationsSetNonlinearSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DIMENSIONS,NUMBER_OF_MATERIALS_COMPONENTS
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_EquationsSetNonlinearSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE.OR. &
        & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Diffusion_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        !-----------------------------------------------------------------
        ! D e p e n d e n t   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)            
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,2,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & "U",err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & "del U/del n",err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, & 
                & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE, & 
                & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                & err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & NUMBER_OF_DIMENSIONS,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, & 
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !component_idx
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! M a t e r i a l s   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field                
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,err,error,*999)
                CALL FIELD_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,"Materials Field",err,error,*999)
                CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                  & GEOMETRIC_FIELD,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & "Materials",err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ENDIF
                !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                !Default the k materials components to the geometric interpolation setup with constant interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDDO !component_idx
                !Default the source materials components to the first component geometric interpolation with constant
                !interpolation                
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                ENDDO !components_idx
                !Default the field scaling to that of the geometric field
                CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              ELSE
                !Check the user specified field
                CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+2                                
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_DIMENSIONS,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                  !Quadratic source. Materials field components are 1 for each dimension and 3 for the quadratic source
                  !i.e., k and a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)u(x)+c(x)u^2(x)
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ELSE
                  !Exponential source. Matierals field components are 1 for each dimension and 3 for the exponential source
                  !i.e., k, a, b and c in del u/del t = div(k.grad(u(x)))+a(x)+b(x)e^[c(x)u(x)]
                  NUMBER_OF_MATERIALS_COMPONENTS=NUMBER_OF_DIMENSIONS+3
                ENDIF
                !First set the k values to 1.0
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
                !Set the source values to 1.0
                DO component_idx=NUMBER_OF_DIMENSIONS+1,NUMBER_OF_MATERIALS_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,err,error,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing put the constant source directly into the RHS
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing put the constant source directly into the RHS
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! A n a l y t i c   t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                  EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                  IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                    IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                      GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                      IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                        CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS, &
                          & err,error,*999)
                        IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 3,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Check that the a parameter is zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,A_PARAM,err,error,*999)
                            IF(ABS(A_PARAM)>ZERO_TOLERANCE)  &
                              & CALL FlagError("The 1st material component must be zero.",err,error,*999)
                            !Check that the b parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,err,error,*999)
                            IF(B_PARAM<ZERO_TOLERANCE)  &
                              & CALL FlagError("The 2nd material component must be greater than zero.",err,error,*999)
                            !Check to ensure we get real solutions
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,err,error,*999)
                            IF((B_PARAM*C_PARAM)>ZERO_TOLERANCE) &
                              & CALL FlagError("The product of the 2nd and 3rd material components must not be positive.", &
                              & err,error,*999)
                            !Set the number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=1
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a nonlinear diffusion equation with a quadratic source."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        ELSE
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_DIMENSIONS/=1) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NumberToVString(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 3,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Check that the a parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,1,A_PARAM,err,error,*999)
                            IF(ABS(A_PARAM)<ZERO_TOLERANCE)  &
                              & CALL FlagError("The 1st material component must not be zero.",err,error,*999)
                            !Check that the c parameter is not zero.
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,3,C_PARAM,err,error,*999)
                            IF(ABS(C_PARAM)<ZERO_TOLERANCE)  &
                              & CALL FlagError("The 3rd material component must not be zero.",err,error,*999)
                            !Check to ensure we get real solutions
                            CALL FIELD_PARAMETER_SET_GET_CONSTANT(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VALUES_SET_TYPE,2,B_PARAM,err,error,*999)
                            IF((A_PARAM*B_PARAM)>ZERO_TOLERANCE) &
                              & CALL FlagError("The product of the 1st and 2nd material components must not be positive.", &
                              & err,error,*999)
                            IF((A_PARAM*C_PARAM)<ZERO_TOLERANCE) &
                              & CALL FlagError("The product of the 1st and 3rd material components must not be negative.", &
                              & err,error,*999)
                            !Set the number of analytic field components
                            NUMBER_OF_ANALYTIC_COMPONENTS=0
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE= &
                              & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a nonlinear diffusion equation with an exponential source."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        ENDIF
                        !Create analytic field if required
                        IF(NUMBER_OF_ANALYTIC_COMPONENTS>=1) THEN
                          IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                            !Create the auto created source field
                            CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION, &
                              & EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                            CALL FIELD_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,"Analytic Field",err,error,*999)
                            CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                            CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_INDEPENDENT_TYPE, &
                              & err,error,*999)
                            CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION, &
                              & err,error,*999)
                            CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD, &
                              & GEOMETRIC_DECOMPOSITION,err,error,*999)
                            CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,EQUATIONS_SET%GEOMETRY% &
                              & GEOMETRIC_FIELD,err,error,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,1,err,error,*999)
                            CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,[FIELD_U_VARIABLE_TYPE], &
                              & err,error,*999)
                            CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & "Analytic",err,error,*999)
                            CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                            CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_DP_TYPE,err,error,*999)
                            !Set the number of analytic components
                            CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
                            CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE,1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                            DO component_idx=1,NUMBER_OF_ANALYTIC_COMPONENTS
                              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                                & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            ENDDO !component_idx
                            !Default the field scaling to that of the geometric field
                            CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE, &
                              & err,error,*999)
                            CALL FIELD_SCALING_TYPE_SET(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                          ELSE
                            !Check the user specified field
                            CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                            CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                            CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                            CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                            IF(NUMBER_OF_ANALYTIC_COMPONENTS==1) THEN
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                            ELSE
                              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                            ENDIF
                            CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
                              & err,error,*999)
                            CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                              & NUMBER_OF_ANALYTIC_COMPONENTS,err,error,*999)
                          ENDIF
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set materials is not finished.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations set materials is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations analytic is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_ANALYTIC=>EQUATIONS_SET%ANALYTIC
            IF(ASSOCIATED(EQUATIONS_ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  !Finish creating the analytic field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,err,error,*999)
                  !Set the default values for the analytic field
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
                    !Do nothing
                  CASE DEFAULT
                    localError="The equation set subtype of "// &
                      & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                      & " is invalid for an analytical nonlinear diffusion equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_NONLINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              ! Use the analytic Jacobian calculation
              CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                & err,error,*999)
              !Set up matrix storage and structure
              IF(EQUATIONS%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                    & err,error,*999)
                 CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ELSE
                SELECT CASE(EQUATIONS%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
              ENDIF
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
              localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a nonlinear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not a nonlinear diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Diffusion_EquationsSetNonlinearSetup")
    RETURN
999 ERRORS("Diffusion_EquationsSetNonlinearSetup",err,error)
    EXITS("Diffusion_EquationsSetNonlinearSetup")
    RETURN 1
    
  END SUBROUTINE Diffusion_EquationsSetNonlinearSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the diffusion problem.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a diffusion equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DIFFUSION_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DIFFUSION_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Performs pre-solve operations for the a diffusion problem.
  SUBROUTINE Diffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    LOGICAL :: updateMaterials
    LOGICAL :: updateBoundaryConditions
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem    
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolve",err,error,*999)

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for a diffusion problem.",err,error,*999)
    
    updateMaterials = .FALSE.    
    updateBoundaryConditions = .TRUE.

    IF(updateMaterials) THEN
      !CALL DIFFUSION_EQUATION_PRE_SOLVE_UPDATE_MATERIALS_FIELD(controlLoop,solver,err,error,*999)
    ENDIF
    
    !IF(updateBoundaryConditions) THEN
    !  CALL Diffusion_PreSolveUpdateBoundaryConditions(controlLoop,solver,err,error,*999)
    !ENDIF

    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
      ! do nothing ???
      CALL Diffusion_PreSolveUpdateAnalyticValues(solver,err,error,*999)
    CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"ALE diffusion pre solve... ",err,error,*999)
      IF(solver%DYNAMIC_SOLVER%ALE) THEN
        !First update mesh and calculate boundary velocity values
        CALL Diffusion_PreSolveALEUpdateMesh(controlLoop,solver,err,error,*999)
        !Then apply both normal and moving mesh boundary conditions
        !CALL DIFFUSION_PRE_SOLVE_UPDATE_BOUNDARY_CONDITIONS(controlLoop,solver,err,error,*999)
      ELSE  
        CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
      END IF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Diffusion_PreSolve")
    RETURN
999 ERRORSEXITS("Diffusion_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolve
      
  !   
  !================================================================================================================================
  !
  !>Within the diffusion pre-solve, update the boundary conditions
  SUBROUTINE Diffusion_PreSolveUpdateBoundaryConditions(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
!     TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD
! !    TYPE(FIELD_TYPE), POINTER :: FIELD !<A pointer to the field
!     TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
!     TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
!     TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
!     TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
!     TYPE(EquationsType), POINTER :: EQUATIONS
!     TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
!     TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
! !    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: DOMAIN_TOPOLOGY
!     TYPE(VARYING_STRING) :: localError
! !    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
! !    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
! !    REAL(DP), POINTER :: BOUNDARY_VALUES(:)
!     REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
!     INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE
! 
!     REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
!     REAL(DP) :: VALUE,X(3) !<The value to add
! !     REAL(DP) :: k_xx, k_yy, k_zz
!     INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,variable_idx
!     INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
!     INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
!     INTEGER(INTG) :: GLOBAL_DERIV_INDEX
! !    INTEGER(INTG) :: FIELD_SET_TYPE !<The field parameter set identifier \see FIELD_ROUTINES_ParameterSetTypes,FIELD_ROUTINES
! !    INTEGER(INTG) :: DERIVATIVE_NUMBER !<The node derivative number
! !    INTEGER(INTG) :: COMPONENT_NUMBER !<The field variable component number
! !    INTEGER(INTG) :: TOTAL_NUMBER_OF_NODES !<The total number of (geometry) nodes
! !    INTEGER(INTG) :: LOCAL_NODE_NUMBER
! !    INTEGER(INTG) :: EQUATIONS_SET_IDX
! !    INTEGER(INTG) :: equations_row_number

    ENTERS("Diffusion_PreSolveUpdateBoundaryConditions",err,error,*999)

    CALL FlagError("Not implemented.",err,error,*999)!This routine previously set analytic BCs, but this has been moved. Needs rewriting to set
    !boundary conditions from file, time varying if appropriate.

!     IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
!        !write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!        !write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
!       IF(ASSOCIATED(SOLVER)) THEN
!         IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
!           SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
!             CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
!                 SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
!                 IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
!                   SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
!                   EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                   IF(ASSOCIATED(EQUATIONS)) THEN
!                     EQUATIONS_SET=>equations%equationsSet
!                     IF(ASSOCIATED(EQUATIONS_SET)) THEN
!                       IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                         DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
!                         IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
!                           GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
!                           IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
!                             CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
!                               & NUMBER_OF_DIMENSIONS,err,error,*999)
!                             NULLIFY(GEOMETRIC_VARIABLE)
!                             NULLIFY(GEOMETRIC_PARAMETERS)
!                             CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
!                             CALL Field_ParameterSetDataGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
!                               & GEOMETRIC_PARAMETERS,err,error,*999)
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
!                               variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
!                               FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
!                               IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                                 DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
!                                     & FIELD_NODE_BASED_INTERPOLATION) THEN
!                                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                                     IF(ASSOCIATED(DOMAIN)) THEN
!                                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                           !Loop over the local nodes excluding the ghosts.
!                                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                             !!TODO \todo We should interpolate the geometric field here and the node position.
!                                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
!                                               local_ny= & 
!                                           & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
!                                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
!                                             ENDDO !dim_idx
!                                             !Loop over the derivatives
!                                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
!                                               ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
!                                               GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%globalDerivativeIndex_INDEX(deriv_idx)
!                                               CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
!                                                 & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
!                                                 & ANALYTIC_FUNCTION_TYPE,err,error,*999)
!                                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
!                                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
!                                               CALL Field_ParameterSetUpdateLocalDof(DEPENDENT_FIELD,variable_type, &
!                                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
!                                               BOUNDARY_CONDITION_CHECK_VARIABLE=SOLVER_EQUATIONS%BOUNDARY_CONDITIONS% &
!                                                 & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr% & 
!                                                 & CONDITION_TYPES(local_ny)
!                                               IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
!                                                CALL Field_ParameterSetUpdateLocalDof(DEPENDENT_FIELD, & 
!                                                  & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
!                                                  & VALUE,err,error,*999)
!                                               ENDIF
! !                                              IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
! !                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
!                                                   !If we are a boundary node then set the analytic value on the boundary
! !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
! !                                                ENDIF
! !                                              ENDIF
!                                             ENDDO !deriv_idx
!                                           ENDDO !node_idx
!                                         ELSE
!                                           CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
!                                         ENDIF
!                                       ELSE
!                                         CALL FlagError("Domain topology is not associated.",err,error,*999)
!                                       ENDIF
!                                     ELSE
!                                       CALL FlagError("Domain is not associated.",err,error,*999)
!                                     ENDIF
!                                   ELSE
!                                     CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
!                                   ENDIF
!                                 ENDDO !component_idx
!                                 CALL CALL Field_ParameterSetUpdateStart(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
!                                 CALL Field_ParameterSetUpdateFinish(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
!                                 CALL CALL Field_ParameterSetUpdateStart(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
!                                 CALL Field_ParameterSetUpdateFinish(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
!                               ELSE
!                                 CALL FlagError("Field variable is not associated.",err,error,*999)
!                               ENDIF
! 
!                              ENDDO !variable_idx
!                              CALL Field_ParameterSetDataRestore(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
!                               & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
!                           ELSE
!                             CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
!                           ENDIF            
!                         ELSE
!                           CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
!                         ENDIF
!                       ELSE
!                         !CALL FlagError("Equations set analytic is not associated.",err,error,*999)
!                       ENDIF
!                     ELSE
!                       CALL FlagError("Equations set is not associated.",err,error,*999)
!                     ENDIF
!                   ELSE
!                     CALL FlagError("Equations are not associated.",err,error,*999)
!                   END IF                
!                 ELSE
!                   CALL FlagError("Solver equations are not associated.",err,error,*999)
!                 END IF  
!                 CALL CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,err,error,*999)
!                 CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,err,error,*999)
!             !do nothing?! 
!             CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
!             !do nothing?! 
!             CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
!                 SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
!                 IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
!                   SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
!                   EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(1)%EQUATIONS
!                   IF(ASSOCIATED(EQUATIONS)) THEN
!                     EQUATIONS_SET=>equations%equationsSet
!                     IF(ASSOCIATED(EQUATIONS_SET)) THEN
!                      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
!                         DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
!                         IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
!                           GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
!                           IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
!                             CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,&
!                               & NUMBER_OF_DIMENSIONS,err,error,*999)
!                             NULLIFY(GEOMETRIC_VARIABLE)
!                             NULLIFY(GEOMETRIC_PARAMETERS)
!                             CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
!                             CALL Field_ParameterSetDataGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,& 
!                               & GEOMETRIC_PARAMETERS,err,error,*999)
!                              DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
!                               variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
!                               FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%ptr
!                               IF(ASSOCIATED(FIELD_VARIABLE)) THEN
!                                 DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
!                                   IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== & 
!                                     & FIELD_NODE_BASED_INTERPOLATION) THEN
!                                     DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
!                                     IF(ASSOCIATED(DOMAIN)) THEN
!                                       IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
!                                         DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
!                                         IF(ASSOCIATED(DOMAIN_NODES)) THEN
!                                           !Loop over the local nodes excluding the ghosts.
!                                           DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!                                             !!TODO \todo We should interpolate the geometric field here and the node position.
!                                             DO dim_idx=1,NUMBER_OF_DIMENSIONS
!                                               local_ny= & 
!                                           & GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP(1,node_idx)
!                                               X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
!                                             ENDDO !dim_idx
!                                             !Loop over the derivatives
!                                             DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
!                                               ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
!                                               GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%globalDerivativeIndex_INDEX(deriv_idx)
!                                               CALL DIFFUSION_EQUATION_ANALYTIC_FUNCTIONS(VALUE,X, & 
!                                                 & CURRENT_TIME,variable_type,GLOBAL_DERIV_INDEX, &
!                                                 & ANALYTIC_FUNCTION_TYPE,err,error,*999)
!                                               local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
!                                                 & NODE_PARAM2DOF_MAP(deriv_idx,node_idx)
!                                               CALL Field_ParameterSetUpdateLocalDof(DEPENDENT_FIELD,variable_type, &
!                                                 & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
!                                               BOUNDARY_CONDITION_CHECK_VARIABLE=SOLVER_EQUATIONS%BOUNDARY_CONDITIONS% &
!                                                 & BOUNDARY_CONDITIONS_VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr% & 
!                                                 & CONDITION_TYPES(local_ny)
!                                               IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
!                                                CALL Field_ParameterSetUpdateLocalDof(DEPENDENT_FIELD, & 
!                                                  & variable_type,FIELD_VALUES_SET_TYPE,local_ny, & 
!                                                  & VALUE,err,error,*999)
!                                               ENDIF
! !                                              IF(variable_type==FIELD_U_VARIABLE_TYPE .OR. variable_type==FIELD_V_VARIABLE_TYPE) THEN
! !                                                IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
! !                                                   If we are a boundary node then set the analytic value on the boundary
! !                                                  CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
! !                                                    & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
! !                                                ENDIF
! !                                              ENDIF
!                                             ENDDO !deriv_idx
!                                           ENDDO !node_idx
!                                         ELSE
!                                           CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
!                                         ENDIF
!                                       ELSE
!                                         CALL FlagError("Domain topology is not associated.",err,error,*999)
!                                       ENDIF
!                                     ELSE
!                                       CALL FlagError("Domain is not associated.",err,error,*999)
!                                     ENDIF
!                                   ELSE
!                                     CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
!                                   ENDIF
!                                 ENDDO !component_idx
!                                 CALL CALL Field_ParameterSetUpdateStart(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
!                                 CALL Field_ParameterSetUpdateFinish(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
!                                 CALL CALL Field_ParameterSetUpdateStart(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
!                                 CALL Field_ParameterSetUpdateFinish(DEPENDENT_FIELD,variable_type, &
!                                  & FIELD_VALUES_SET_TYPE,err,error,*999)
!                               ELSE
!                                 CALL FlagError("Field variable is not associated.",err,error,*999)
!                               ENDIF
! 
!                              ENDDO !variable_idx
!                              CALL Field_ParameterSetDataRestore(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,& 
!                               & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
!                           ELSE
!                             CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
!                           ENDIF            
!                         ELSE
!                           CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
!                         ENDIF
!                       ELSE
!                         !CALL FlagError("Equations set analytic is not associated.",err,error,*999)
!                       ENDIF
!                     ELSE
!                       CALL FlagError("Equations set is not associated.",err,error,*999)
!                     ENDIF
!                   ELSE
!                     CALL FlagError("Equations are not associated.",err,error,*999)
!                   END IF                
!                 ELSE
!                   CALL FlagError("Solver equations are not associated.",err,error,*999)
!                 END IF  
!                 CALL CALL Field_ParameterSetUpdateStart(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,err,error,*999)
!                 CALL Field_ParameterSetUpdateFinish(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, & 
!                   & FIELD_VALUES_SET_TYPE,err,error,*999)
!             CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE(PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
!               ! do nothing ???
!             CASE DEFAULT
!               localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
!                 & " is not valid for a diffusion equation type of a classical field problem class."
!             CALL FlagError(localError,err,error,*999)
!           END SELECT
!         ELSE
!           CALL FlagError("Problem is not associated.",err,error,*999)
!         ENDIF
!       ELSE
!         CALL FlagError("Solver is not associated.",err,error,*999)
!       ENDIF
!     ELSE
!       CALL FlagError("Control loop is not associated.",err,error,*999)
    !     ENDIF
    
    EXITS("Diffusion_PreSolveUpdateBoundaryConditions")
    RETURN
999 ERRORS("Diffusion_PreSolveUpdateBoundaryConditions",err,error)
    EXITS("Diffusion_PreSolveUpdateBoundaryConditions")
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveUpdateBoundaryConditions

  !   
  !================================================================================================================================
  !
  
  !>Updates the boundary conditions and source term to the required analytic values
  SUBROUTINE Diffusion_PreSolveUpdateAnalyticValues(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionType,componentIdx,derivativeIdx,dimensionIdx, &
      & dynamicVariableType,equationsSetIdx,globalDerivativeIndex,globalDofIdx,localDofIdx,nodeIdx,numberOfDimensions
    REAL(DP) :: A1,currentTime,D1,timeIncrement,normal(3),tangents(3,3),VALUE,X(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: boundaryConditions
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: boundaryConditionsVariable
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(DOMAIN_TYPE), POINTER :: domain
    TYPE(DOMAIN_NODES_TYPE), POINTER :: domainNodes
    TYPE(DOMAIN_TOPOLOGY_TYPE), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet
    TYPE(FIELD_TYPE), POINTER :: analyticField,dependentField,geometricField,materialsField,sourceField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: dynamicVariable,geometricVariable,sourceVariable
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: solverEquations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolveUpdateAnalyticValues",err,error,*999)

    A1=0.4_DP
    D1=1.0_DP

    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for a diffusion problem.",err,error,*999)
      
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !loop over all the equation sets and set the appropriate field variable type BCs and
      !the source field associated with each equation set
      DO equationsSetIdx=1,solverMapping%NUMBER_OF_EQUATIONS_SETS
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(ASSOCIATED(equationsSet%analytic)) THEN
          CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)          
          CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
          NULLIFY(geometricVariable)
          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
          NULLIFY(geometricParameters)
          CALL Field_ParameterSetDataGet(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,geometricParameters, &
            & err,error,*999)
          NULLIFY(analyticField)
          NULLIFY(analyticParameters)
          CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
          IF(ASSOCIATED(analyticField)) &
            & CALL Field_ParameterSetDataGet(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & analyticParameters,err,error,*999)
          NULLIFY(materialsField)
          NULLIFY(materialsParameters)
          CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
          IF(ASSOCIATED(materialsField)) &
            & CALL Field_ParameterSetDataGet(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999) 
          NULLIFY(equations)
          CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
          NULLIFY(vectorEquations)
          CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
          NULLIFY(vectorMapping)
          CALL EquationsVector_VectorMappingGet(vectorEquations,vectorMapping,err,error,*999)
          NULLIFY(dynamicMapping)
          CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
          NULLIFY(dynamicVariable)
          CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,dynamicVariable,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dynamicVariable,boundaryConditionsVariable,err,error,*999)
          dynamicVariableType=dynamicVariable%VARIABLE_TYPE
          analyticFunctionType=equationsSet%analytic%ANALYTIC_FUNCTION_TYPE          
          DO componentIdx=1,dynamicVariable%NUMBER_OF_COMPONENTS
            IF(dynamicVariable%components(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
              NULLIFY(domain)
              CALL FieldVariable_DomainGet(dynamicVariable,componentIdx,domain,err,error,*999)
              NULLIFY(domainTopology)
              CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
              NULLIFY(domainNodes)
              CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
              !Loop over the local nodes excluding the ghosts.
              DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
                DO dimensionIdx=1,numberOfDimensions
                  !Default to version 1 of each node derivative
                  localDofIdx=geometricVariable%components(dimensionIdx)%PARAM_TO_DOF_MAP% &
                    & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(1)%versions(1)
                  x(dimensionIdx)=geometricParameters(localDofIdx)
                ENDDO !dimensionIdx
                !Loop over the derivatives
                DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                  globalDerivativeIndex=domainNodes%nodes(nodeIdx)%derivatives(derivativeIdx)%GLOBAL_DERIVATIVE_INDEX
                  CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,currentTime, &
                    & dynamicVariableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters, &
                    & VALUE,err,error,*999)
                  !Default to version 1 of each node derivative
                  localDofIdx=dynamicVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                    & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(derivativeIdx)%VERSIONS(1)
                  globalDofIdx=dynamicVariable%DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(localDofIdx)
                  CALL Field_ParameterSetUpdateLocalDof(dependentField,dynamicVariableType, &
                    & FIELD_ANALYTIC_VALUES_SET_TYPE,localDofIdx,VALUE,err,error,*999)
                  boundaryConditionType=boundaryConditionsVariable%DOF_TYPES(globalDofIdx)
                  IF(boundaryConditionType==BOUNDARY_CONDITION_FIXED) THEN
                    CALL Field_ParameterSetUpdateLocalDof(dependentField,dynamicVariableType,FIELD_VALUES_SET_TYPE, &
                      & localDofIdx,value,err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ELSE
              CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
            ENDIF 
          ENDDO !componentIdx
          CALL Field_ParameterSetUpdateStart(dependentField,dynamicVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateStart(dependentField,dynamicVariableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,dynamicVariableType,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL Field_ParameterSetUpdateFinish(dependentField,dynamicVariableType,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(problem%specification(3)==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
            !>Set the source field to a specified analytical function
            NULLIFY(sourceField)
            CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
            NULLIFY(sourceVariable)
            CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
            DO componentIdx=1,sourceVariable%NUMBER_OF_COMPONENTS
              IF(sourceVariable%components(componentIdx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                NULLIFY(domain)
                CALL FieldVariable_DomainGet(dynamicVariable,componentIdx,domain,err,error,*999)
                NULLIFY(domainTopology)
                CALL Domain_TopologyGet(domain,domainTopology,err,error,*999)
                NULLIFY(domainNodes)
                CALL DomainTopology_NodesGet(domainTopology,domainNodes,err,error,*999)
                !Loop over the local nodes excluding the ghosts.
                DO nodeIdx=1,domainNodes%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
                  DO dimensionIdx=1,numberOfDimensions
                    !Default to version 1 of each node derivative
                    localDofIdx=geometricVariable%components(dimensionIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                      & nodes(nodeIdx)%derivatives(1)%versions(1)
                    X(dimensionIdx)=geometricParameters(localDofIdx)
                  ENDDO !dimensionIdx
                  !Loop over the derivatives
                  DO derivativeIdx=1,domainNodes%nodes(nodeIdx)%NUMBER_OF_DERIVATIVES
                    SELECT CASE(analyticFunctionType)
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
                      VALUE=-1*A1*EXP(-1*currentTime)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+6)
                    CASE DEFAULT
                      localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*", err,error))// &
                        & " is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                    !Default to version 1 of each node derivative
                    localDofIdx=sourceVariable%components(componentIdx)%PARAM_TO_DOF_MAP% &
                      & NODE_PARAM2DOF_MAP%nodes(nodeIdx)%derivatives(derivativeIdx)%versions(1)
                    CALL Field_ParameterSetUpdateLocalDof(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                      & localDofIdx,VALUE,err,error,*999)
                  ENDDO !derivativeIdx
                ENDDO !nodeIdx
              ELSE
                CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
              ENDIF
            ENDDO !componentIdx
            CALL Field_ParameterSetUpdateStart(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          ENDIF
          IF(ASSOCIATED(materialsField)) & 
            & CALL Field_ParameterSetDataRestore(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          IF(ASSOCIATED(analyticField)) &
            & CALL Field_ParameterSetDataRestore(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & analyticParameters,err,error,*999)
          CALL Field_ParameterSetDataRestore(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & geometricParameters,err,error,*999)
        ENDIF !Analytic field
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Diffusion_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORS("Diffusion_PreSolveUpdateAnalyticValues",err,error)
    EXITS("Diffusion_PreSolveUpdateAnalyticValues")
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveUpdateAnalyticValues

  !
  !================================================================================================================================
  !
  !>Update mesh position and velocity for ALE diffusion problem
  SUBROUTINE Diffusion_PreSolveALEUpdateMesh(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: GEOMETRIC_FIELD
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_ALE_DIFFUSION !<A pointer to the solvers
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,ALPHA
    REAL(DP), POINTER :: MESH_DISPLACEMENT_VALUES(:)

    INTEGER(INTG) :: dof_number,TOTAL_NUMBER_OF_DOFS,NDOFS_TO_PRINT

    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS
    INTEGER(INTG) :: INPUT_TYPE,INPUT_OPTION
    REAL(DP), POINTER :: INPUT_DATA1(:)

    ENTERS("Diffusion_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(SOLVER_ALE_DIFFUSION)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVER_MAPPING)
    NULLIFY(EQUATIONS_SET)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%LOOP_TYPE==PROBLEM_CONTROL_TIME_LOOP_TYPE) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ELSE IF(CONTROL_LOOP%CONTROL_LOOP_LEVEL>1) THEN
        CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP%PARENT_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      ENDIF
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                  & PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
                SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                  SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
                        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
                      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
                        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
                          & err,error,*999)
                      END IF
                      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
                          ! do nothing ???
                        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
                          & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
                          CALL WriteString(GENERAL_OUTPUT_TYPE,"Diffusion update mesh ... ",err,error,*999)
                          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                            !--- First, read mesh displacement values from file

                           CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                            & NUMBER_OF_DIMENSIONS,err,error,*999)

                           INPUT_TYPE=42
                           INPUT_OPTION=2
                           NULLIFY(INPUT_DATA1)
                           !CALL Field_ParameterSetDataGet(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                            !& FIELD_VALUES_SET_TYPE,INPUT_DATA1,err,error,*999)
!                            CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, &
!                             & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_TIME)

                            NULLIFY(MESH_DISPLACEMENT_VALUES)
                            CALL Field_ParameterSetDataGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                            IF(DIAGNOSTICS1) THEN
                              NDOFS_TO_PRINT = SIZE(MESH_DISPLACEMENT_VALUES,1)
                              CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,&
                                & MESH_DISPLACEMENT_VALUES,'(" MESH_DISPLACEMENT_VALUES = ",3(X,E13.6))','3(3(X,E13.6))', &
                                & err,error,*999)
                            ENDIF

!                            CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,INPUT_DATA1, &
!                             & NUMBER_OF_DIMENSIONS,INPUT_TYPE,INPUT_OPTION,CURRENT_TIME)

                            TOTAL_NUMBER_OF_DOFS = GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & TOTAL_NUMBER_OF_DOFS

                            !--- Second, update geometric field
                            DO dof_number=1,TOTAL_NUMBER_OF_DOFS
                              CALL FIELD_PARAMETER_SET_ADD_LOCAL_DOF(GEOMETRIC_FIELD, &
                                & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,dof_number, &
                                & MESH_DISPLACEMENT_VALUES(dof_number), &
                                & err,error,*999)
                            END DO
                            CALL Field_ParameterSetUpdateStart(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)
                            CALL Field_ParameterSetUpdateFinish(GEOMETRIC_FIELD, &
                              & FIELD_U_VARIABLE_TYPE, FIELD_VALUES_SET_TYPE,err,error,*999)

                            !--- Third, use displacement values to calculate velocity values
                            ALPHA=1.0_DP/TIME_INCREMENT
                            CALL FIELD_PARAMETER_SETS_COPY(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE,ALPHA,err,error,*999)
                            CALL Field_ParameterSetDataRestore(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & FIELD_MESH_DISPLACEMENT_SET_TYPE,MESH_DISPLACEMENT_VALUES,err,error,*999)
                          ELSE
                            CALL FlagError("Geometric field is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="Equations set subtype " &
                            & //TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is not valid for a diffusion equation type of a classical field problem class."
                          CALL FlagError(localError,err,error,*999)
                      END SELECT
                    ELSE
                      CALL FlagError("Equations set is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Solver mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Solver equations is not associated.",err,error,*999)
                ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
              CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Diffusion_PreSolveALEUpdateMesh")
    RETURN
999 ERRORSEXITS("Diffusion_PreSolveALEUpdateMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveALEUpdateMesh
  !   
  !================================================================================================================================
  !
  SUBROUTINE Diffusion_PreSolveStoreCurrentSolution(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: I

    ENTERS("Diffusion_PreSolveStoreCurrentSolution",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion-one equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value diffusion-one dependent field at time, t ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      DEPENDENT_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_ONE)) THEN
                        CALL Field_NumberOfComponentsGet(DEPENDENT_FIELD_DIFFUSION_ONE, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_DIFFUSION_ONE is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion-one equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion-one solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion-one solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,err,error,*999)
                  END DO


!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL Field_ParameterSetDataGet(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL Field_ParameterSetDataRestore(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==2) THEN
                !--- Get the dependent field of the diffusion equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value of diffusion solution &
                   & (dependent field - V variable_type) at time, t ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION_ONE,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      DEPENDENT_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_ONE)) THEN
                        CALL Field_NumberOfComponentsGet(DEPENDENT_FIELD_DIFFUSION_ONE, &
                          & FIELD_V_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_DIFFUSION_ONE is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the current time value parameters set from diffusion-one's dependent field 
                  DO I=1,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_ONE
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,DEPENDENT_FIELD_DIFFUSION_ONE, & 
                      & FIELD_V_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,I,err,error,*999)
                  END DO


!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL Field_ParameterSetDataGet(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL Field_ParameterSetDataRestore(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Diffusion_PreSolveStoreCurrentSolution")
    RETURN
999 ERRORS("Diffusion_PreSolveStoreCurrentSolution",err,error)
    EXITS("Diffusion_PreSolveStoreCurrentSolution")
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveStoreCurrentSolution
  
  !   
  !================================================================================================================================
  !
  SUBROUTINE Diffusion_PreSolveGetSourceValue(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER_DIFFUSION_ONE, SOLVER_DIFFUSION_TWO  !<A pointer to the solvers
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD_DIFFUSION_TWO, SOURCE_FIELD_DIFFUSION_ONE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS_DIFFUSION_ONE, SOLVER_EQUATIONS_DIFFUSION_TWO  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING_DIFFUSION_ONE, SOLVER_MAPPING_DIFFUSION_TWO !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET_DIFFUSION_ONE, EQUATIONS_SET_DIFFUSION_TWO !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    INTEGER(INTG) :: NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
    INTEGER(INTG) :: I


    ENTERS("Diffusion_PreSolveGetSourceValue",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      NULLIFY(SOLVER_DIFFUSION_ONE)
      NULLIFY(SOLVER_DIFFUSION_TWO)

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_COUPLED_SOURCE_DIFFUSION_DIFFUSION_SUBTYPE)
              IF(SOLVER%GLOBAL_NUMBER==1) THEN
                !--- Get the dependent field of the diffusion_two equations
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Update diffusion-one source field ... ",err,error,*999)
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,2,SOLVER_DIFFUSION_TWO,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION_TWO=>SOLVER_DIFFUSION_TWO%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_TWO)) THEN
                  SOLVER_MAPPING_DIFFUSION_TWO=>SOLVER_EQUATIONS_DIFFUSION_TWO%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_TWO)) THEN
                    EQUATIONS_SET_DIFFUSION_TWO=>SOLVER_MAPPING_DIFFUSION_TWO%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_TWO)) THEN
                      DEPENDENT_FIELD_DIFFUSION_TWO=>EQUATIONS_SET_DIFFUSION_TWO%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(DEPENDENT_FIELD_DIFFUSION_TWO)) THEN
                        CALL Field_NumberOfComponentsGet(DEPENDENT_FIELD_DIFFUSION_TWO, &
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO,err,error,*999)
                      ELSE
                        CALL FlagError("DEPENDENT_FIELD_DIFFUSION_TWO is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion-two equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion-two solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion-two solver equations are not associated.",err,error,*999)
                END IF


                !--- Get the source field for the diffusion_one equations
                CALL SOLVERS_SOLVER_GET(SOLVER%SOLVERS,1,SOLVER_DIFFUSION_ONE,err,error,*999)
                SOLVER_EQUATIONS_DIFFUSION_ONE=>SOLVER_DIFFUSION_ONE%SOLVER_EQUATIONS
                IF(ASSOCIATED(SOLVER_EQUATIONS_DIFFUSION_ONE)) THEN
                  SOLVER_MAPPING_DIFFUSION_ONE=>SOLVER_EQUATIONS_DIFFUSION_ONE%SOLVER_MAPPING
                  IF(ASSOCIATED(SOLVER_MAPPING_DIFFUSION_ONE)) THEN
                    EQUATIONS_SET_DIFFUSION_ONE=>SOLVER_MAPPING_DIFFUSION_ONE%EQUATIONS_SETS(1)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET_DIFFUSION_ONE)) THEN
                      SOURCE_FIELD_DIFFUSION_ONE=>EQUATIONS_SET_DIFFUSION_ONE%SOURCE%SOURCE_FIELD
                      IF(ASSOCIATED(SOURCE_FIELD_DIFFUSION_ONE)) THEN
                        CALL Field_NumberOfComponentsGet(SOURCE_FIELD_DIFFUSION_ONE, & 
                          & FIELD_U_VARIABLE_TYPE,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE,err,error,*999)
                      ELSE
                        CALL FlagError("SOURCE_FIELD_DIFFUSION_ONE is not associated.",err,error,*999)
                      END IF
                    ELSE
                      CALL FlagError("Diffusion-one equations set is not associated.",err,error,*999)
                    END IF
                  ELSE
                    CALL FlagError("Diffusion-one solver mapping is not associated.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("Diffusion-one solver equations are not associated.",err,error,*999)
                END IF

                !--- Copy the result from diffusion-two's dependent field to diffusion-one's source field
                IF(NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE==NUMBER_OF_COMPONENTS_DEPENDENT_FIELD_DIFFUSION_TWO) THEN
                  DO I=1,NUMBER_OF_COMPONENTS_SOURCE_FIELD_DIFFUSION_ONE
                    CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD_DIFFUSION_TWO, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,SOURCE_FIELD_DIFFUSION_ONE, & 
                      & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,I,err,error,*999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIELDMESHDISPLACEMENTTYPE needs to be changed to appropriate type for this problem
                  END DO
                ELSE
                  localError="Number of components of diffusion-two dependent field "// &
                    & "is not consistent with diffusion-one-equation source field."
                  CALL FlagError(localError,err,error,*999)
                END IF

!                 IF(DIAGNOSTICS3) THEN
!                   NULLIFY( DUMMY_VALUES2 )
!                   CALL Field_ParameterSetDataGet(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                   NDOFS_TO_PRINT = SIZE(DUMMY_VALUES2,1)
!                   CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,NDOFS_TO_PRINT,NDOFS_TO_PRINT,NDOFS_TO_PRINT,DUMMY_VALUES2, &
!                     & '(" DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))',&
!                     & '4(4(X,E13.6))',err,error,*999)
!                   CALL Field_ParameterSetDataRestore(DEPENDENT_FIELD_FINITE_ELASTICITY,FIELD_U_VARIABLE_TYPE, &
!                     & FIELD_VALUES_SET_TYPE,DUMMY_VALUES2,err,error,*999)
!                 ENDIF

              END IF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("Diffusion_PreSolveGetSourceValue")
    RETURN
999 ERRORSEXITS("Diffusion_PreSolveGetSourceValue",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveGetSourceValue
  
  !   
  !================================================================================================================================
  !
  
  !>Performs post-solve operations for a diffusion problem.
  SUBROUTINE Diffusion_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: controlLoop
    TYPE(PROBLEM_TYPE), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PostSolve",err,error,*999)
    
    IF(.NOT.ASSOCIATED(solver)) CALL FlagError("Solver is not associated.",err,error,*999)
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    IF(.NOT.ALLOCATED(problem%specification)) CALL FlagError("Problem specification is not allocated.",err,error,*999)
    IF(SIZE(problem%specification,1)<3) &
      & CALL FlagError("Problem specification must have at least three entries for a diffusion problem.",err,error,*999)
    
    SELECT CASE(problem%specification(3))
    CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
      !CALL Diffusion_PostSolveOuputData(controlLoop,solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(problem%specification(3),"*",err,error))// &
        & " is not valid for a diffusion type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Diffusion_PostSolve")
    RETURN
999 ERRORSEXITS("Diffusion_PostSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PostSolve
  
  !   
  !================================================================================================================================
  !
  
  !>Output data post solve
  SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(VARYING_STRING) :: localError

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(14) :: FILE
    CHARACTER(14) :: OUTPUT_FILE

    ENTERS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
!       write(*,*)'CURRENT_TIME = ',CURRENT_TIME
!       write(*,*)'TIME_INCREMENT = ',TIME_INCREMENT
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE,PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr

                    CURRENT_LOOP_ITERATION=CONTROL_LOOP%TIME_LOOP%ITERATION_NUMBER
                    OUTPUT_ITERATION_NUMBER=CONTROL_LOOP%TIME_LOOP%OUTPUT_NUMBER

                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CONTROL_LOOP%TIME_LOOP%CURRENT_TIME<=CONTROL_LOOP%TIME_LOOP%STOP_TIME) THEN
                        IF(CURRENT_LOOP_ITERATION<10) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_000",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<100) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_00",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<1000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_0",I0)') CURRENT_LOOP_ITERATION
                        ELSE IF(CURRENT_LOOP_ITERATION<10000) THEN
                          WRITE(OUTPUT_FILE,'("TIME_STEP_",I0)') CURRENT_LOOP_ITERATION
                        END IF
                        FILE=OUTPUT_FILE
  !          FILE="TRANSIENT_OUTPUT"
!!!!!!!!ADAPT THIS TO WORK WITH DIFFUSION AND NOT JUST FLUID MECHANICS
!                         METHOD="FORTRAN"
!                         EXPORT_FIELD=.TRUE.
!                         IF(EXPORT_FIELD) THEN          
!                           IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN   
!                             CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
!                             CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
!                             CALL FLUID_MECHANICS_IO_WRITE_CMGUI(EQUATIONS_SET%REGION,EQUATIONS_SET%GLOBAL_NUMBER,FILE, &
!                               & err,error,*999)
!                             CALL WriteString(GENERAL_OUTPUT_TYPE,OUTPUT_FILE,err,error,*999)
!                             CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
!                           ENDIF
!                         ENDIF 

                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE==EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1 .OR. &
                            & EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE== &
                            & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1) THEN
                            CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,err,error,*999)
                          ENDIF
                        ENDIF
                      ENDIF 
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE)
              ! do nothing ???
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a diffusion equation type of a classical field problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF
    EXITS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !
  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,mh,mhs,ms,ng,nh,nhs,ni,nj,ns,my_compartment,Ncompartments,imatrix,num_var_count
    INTEGER(INTG) :: MESH_COMPONENT_1, MESH_COMPONENT_2
    REAL(DP) :: C_PARAM,K_PARAM,RWG,SUM,PGMJ(3),PGNJ(3),A_PARAM,COUPLING_PARAM,PGM,PGN
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS_1, DEPENDENT_BASIS_2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD,SOURCE_FIELD,EQUATIONS_SET_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(FIELD_VARIABLE_PTR_TYPE) :: FIELD_VARIABLES(99)
    TYPE(EquationsMatrixPtrType) :: COUPLING_MATRICES(99) 
    INTEGER(INTG) :: FIELD_VAR_TYPES(99)
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME_1, QUADRATURE_SCHEME_2
    TYPE(VARYING_STRING) :: localError
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS, &
      & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT, &
      & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT
    INTEGER(INTG), POINTER :: EQUATIONS_SET_FIELD_DATA(:)
     
    ENTERS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE, EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            SOURCE_FIELD=>equations%interpolation%sourceField
          ENDIF
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          rhsVector=>vectorMatrices%rhsVector
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            sourceVector=>vectorMatrices%sourceVector
          ENDIF
          vectorMapping=>vectorEquations%vectorMapping
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
            & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN  
            ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS=> &
              & equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATION_PARAMETERS,err,error,*999)
            ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT=> &
              & equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
              & equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,err,error,*999)
            DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
              & equations%interpolation%dependentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr
          ENDIF
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            ENDIF
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                        ELSEIF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                          A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS,NO_PART_DERIV)
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG- &
                            & A_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
              & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
              IF(sourceVector%updateVector) THEN
                C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ELSEIF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
              IF(sourceVector%updateVector) THEN
                !The value of the source term is +0.5*(C_1^{t}+C_1_{t+1}-C_2^{t}) 
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT,err,error,*999)
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,err,error,*999)
                write(*,*) ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                write(*,*) DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                C_PARAM=0.5_DP*ADVEC_DIFF_DEPENDENT_CURRENT_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)- &
                  & DIFFUSION_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                !                     C_PARAM_1_T= equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)!<This is the value of the solution from the advection-diffusion equation at time T
                !                     C_PARAM_1_TPLUSONE= !<This is the value of the solution from the advection-diffusion equation at time T+deltaT
                !                     C_PARAM_2_T= equations%interpolation%dependentInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)!<This is the value of the solution from the diffusion equation at time T
                !                     C_PARAM=C_PARAM_1_T+C_PARAM_1_TPLUSONE+C_PARAM_2_T
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP 
          ENDDO !ng

          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
            ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS=> &
              & equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_PREVIOUS_VALUES_SET_TYPE,ELEMENT_NUMBER, &
              & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATION_PARAMETERS,err,error,*999)
            ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT=> &
              & equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              IF(sourceVector%updateVector) THEN
                !The value of the source term is +0.5*(C_1^{t}+C_1_{t+1}-C_2^{t}) 
                CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
                  & ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT,err,error,*999)
                write(*,*) ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                C_PARAM=0.5_DP*ADVEC_DIFF_DEPENDENT_PREVIOUS_INTERPOLATED_POINT%VALUES(1,NO_PART_DERIV)
                mhs=0
                DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                      & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh   
              ENDIF
              IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP 
            ENDDO !ng
          ENDIF

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
                  & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE) THEN
                  IF(sourceVector%updateVector) sourceVector%elementVector%vector(mhs)= & 
                    & sourceVector%elementVector%vector(mhs)* &
                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
!!!!!!!!!!!!!!!MULTI-COMPARTMENT DIFFUSION - PROTOTYPE FOR OTHER MULTI-COMPARTMENT MODELS IN FUTURE
!!!!!!!!!!!!!!!HAS BEEN SEPARATED HERE FOR EASE OF DEVELOPMENT & READABILITY OF THIS NEW FEATURE
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)

          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          SOURCE_FIELD=>equations%interpolation%sourceField
          EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD

          vectorEquations=>equations%vectorEquations
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          rhsVector=>vectorMatrices%rhsVector
          sourceVector=>vectorMatrices%sourceVector
          stiffnessMatrix%elementMatrix%matrix=0.0_DP
          dampingMatrix%elementMatrix%matrix=0.0_DP
          vectorMapping=>vectorEquations%vectorMapping


          CALL Field_ParameterSetDataGet(EQUATIONS_SET_FIELD,FIELD_U_VARIABLE_TYPE, &
            & FIELD_VALUES_SET_TYPE,EQUATIONS_SET_FIELD_DATA,err,error,*999)

          my_compartment = EQUATIONS_SET_FIELD_DATA(1)
          Ncompartments  = EQUATIONS_SET_FIELD_DATA(2)

          linearMatrices=>vectorMatrices%linearMatrices
          linearMapping=>vectorMapping%linearMapping


          num_var_count=0
          DO imatrix = 1,Ncompartments
            IF(imatrix/=my_compartment)THEN
              num_var_count=num_var_count+1
              COUPLING_MATRICES(num_var_count)%ptr=>linearMatrices%matrices(num_var_count)%ptr
              FIELD_VARIABLES(num_var_count)%ptr=>linearMapping%equationsMatrixToVarMaps(num_var_count)%VARIABLE
              FIELD_VAR_TYPES(num_var_count)=FIELD_VARIABLES(num_var_count)%ptr%VARIABLE_TYPE
              COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix=0.0_DP
            ENDIF
          END DO


          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & sourceInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        SUM=0.0_DP
                        DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                          PGMJ(nj)=0.0_DP
                          PGNJ(nj)=0.0_DP
                          DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI                          
                            PGMJ(nj)=PGMJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            PGNJ(nj)=PGNJ(nj)+ &
                              & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                              & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                          ENDDO !ni
                          K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                            & VALUES(nj,NO_PART_DERIV)
                          SUM=SUM+K_PARAM*PGMJ(nj)*PGNJ(nj)
                        ENDDO !nj
                        COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                          & VALUES(my_compartment,NO_PART_DERIV)
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+ & 
                          & SUM*RWG + QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG*COUPLING_PARAM
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
            IF(sourceVector%updateVector) THEN
              C_PARAM=equations%interpolation%sourceInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1, NO_PART_DERIV)
              mhs=0
              DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  sourceVector%elementVector%vector(mhs)=sourceVector%elementVector%vector(mhs)+ &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
            IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP 
            !Calculate the coupling matrices

            !Loop over element rows
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !field_variable is the variable associated with the equations set under consideration

              MESH_COMPONENT_1 = FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS_1 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_1)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME_1 => DEPENDENT_BASIS_1%QUADRATURE% &
                & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              RWG = equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN * &
                & QUADRATURE_SCHEME_1%GAUSS_WEIGHTS(ng)

              DO ms=1,DEPENDENT_BASIS_1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1

                num_var_count=0
                DO imatrix = 1,Ncompartments
                  IF(imatrix/=my_compartment)THEN
                    num_var_count=num_var_count+1

                    !need to test for the case where imatrix==mycompartment
                    !the coupling terms then needs to be added into the stiffness matrix
                    IF(COUPLING_MATRICES(num_var_count)%ptr%updateMatrix) THEN

                      !                       !Loop over element columns
                      nhs=0
                      ! !                       DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      DO nh=1,FIELD_VARIABLES(num_var_count)%ptr%NUMBER_OF_COMPONENTS

                        MESH_COMPONENT_2 = FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                        DEPENDENT_BASIS_2 => DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(MESH_COMPONENT_2)%ptr% &
                          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                        !--- We cannot use two different quadrature schemes here !!!
                        QUADRATURE_SCHEME_2 => DEPENDENT_BASIS_2%QUADRATURE% &
                          & QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                        !RWG = equations%interpolation%geometricInterpPointMetrics%JACOBIAN * &
                        !  & QUADRATURE_SCHEME_2%GAUSS_WEIGHTS(ng)

                        DO ns=1,DEPENDENT_BASIS_2%NUMBER_OF_ELEMENT_PARAMETERS
                          nhs=nhs+1

                          !                           !-------------------------------------------------------------------------------------------------------------
                          !                           !concentration test function, concentration trial function
                          !                           !For now, this is only a dummy implementation - this still has to be properly set up.
                          !                           IF(mh==nh.AND.nh<NUMBER_OF_VEL_PRESS_COMPONENTS) THEN ! don't need this for diffusion equation

                          !                             SUM = 0.0_DP

                          PGM=QUADRATURE_SCHEME_1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PGN=QUADRATURE_SCHEME_2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)

                          !Get the coupling coefficients 
                          COUPLING_PARAM=equations%interpolation%materialsInterpPoint(FIELD_V_VARIABLE_TYPE)%ptr% &
                            & VALUES(imatrix,NO_PART_DERIV)

                          !                              SUM = SUM + COUPLING_PARAM * PGM * PGN

                          COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) = &
                            & COUPLING_MATRICES(num_var_count)%ptr%elementMatrix%matrix(mhs,nhs) + & 
                            & COUPLING_PARAM * PGM * PGN * RWG
                          !                           ENDIF

                        ENDDO !ns
                      ENDDO !nh
                    ENDIF
                  ENDIF
                ENDDO !imatrix
              ENDDO !ms
            ENDDO !mh


          ENDDO !ng

          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                IF(sourceVector%updateVector) sourceVector%elementVector%vector(mhs)= & 
                  & sourceVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF

        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE,EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE,EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a diffusion equation problem.
  SUBROUTINE Diffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Diffusion_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE, &
          & PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
          & PROBLEM_NONLINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          !All ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a diffusion type of a classical field problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_DIFFUSION_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Diffusion equation problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("Diffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("Diffusion_ProblemSpecificationSet",err,error)
    EXITS("Diffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Diffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: PROBLEM_SUBTYPE
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",err,error,*999)


    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
      END IF
      PROBLEM_SUBTYPE=PROBLEM%SPECIFICATION(3)
      IF(PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_NO_SOURCE_ALE_DIFFUSION_SUBTYPE .OR. &
         & PROBLEM_SUBTYPE==PROBLEM_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Dynamic solver",err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a linear diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a linear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM_SUBTYPE,"*",err,error))// &
          & " does not equal a linear diffusion equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_LINEAR_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the nonlinear diffusion problem
  SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Diffusion problem.",err,error,*999)
      END IF
      IF( PROBLEM%SPECIFICATION(3)==PROBLEM_NONLINEAR_SOURCE_DIFFUSION_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a time control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,PROBLEM_CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,err,error,*999)
            !Set the solver to be a first order dynamic solver 
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Dynamic solver",err,error,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear diffusion problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a nonlinear diffusion problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a nonlinear diffusion problem subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP",err,error)
    RETURN 1
  END SUBROUTINE DIFFUSION_EQUATION_PROBLEM_NONLINEAR_SETUP
  
  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a diffusion equation finite element equations set.
  SUBROUTINE Diffusion_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)
    
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element Jacobian evaluation on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nh,nhs,ns
    REAL(DP) :: B_PARAM,C_PARAM,RWG,U_VALUE,VALUE
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Diffusion_FiniteElementJacobianEvaluate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      NULLIFY(vectorEquations)
      CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
      IF(ASSOCIATED(EQUATIONS)) THEN
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with no source.",err,error,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with a constant source.",err,error,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with a linear source.",err,error,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
          vectorMatrices=>vectorEquations%vectorMatrices
          nonlinearMatrices=>vectorMatrices%nonlinearMatrices
          jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
          IF(jacobianMatrix%updateJacobian) THEN
            !Store all these in equations matrices/somewhere else?????
            DEPENDENT_FIELD=>equations%interpolation%dependentField
            GEOMETRIC_FIELD=>equations%interpolation%geometricField
            MATERIALS_FIELD=>equations%interpolation%materialsField
            vectorMapping=>vectorEquations%vectorMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            DEPENDENT_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
            GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              !Calculate RWG.
              RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              !Find material parameters and u value at this Gauss point
              B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              !Loop over field components
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      VALUE=-2.0_DP*B_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*U_VALUE
                      jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+VALUE*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh              
            ENDDO !ng
          ENDIF
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)                 
          vectorMatrices=>vectorEquations%vectorMatrices
          nonlinearMatrices=>vectorMatrices%nonlinearMatrices
          jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
          IF(jacobianMatrix%updateJacobian) THEN
            !Store all these in equations matrices/somewhere else?????
            DEPENDENT_FIELD=>equations%interpolation%dependentField
            GEOMETRIC_FIELD=>equations%interpolation%geometricField
            MATERIALS_FIELD=>equations%interpolation%materialsField
            vectorMapping=>vectorEquations%vectorMapping
            nonlinearMapping=>vectorMapping%nonlinearMapping
            DEPENDENT_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
            FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
            GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
            DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
              & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
            QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
              & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              !Calculate RWG.
              RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
              !Find material parameter and u value at this Gauss point
              B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              C_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              !Loop over field components
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      VALUE=-B_PARAM*C_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                        & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*EXP(C_PARAM*U_VALUE)
                      jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+VALUE*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDDO !ms
              ENDDO !mh              
            ENDDO !ng
          ENDIF
        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with no source.",err,error,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with a constant source.",err,error,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with a linear source.",err,error,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a multi component transport diffusion equation.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Diffusion_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("Diffusion_FiniteElementJacobianEvaluate",err,error)
    EXITS("Diffusion_FiniteElementJacobianEvaluate")
    RETURN 1
    
  END SUBROUTINE Diffusion_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Diffusion equation finite element equations set.
  SUBROUTINE Diffusion_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nj,nh,nhs,ni,ns
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,K_PARAM,RWG,SUM1,SUM2,PGMJ(3),PGNJ(3),U_VALUE
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Diffusion_FiniteElementResidualEvaluate",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a diffusion type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_NO_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with no source.",err,error,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with a constant source.",err,error,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a diffusion equation with a linear source.",err,error,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          nonlinearMatrices=>vectorMatrices%nonlinearMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(1)%ptr
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping
          dynamicMapping=>vectorMapping%dynamicMapping
          nonlinearMapping=>vectorMapping%nonlinearMapping
          DEPENDENT_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            IF(stiffnessMatrix%firstAssembly.OR.dampingMatrix%firstAssembly.OR.rhsVector%firstAssembly) THEN
              B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(stiffnessMatrix%updateMatrix) THEN
                          SUM1=0.0_DP
                          DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                            PGMJ(nj)=0.0_DP
                            PGNJ(nj)=0.0_DP
                            DO ni=1,GEOMETRIC_BASIS%NUMBER_OF_XI
                              PGMJ(nj)=PGMJ(nj)+ &
                                & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                              PGNJ(nj)=PGNJ(nj)+ &
                                & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            ENDDO !ni
                            K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & VALUES(nj,NO_PART_DERIV)
                            SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          ENDDO !nj
                          SUM2=B_PARAM*QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)                        
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+ &
                            & (SUM1+SUM2)*RWG
                        ENDIF
                        IF(dampingMatrix%updateMatrix) THEN
                          dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                ENDDO !ms
              ENDDO !mh
            ENDIF
            IF(rhsVector%firstAssembly) THEN
              IF(rhsVector%updateVector) THEN
                A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                  & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1,NO_PART_DERIV)
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*A_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(nonlinearMatrices%updateResidual) THEN
              C_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(1,NO_PART_DERIV)
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)- &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*C_PARAM*U_VALUE**2*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ENDDO !ng
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)                 
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          nonlinearMatrices=>vectorMatrices%nonlinearMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(1)%ptr
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping
          dynamicMapping=>vectorMapping%dynamicMapping
          nonlinearMapping=>vectorMapping%nonlinearMapping
          DEPENDENT_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            IF(stiffnessMatrix%firstAssembly.OR.dampingMatrix%firstAssembly.OR.rhsVector%firstAssembly) THEN
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        IF(stiffnessMatrix%updateMatrix) THEN
                          SUM1=0.0_DP
                          DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                            PGMJ(nj)=0.0_DP
                            PGNJ(nj)=0.0_DP
                            DO ni=1,GEOMETRIC_BASIS%NUMBER_OF_XI
                              PGMJ(nj)=PGMJ(nj)+ &
                                & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                              PGNJ(nj)=PGNJ(nj)+ &
                                & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                                & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                            ENDDO !ni
                            K_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                              & VALUES(nj,NO_PART_DERIV)
                            SUM1=SUM1+K_PARAM*PGMJ(nj)*PGNJ(nj)
                          ENDDO !nj
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM1*RWG
                        ENDIF
                        IF(dampingMatrix%updateMatrix) THEN
                          dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)*RWG
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF
                  IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
                ENDDO !ms
              ENDDO !mh
            ENDIF
            IF(rhsVector%firstAssembly) THEN
              IF(rhsVector%updateVector) THEN
                A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                  & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+1,NO_PART_DERIV)
                mhs=0
                DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                  !Loop over element rows
                  DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                    mhs=mhs+1
                    rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)+ &
                       & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*A_PARAM*RWG
                  ENDDO !ms
                ENDDO !mh
              ENDIF
            ENDIF
            IF(nonlinearMatrices%updateResidual) THEN
              B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+2,NO_PART_DERIV)
              C_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS+3,NO_PART_DERIV)
              U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr% &
                & VALUES(1,NO_PART_DERIV)
!!TODO: Handle floating point exceptions better
              IF((C_PARAM*U_VALUE)>20000.0_DP) THEN
                localError="The value of "//TRIM(NumberToVString(C_PARAM*U_VALUE,"*",err,error))// &
                  & " is out of range for an exponential function."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              mhs=0
              DO mh=1,DEPENDENT_VARIABLE%NUMBER_OF_COMPONENTS
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                  mhs=mhs+1
                  nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)- &
                    & QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)*B_PARAM*EXP(C_PARAM*U_VALUE)*RWG
                ENDDO !ms
              ENDDO !mh
            ENDIF
          ENDDO !ng
        CASE(EQUATIONS_SET_NO_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with no source.",err,error,*999)
        CASE(EQUATIONS_SET_CONSTANT_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with a constant source.",err,error,*999)
        CASE(EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for an ALE diffusion equation with a linear source.",err,error,*999)
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          CALL FlagError("Can not evaluate a residual for a multi component transport diffusion equation.",err,error,*999)
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a diffusion equation type of a classical field equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Diffusion_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Diffusion_FiniteElementResidualEvaluate",err,error)
    EXITS("Diffusion_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Diffusion_FiniteElementResidualEvaluate
 
  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(CONTROL_LOOP_TIME_TYPE), POINTER :: TIME_LOOP,TIME_LOOP_PARENT
    TYPE(CONTROL_LOOP_TYPE), POINTER :: PARENT_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM
    TYPE(REGION_TYPE), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,localError,METHOD
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER,CURRENT_LOOP_ITERATION

    ENTERS("DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP",err,error,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        SELECT CASE(CONTROL_LOOP%LOOP_TYPE)
        CASE(PROBLEM_CONTROL_SIMPLE_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_FIXED_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_TIME_LOOP_TYPE)
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%TIME_LOOP
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              !Get the solver. 
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)            
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              !Loop over the equations sets associated with the solver
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      NULLIFY(DEPENDENT_REGION)
                      CALL FIELD_REGION_GET(DEPENDENT_FIELD,DEPENDENT_REGION,err,error,*999)
                      NULLIFY(PARENT_LOOP)
                      PARENT_LOOP=>CONTROL_LOOP%PARENT_LOOP
                      IF(ASSOCIATED(PARENT_LOOP)) THEN
                        !add the iteration number of the parent loop to the filename
                        NULLIFY(TIME_LOOP_PARENT)
                        TIME_LOOP_PARENT=>PARENT_LOOP%TIME_LOOP
                        IF(ASSOCIATED(TIME_LOOP_PARENT)) THEN
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP_PARENT%OUTPUT_NUMBER
                          CURRENT_LOOP_ITERATION=TIME_LOOP_PARENT%GLOBAL_ITERATION_NUMBER
                          FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%USER_NUMBER,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP_PARENT%GLOBAL_ITERATION_NUMBER,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP%ITERATION_NUMBER,"*",err,error))
                        ELSE
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP%OUTPUT_NUMBER
                          CURRENT_LOOP_ITERATION=TIME_LOOP%GLOBAL_ITERATION_NUMBER
                          FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%USER_NUMBER,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",err,error))
                        ENDIF
                      ELSE
                        OUTPUT_ITERATION_NUMBER=TIME_LOOP%OUTPUT_NUMBER
                        CURRENT_LOOP_ITERATION=TIME_LOOP%GLOBAL_ITERATION_NUMBER
                        FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%USER_NUMBER,"*",err,error))// &
                          & "_"//TRIM(NumberToVString(TIME_LOOP%GLOBAL_ITERATION_NUMBER,"*",err,error))
                      ENDIF
                      METHOD="FORTRAN"
                      IF(OUTPUT_ITERATION_NUMBER/=0.AND.MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0) THEN
                      !IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0) THEN
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,err,error,*999)
                      ENDIF
                    ELSE
                      localError="Equations set is not associated for equations set index "// &
                        & TRIM(NumberToVString(equations_set_idx,"*",err,error))// &
                        & " in the solver mapping."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FlagError("Solver equations solver mapping is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver solver equations are not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Control loop problem is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Time loop is not associated.",err,error,*999)
          ENDIF
        CASE(PROBLEM_CONTROL_WHILE_LOOP_TYPE)
          !do nothing
        CASE(PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE DEFAULT
          localError="The control loop type of "//TRIM(NumberToVString(CONTROL_LOOP%LOOP_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP",err,error)
    RETURN 1
    
  END SUBROUTINE DIFFUSION_EQUATION_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

END MODULE DIFFUSION_EQUATION_ROUTINES 
