!> \file
!> \author 
!> \brief This module handles all Burgers equation routines.
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
!> Contributor(s): David Ladd, Chris Bradley
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

!>This module handles all Burgers equation routines.
MODULE BurgersEquationsRoutines

  USE AnalyticAnalysisRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetAccessRoutines
  USE FieldRoutines
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
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

  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC Burgers_AnalyticFunctionsEvaluate

  PUBLIC Burgers_BoundaryConditionsAnalyticCalculate

  PUBLIC Burgers_EquationsSetSetup

  PUBLIC Burgers_EquationsSetSolutionMethodSet

  PUBLIC Burgers_EquationsSetSpecificationSet

  PUBLIC Burgers_FiniteElementJacobianEvaluate

  PUBLIC Burgers_FiniteElementResidualEvaluate

  PUBLIC Burgers_ProblemSpecificationSet

  PUBLIC Burgers_ProblemSetup

  PUBLIC Burgers_PreSolve,Burgers_PostSolve

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a one-dimensional dynamic solution to the burgers equation
  SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: errOR !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,derivativeIdx,dimensionIdx,localDOFIdx,nodeIdx,numberOfDimensions,variableIdx,variableType, &
      & versionIdx
    REAL(DP) :: dependentValue,X(3),initialValue
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,geometricVariable,materialsVariable
    INTEGER(INTG) :: globalDerivativeIndex,analyticFunctionType
    !THESE ARE TEMPORARY VARIABLES - they need to be replace by constant field values and the current simulation time
    REAL(DP) :: time,normal(3),tangents(3,3)
    !CURRENT_TIME = 1.2_DP

    ENTERS("Burgers_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(materialsField)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    NULLIFY(analyticField)
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
    NULLIFY(geometricParameters)
    CALL FieldVarible_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(analyticVariable)
    NULLIFY(analyticParameters)
    IF(ASSOCIATED(analyticField)) THEN
      CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    ENDIF
    NULLIFY(materialsVariable)
    NULLIFY(materialsParameters)
    IF(ASSOCIATED(materialsField)) THEN
      CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
      CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
    ENDIF
    CALL EquationsSet_AnalyticTimeGet(equationsSet,time,err,error,*999)
    CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
    DO variableIdx=1,numberOfVariables
      NULLIFY(dependentVariable)
      CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
      CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
      DO componentIdx=1,numberOfComponents
        CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
          & err,error,*999)
        NULLIFY(domain)
        CALL FieldVariable_ComponentDomainGet(dependentVariable,componentIdx,domain,err,error,*999)
        NULLIFY(domainTopology)
        CALL Domain_DomainTopologyGet(domain,domainTopology,err,error,*999)
        NULLIFY(domainNodes)
        CALL DomainTopology_DomainNodesGet(domainTopology,domainNodes,err,error,*999)
        !Loop over the local nodes excluding the ghosts.
        CALL DomainNode_NumberOfNodesGet(domainNodes,numberOfNodes,err,error,*999)
        DO nodeIdx=1,numberOfNodes
!!TODO \todo We should interpolate the geometric field here and the node position.
          DO dimensionIdx=1,numberOfDimensions
            !Default to version 1 of each node derivative
            CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
            x(dimensionIdx)=geometricParameters(localDOFIdx)
          ENDDO !dimensionIdx
          CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
          !Loop over the derivatives
          CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
          DO derivativeIdx=1,numberOfNodeDerivatives
            CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
            CALL Burgers_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,0.0_DP,variableType, &
              & globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,initialValue,err,error,*999)
            CALL Burgers_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
              & globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,dependentValue,err,error,*999)
            CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
            DO versionIdx=1,numberOfVersions
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,versionIdx,derivativeIdx,nodeIdx,dimensionIdx,localDOFIdx, &
                & err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                & dependentValue,err,error,*999)
              IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                IF(boundaryNode) THEN
                  !If we are a boundary node then set the analytic value on the boundary
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentField,variableType,localDOFIdx, &
                    & BOUNDARY_CONDITION_FIXED,dependentValue,err,error,*999)
                ELSE
                  !Set the initial condition.
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                    & initialValue,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !versionIdx
          ENDDO !derivativeIdx
        ENDDO !nodeIdx
      ENDDO !componentIdx
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
      CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
    ENDDO !variableIdx
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)

    EXITS("Burgers_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Burgers_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1

  END SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate


  !
  !================================================================================================================================
  !
  !>Evaluate the analytic solutions for a Burgers equation
  SUBROUTINE Burgers_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
    & globalDerivative,componentNumber,analyticParameters,materialsParameters,dependentValue,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to evaluate
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: x(:) !<x(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivative !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: dependentValue !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: equationsSubtype,esSpecification(3)
    REAL(DP) :: aParam,bParam,cParam,dParam,eParam,x0Param
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_AnalyticFunctionsEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    equationsSubtype=esSpecification(3)
    SELECT CASE(equationsSubtype)
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
        !For del[u]/del[t] + u.(del[u]/del[x]) = nu.(del^2[u]/del[x]^2)
        !u(x,t)=1-tanh(x-x0-t)/(2.nu))   with BCs,
        !u(0,t) = 2, u_{n} = 2.u_{n-1} - u_{n-2}
        !see http://www.cfd-online.com/Wiki/Burgers_equation
        !OpenCMISS has del[u]/del[t] + K.(del^2[u]/del[x]^2) + u.(del[u]/del[x]) = 0,
        !u(x,t)= 1 - tanh(x-x0 - t)/(2.K)
        bParam=materialsParameters(1)  !nu
        x0Param=analyticParameters(1)  !x0
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=1.0_DP - TANH((x(1)-x0Param-time)/(2.0_DP*bParam))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
          & " is invalid for a Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
      !a.del u/del t + b.del^2 u/del x^2 + c.u.del u/del x = 0
      aParam=materialsParameters(1)
      bParam=materialsParameters(2)
      cParam=materialsParameters(3)
      SELECT CASE(analyticFunctionType)
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1)
        !Analytic solution is u(x,t)=(d+a.x)/(e+c.t)
        dParam = analyticParameters(1)
        eParam = analyticParameters(2)
        SELECT CASE(variableType)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=(dParam+aParam*x(1))/(eParam+cParam*time)
          CASE(GLOBAL_DERIV_S1)
            dependentValue=dParam/(eParam+cParam*time)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            dependentValue=0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
        !Analytic_solution=a.d+2.b/c(x-c.d.t+e)
        dParam = analyticParameters(1)
        eParam = analyticParameters(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=aParam*dParam+2.0_DP*bParam/(cParam*(x(1)-cParam*dParam*time+eParam))
          CASE(GLOBAL_DERIV_S1)
            dependentValue=-2.0_DP*bParam/(cParam*(x(1)-cParam*dParam*time+eParam)**2)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(globalDerivative)
          CASE(NO_GLOBAL_DERIV)
            dependentValue=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            dependentValue=0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NumberToVString(globalDerivative,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NumberToVString(variableType,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))// &
          & " is invalid for a generalised Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(equationsSubtype,"*",err,error))//" is invalid."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_AnalyticFunctionsEvaluate")
    RETURN
999 ERRORSEXITS("Burgers_AnalyticFunctionsEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE Burgers_AnalyticFunctionsEvaluate

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a burgers equation type of an fluid mechanics equations set class.
  SUBROUTINE Burgers_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_EquationsSetSolutionMethodSet",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
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
        & " is not valid for a Burgers equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Burgers_EquationsSetSolutionMethodSet",err,error)
    EXITS("Burgers_EquationsSetSolutionMethodSet")
    RETURN 1

  END SUBROUTINE Burgers_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Burgers type of a fluid mechanics equations set.
  SUBROUTINE Burgers_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
        & " is not valid for a Burgers type of fluid mechanics equations set."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    !Set full specification
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_BURGERS_EQUATION_TYPE,subtype]

    EXITS("Burgers_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Burgers_EquationsSetSpecificationSet",err,error)
    EXITS("Burgers_EquationsSetSpecificationSet")
    RETURN 1

  END SUBROUTINE Burgers_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Burgers equation type of a fluid mechanics equations set class.
  SUBROUTINE Burgers_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup a Burgers equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,esSpecification(3),geometricMeshComponent,geometricScalingType,lumpingType, &
      & numberOfAnalyticComponents,numberOfDependentComponents,numberOfDimensions,numberOfMaterialsComponents, &
      & solutionMethod,sparsityType
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
      EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is invalid for a Burgers equation subtype."
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
        CALL Burgers_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        !Do nothing
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Burgers equation."
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
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
          CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
            & err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
            & err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            numberOfDependentComponents=numberOfDimensions
          ELSE
            numberOfDependentComponents=1
          ENDIF
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
            & numberOfDependentComponents,err,error,*999)
          !Default to the geometric interpolation setup
          DO componentIdx=1,numberOfDependentComponents
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,componentIdx,geometricMeshComponent, &
              & err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
            CALL Field_ComponentMeshComponentSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & componentIdx,geometricMeshComponent,err,error,*999)
          ENDDO !componentIdx
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
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
            localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE, &
            & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            numberOfDependentComponents=numberOfDimensions
          ELSE
            numberOfDependentComponents=1
          ENDIF
          IF(numberOfDependentComponents==1) THEN
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
          ELSE
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
          ENDIF
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,numberOfDependentComponents, &
            & err,error,*999)
          CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            DO componentIdx=1,numberOfDependentComponents
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, &
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
          & " is invalid for a nonlinear Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
      !-----------------------------------------------------------------
      ! M a t e r i a l s   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(equationsMaterials)
      CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(esSpecification(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
          !Not an inviscid Burgers equation
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
            CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_DecompositionGetgeometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
              !1 materials field component i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
              numberOfMaterialsComponents=1
            CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
              !3 materials field components i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
              numberOfMaterialsComponents=3
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            !Set the number of materials components
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & numberOfMaterialsComponents,err,error,*999)
            !Default the materials components to the 1st geometric component interpolation setup with constant interpolation
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfMaterialsComponents
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
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
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
              !1 materials field component i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
              numberOfMaterialsComponents=1
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
              !3 materials field components i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
              numberOfMaterialsComponents=3
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfMaterialsComponents, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(esSpecification(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
          !Not an inviscid Burgers equation
          IF(equationsMaterials%materialsFieldAutoCreated) THEN
            !Finish creating the materials field
            CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
            !Set the default values for the materials field
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
              !1 materials field component. Default to du/dt - d^2u/dx^2 + u*(du/dx) = 0
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 1,-1.0_DP,err,error,*999)
            CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
              !3 materials field components. Default to du/dt - d^2u/dx^2 + u*(du/dx) = 0
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 1,1.0_DP,err,error,*999)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 2,-1.0_DP,err,error,*999)              
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & 3,1.0_DP,err,error,*999)
            CASE DEFAULT
              localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
                & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a nonlinear Burgers equation."
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
          & " is invalid for a nonlinear Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c   t y p e
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      CALL EquationsSet_AnalyticFunctionType(equationsSet,analyticFunctionType,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
          SELECT CASE(analyticFunctionType)
          CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
            !Check that domain is 1D
            IF(numberOfDimensions/=1) THEN
              localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " requires that there be 1 geometric dimension."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Check the materials values are constant
            CALL Field_ComponentInterpolationCheck(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE,1, &
              & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1
            numberOfAnalyticComponents=1
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))//" is invalid for a Burgers equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
          SELECT CASE(analyticFunctionType)
          CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1, &
            & EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
            !Check that domain is 1D
            IF(numberOfDimensions/=1) THEN
              localError="The number of geometric dimensions of "//TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " requires that there be 1 geometric dimension."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Check the materials values are constant
            CALL Field_ComponentInterpolationCheck(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & 2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            CALL Field_ComponentInterpolationCheck(equationsMaterials%materialsField,FIELD_U_VARIABLE_TYPE, &
              & 3,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            !Set analytic function type
            equationsSet%analytic%analyticFunctionType=equationsSetSetup%analyticFunctionType
            numberOfAnalyticComponents=2
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a generalised Burgers equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid for an analytical Burgers equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Create analytic field if required
        IF(numberOfAnalyticComponents>=1) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            !Create the auto created source field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsAnalytic%analyticField,err,error,*999)
            CALL Field_LabelSet(equationsAnalytic%analyticField,"Analytic Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsAnalytic%analyticField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsAnalytic%analyticField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsAnalytic%analyticField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsAnalytic%analyticField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsAnalytic%analyticField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsAnalytic%analyticField,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_VariableLabelSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,"Analytic",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            !Set the number of analytic components
            CALL Field_NumberOfComponentsSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
              & numberOfAnalyticComponents,err,error,*999)
            !Default the analytic components to the 1st geometric interpolation setup with constant interpolation
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfAnalyticComponents
              CALL Field_ComponentMeshComponentSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricMeshComponent,err,error,*999)
              CALL Field_ComponentInterpolationSet(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsAnalytic%analyticField,geometricScalingType,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,1,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            IF(numberOfAnalyticComponents==1) THEN
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            ELSE
              CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            ENDIF
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfAnalyticComponents, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        NULLIFY(analyticField)
        CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
        IF(ASSOCIATED(analyticField)) THEN
          IF(equationsAnalytic%analyticFieldAutoCreated) THEN
            !Finish creating the analytic field
            CALL Field_CreateFinish(equationsAnalytic%analyticField,err,error,*999)
            !Set the default values for the analytic field
            SELECT CASE(esSpecification(3))
            CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
              !Default the analytic parameter value to 0.0
              CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
            CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
              !Default the analytic parameter values to 1.0
              CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
              CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
            CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                & " is invalid for an analytical Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------
      CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Finish the equations
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            CALL Equations_CreateFinish(equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            NULLIFY(vectorMapping)
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.FALSE.,err,error,*999)
            ELSE
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
            ENDIF
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
            CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            NULLIFY(vectorMatrices)
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            !Set up matrix storage and structure
            CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
              !Set up lumping
              IF(esSpecification(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_LUMPED],err,error,*999)
              ELSE
                CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_UNLUMPED, &
                  & EQUATIONS_MATRIX_LUMPED],err,error,*999)
              ENDIF
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                IF(esSpecfication(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_DIAGONAL_STRUCTURE], &
                    & err,error,*999)
                ELSE
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                    & EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                ENDIF
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                IF(esSpecification(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_DIAGONAL_STRUCTURE], &
                    & err,error,*999)
                ELSE
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                ENDIF
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1, &
                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//RIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                IF(esSpecification(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                ELSE
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                    & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                ENDIF
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                IF(esSpecification(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                    & err,error,*999)
                ELSE
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                ENDIF
                CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1, &
                  & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "// TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
            CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            !Use the analytic Jacobian calculation
            CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
              & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
            & " is invalid for a Burgers equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
          NULLIFY(equations)
          CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
          CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          SELECT CASE(solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Finish the equations
            NULLIFY(equations)
            CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
            CALL Equations_CreateFinish(equations,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            NULLIFY(vectorMapping)
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_U_VARIABLE_TYPE,vectorMapping,err,error,*999)
            CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
            CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,[FIELD_U_VARIABLE_TYPE],err,error,*999)
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
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
              CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CALL EquationsMatricesVector_NonlinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                & err,error,*999)
              CALL EquationsMatricesVector_NonlinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,err,error,*999)
            !Use the analytic Jacobian calculation
            CALL EquationsMatricesVector_JacobianCalculationTypeSet(vectorMatrices,FIELD_U_VARIABLE_TYPE,1, &
              & EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED,err,error,*999)
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
            & " is invalid for a Burgers equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is invalid for an analytical Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a Burgers equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Burgers_EquationsSetSetup",err,error)
    RETURN 1

  END SUBROUTINE Burgers_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets up the Burgers problem pre-solve.
  SUBROUTINE Burgers_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: pSpecification(3)
    LOGICAL :: solverInitialised
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      CALL Solver_DynamicSolverInitialisedGet(solver,solverInitialised,err,error,*999)
      IF(solverInitialised) CALL Burgers_PreSolveUpdateAnalyticValues(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Burgers equation type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_PreSolve")
    RETURN
999 ERRORSEXITS("Burgers_PreSolve",err,error)
    RETURN 1
    
  END SUBROUTINE Burgers_PreSolve


  !
  !================================================================================================================================
  !
  
  !>Updates the boundary conditions and source term to the required analytic values
  SUBROUTINE Burgers_PreSolveUpdateAnalyticValues(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionCheckVariable,componentIdx,derivativeIdx,dimensionIdx, &
      & equationsSetIdx,globalDerivativeIndex,localDOFIdx,nodeIdx,numberOfComponents,numberOfDimensions,numberOfEquationsSets, &
      & numberOfNodes,numberOfNodeDerivatives,numberOfVariables,pSpecification(3),variableIdx,variableType
    REAL(DP) :: currentTime,dependentValue,normal(3),tangents(3,3),timeIncrement,x(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    TYPE(BoundaryConditionType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: analyticVariable,dependentVariable,geometricVariable,materialsVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations 
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_PreSolveUpdateAnalyticValues",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE, &
      & PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(boundaryConditions)
      CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
      !Loop over all the equation sets and set the appropriate field variable type BCs and
      !the source field associated with each equation set
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
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
          NULLIFY(analyticField)
          CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
          NULLIFY(analyticVariable)
          NULLIFY(analyticParameters)
          IF(ASSOCIATED(analyticField)) THEN
            CALL Field_VariableGet(analyticField,FIELD_U_VARIABLE_TYPE,analyticVariable,err,error,*999)
            CALL FieldVariable_ParameterSetDataGet(analyticVariable,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
          ENDIF
          NULLIFY(materialsField)
          CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
          NULLIFY(materialsVariable)
          NULLIFY(materialsParameters)
          IF(ASSOCIATED(materialsField)) THEN
            CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,materialsVariable,err,error,*999)
            CALL FieldVariable_ParameterSetDataGet(materialsVariable,FIELD_VALUES_SET_TYPE,materialsParameters,err,error,*999)
          ENDIF
          CALL Field_NumberOfVariablesGet(dependentField,numberOfVariables,err,error,*999)
          DO variableIdx=1,numberOfVariables
            NULLIFY(dependentVariable)
            CALL Field_VariableIndexGet(dependentField,variableIdx,dependentVariable,variableType,err,error,*999)
            CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
            DO componentIdx=1,numberOfComponents
              CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
                & err,error,*999)
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
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDOFIdx,err,error,*999)
                  x(dimensionIdx)=geometricParameters(localDOFIdx)
                ENDDO !dimensionIdx
                !Loop over the derivatives
                CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                DO derivativeIdx=1,numberOfNodeDerivatives
                  CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
                  CALL Burgers_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,currentTime, &
                    & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,dependentValue, &
                    & err,error,*999)
                  !Default to version 1 of each node derivative
                  CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                    & err,error,*999)
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDOFIdx, &
                    & dependentValue,err,error,*999)
                  NULLIFY(boundaryConditionsVariable)
                  CALL BoundaryConditions_VariableGet(boundaryConditions,dependentVariable,boundaryConditionsVariable, &
                    & err,error,*999)
                  boundaryConditionCheckVariable=boundaryConditionsVariable%conditionTypes(localDOFIdx)
                  IF(boundaryConditionCheckVariable==BOUNDARY_CONDITION_FIXED) THEN
                    CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_VALUES_SET_TYPE,localDOFIdx, &
                      & dependentValue,err,error,*999)
                  ENDIF
                ENDDO !derivativeIdx
              ENDDO !nodeIdx
            ENDDO !componentIdx
            CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
            CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)            
          ENDDO !variableIdx
          CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Burgers equation type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORS("Burgers_PreSolveUpdateAnalyticValues",err,error)
    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN 1

  END SUBROUTINE Burgers_PreSolveUpdateAnalyticValues


  !
  !================================================================================================================================
  !
  SUBROUTINE Burgers_PreSolveStoreCurrentSolution(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_PreSolveStoreCurrentSolution",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
      !Do nothing ???
    CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      !Do nothing ???
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Burgers equation type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_PreSolveStoreCurrentSolution")
    RETURN
999 ERRORS("Burgers_PreSolveStoreCurrentSolution",err,error)
    EXITS("Burgers_PreSolveStoreCurrentSolution")
    RETURN 1

  END SUBROUTINE Burgers_PreSolveStoreCurrentSolution

  !
  !================================================================================================================================
  !

  !>Sets up the Burgers problem post solve.
  SUBROUTINE Burgers_PostSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_PostSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE, &
      & PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      CALL Burgers_PostSolveOutputData(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a BURGERS type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_PostSolve")
    RETURN
999 ERRORSEXITS("Burgers_PostSolve",err,error)
    RETURN 1

  END SUBROUTINE Burgers_PostSolve

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE Burgers_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,numberOfEquationsSets,outputIteration,outputType,pSpecification(3)
    REAL(DP) :: currentTime,timeIncrement,startTime,stopTime
    CHARACTER(20) :: file,outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField
    TYPE(FieldsType), POINTER :: fields
    TYPE(RegionType), POINTER :: region
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError,method,filename

    ENTERS("Burgers_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
    CALL SYSTEM('mkdir -p ./output')
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
      CALL ControlLoop_TimesGet(controlLoop,startTime,stopTime,currentTime,timeIncrement,currentIteration,outputIteration, &
        & err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(region)
        CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
        NULLIFY(fields)
        CALL Region_FieldsGet(region,fields,err,error,*999)
        filename="./output/"//"StaticSolution"
        method="FORTRAN"
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
          CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
          CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
        ENDIF
        CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
        CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
      ENDDO !equationsSetIdx
    CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      CALL ControlLoop_TimesGet(controlLoop,startTime,stopTime,currentTime,timeIncrement,currentIteration,outputIteration, &
        & err,error,*999)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !Make sure the equations sets are up to date
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        IF(outputIteration/=0) THEN
          IF(currentTime<=stopTime) THEN
            IF(currentIteration<10) THEN
              WRITE(outputFile,'("TimeStep_000",I0)') currentIteration
            ELSE IF(currentIteration<100) THEN
              WRITE(outputFile,'("TimeStep_00",I0)') currentIteration
            ELSE IF(currentIteration<1000) THEN
              WRITE(outputFile,'("TimeStep_0",I0)') currentIteration
            ELSE IF(currentIteration<10000) THEN
              WRITE(outputFile,'("TimeStep_",I0)') currentIteration
            END IF
            file=outputFile
            filename="./output/"//"MainTime_"//TRIM(NumberToVString(currentIteration,"*",err,error))
            method="FORTRAN"
            IF(MOD(currentIteration,outputIteration)==0)  THEN
              NULLIFY(region)
              CALL EquationsSet_RegionGet(equationsSet,region,err,error,*999)
              NULLIFY(fields)
              CALL Region_FieldsGet(region,fields,err,error,*999)
              IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
              ENDIF
              CALL FIELD_IO_NODES_EXPORT(fields,filename,method,err,error,*999)
              CALL FIELD_IO_ELEMENTS_EXPORT(fields,filename,method,err,error,*999)
              IF(outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                CALL WriteString(GENERAL_OUTPUT_TYPE,filename,err,error,*999)
                CALL WriteString(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
              ENDIF
            ENDIF
            NULLIFY(equationsAnalytic)
            CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
            IF(ASSOCIATED(equationsAnalytic)) THEN
              NULLIFY(dependentField)
              CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
              CALL AnalyticAnalysis_Output(dependentField,file,err,error,*999)
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a Burgers equation type of a fluid mechanics problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Burgers_PostSolveOutputData",err,error)
    RETURN 1

  END SUBROUTINE Burgers_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Burgers problem.
  SUBROUTINE Burgers_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_ProblemSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(problem)) CALL FlagError("Problem is not associated.",err,error,*999)
    IF(ALLOCATED(problem%specification)) CALL FlagError("Problem specification is already allocated.",err,error,*999)
    IF(SIZE(problemSpecification,1)<3) THEN
      localError="The size of the specified problem specification array of "// &
        & TRIM(NumberToVString(SIZE(problemSpecification,1),"*",err,error)//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    ENDIF
    
    problemSubtype=problemSpecification(3)
    SELECT CASE(problemSubtype)
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE, &
      & PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a Burgers type of a fluid mechanics problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_BURGERS_EQUATION_TYPE,problemSubtype]

    EXITS("Burgers_ProblemSpecificationSet")
    RETURN
999 ERRORS("Burgers_ProblemSpecificationSet",err,error)
    EXITS("Burgers_ProblemSpecificationSet")
    RETURN 1

  END SUBROUTINE Burgers_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the Burgers problem.
  SUBROUTINE Burgers_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem set to setup a Burgers problem on.
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
    
    ENTERS("Burgers_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_STATIC_BURGERS_SUBTYPE, &
      & PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is invalid for a Burgers problem subtype."
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
          & " is invalid for a Burgers problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a simple control loop
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        IF(pSpecification(3)==PROBLEM_DYNAMIC_BURGERS_SUBTYPE) THEN
          CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
        ENDIF
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
          & " is invalid for a Burgers problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVERS_TYPE)
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
        IF(pSpecification(3)==PROBLEM_DYNAMIC_BURGERS_SUBTYPE) THEN
            !Set the solver to be a dynamic nonlinear solver
          CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Nonlinear dynamic solver",err,error,*999)
          CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
          CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
          !Set solver defaults
          CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
          CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
        ELSE
          !Set the solver to be a static nonlinear solver
          CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
          CALL Solver_TypeSet(solver,SOLVER_NONLINEAR_TYPE,err,error,*999)
          CALL Solver_LabelSet(solver,"Nonlinear solver",err,error,*999)
          CALL Solver_LibraryTypeSet(solver,SOLVER_PETSC_LIBRARY,err,error,*999)
        ENDIF
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a Burgers problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
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
        CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
        IF(pSpecification(3)==PROBLEM_DYNAMIC_BURGERS_SUBTYPE) THEN
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
         ELSE
          CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_STATIC,err,error,*999)
        ENDIF
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
          & " is invalid for a Burgers problem."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a Burgers problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Burgers_ProblemSetup")
    RETURN
999 ERRORSEXITS("Burgers_ProblemSetup",err,error)
    RETURN 1

  END SUBROUTINE Burgers_ProblemSetup

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a Burgers equation finite element equations set.
  SUBROUTINE Burgers_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error   !<The error string
    !Local Variables 
    INTEGER(INTG) columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,gaussPointIdx,meshComponent1, &
      & meshComponent2,rowComponentIdx,rowElementDOFIdx,componentIdx,rowElementParameterIdxxiIdx,variableType     
    REAL(DP) :: cParam,colsdPhidXi,colsPhi,dependentValue,dXidX,jacobian,jacobianGaussWeight,gaussWeight,sum1,sum2,sum3, &
      & rowsPhi,uValue,uDeriv
    LOGICAL :: evaluateJacobian,updateJacobianMatrix
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS1,DEPENDENT_BASIS2,geometricBasis
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme1,quadratureScheme2

    ENTERS("Burgers_FiniteElementJacobianEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      !OK
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
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) &
        & CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      
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
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType,dependentInterpParameters, &
        & err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,colsVariableType,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
      
      NULLIFY(materialsInterpParamters)
      NULLIFY(materialsInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      ENDIF
      cParam=1.0_DP
      
      !Loop over all Gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
            & err,error,*999)
          cParam=materialsInterpPoint%values(3,NO_PART_DERIV)
        ENDIF

        DO componentIdx=1,numberOfColsComponents
          uValue(componentIdx)=dependentInterpPoint%values(componentIdx,NO_PART_DERIV)
          dudX(componentIdx)=0.0_DP
          DO xiIdx=1,numberOfXi
            uDeriv(componentIdx,xiIdx)=dependentInterpPoint%values(componentIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx))
            dXidX(xiIdx,componentIdx)=geometricInterpPointMetrics%dXidX(xiIdx,columnComponentIdx)
            dudX(componentIdx)=uDeriv(componentIdx,xiIdx)*dXidX(xiIdx,componentIdx)
          ENDDO !xiIdx            
        ENDDO !componentIdx
        
        !Calculate jacobianGaussWeight.
        CALL FieldInterpolatedPointsMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(geometricQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
        !Loop over rows
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
              & rowsPhi,err,error,*999)
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
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,colElementParameterIdx,NO_PART_DERIV, &
                  & colsPhi,err,error,*999)
                !Loop over xi directions
                sum1=dudX(columnComponentIdx)*colsPhi
                sum2=0.0_DP
                IF(columnComponentIdx==rowComponentIdx) THEN
                  !Loop over spatial directions
                  DO componentIdx=1,geometricVariable%numberOfComponents
                    colsdPhidX=0.0_DP
                    DO xiIdx=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx1),gaussPointIdx,colsdPhidXi,err,error,*999)
                      colsdPhidX=colsdPhidX+colsdPhidXi*dXidX(xiIdx,componentIdx)
                    ENDDO !xiIdx
                    sum2=sum2+uValue(componentIdx)*colsdPhidX
                  ENDDO !componentIdx
                ENDIF
                matrixValue=cParam*(sum1+sum2)*rowsPhi
                jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                  & matrixValue*jacobianGaussWeight
              ENDDO !columnElementParameterIdx
            ENDDO !columnComponentIdx
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

    EXITS("Burgers_FiniteElementJacobianEvaluate")
    RETURN
999 ERRORS("Burgers_FiniteElementJacobianEvaluate",err,error)
    EXITS("Burgers_FiniteElementJacobianEvaluate")
    RETURN 1

  END SUBROUTINE Burgers_FiniteElementJacobianEvaluate

  !
  !================================================================================================================================
  !

  !>Evaluates the residual element stiffness matrices and RHS for a Burgers equation finite element equations set.
  SUBROUTINE Burgers_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) colsXiIdx,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,componentIdx,gaussPointIdx, &
      & meshComponent1,meshComponent2,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowsXiIdx,variableType
    REAL(DP) :: aParam,bParam,cParam,JGW,sum,sum1,rowsPhi,colsPhi,rowsdPhidXi(3),colsdPhidXi(3),uValue
    LOGICAL :: evaluateAny,evaluateDamping,evaluateLinearDynamic,evaluateResidual,evaluateRHS, &
      & evaluateStiffness,firstDamping,firstRHS,firstStiffness,updateStiffness,updateDamping,updateRHS,updateResidual
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS,geometricBasis,DEPENDENT_BASIS1,DEPENDENT_BASIS2
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: geometricVariable,dependentVariable
    TYPE(QuadratureSchemeType), POINTER :: quadratureScheme,quadratureScheme1,quadratureScheme2
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_FiniteElementResidualEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(esSpecification,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      !OK
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
    updateResidual=residualVector%updateVector
    NULLIFY(dynamicMapping)
    NULLIFY(dynamicMatrices)
    NULLIFY(linearMapping)
    NULLIFY(linearMatrices)
    NULLIFY(stiffnessMatrix)
    NULLIFY(dampingMatrix)
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
      CALL EquationsMapping_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMatrices_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
      updateStiffness=stiffnessMatrix%updateMatrix
      updateDamping=dampingMatrix%updateMatrix
    CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      CALL EquationsMappingVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)
      CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,1,stiffnessMatrix,err,error,*999)
      updateStiffness=stiffness%updateMatrix
      updateDamping=.FALSE.
    CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      CALL EquationsMapping_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
      CALL EquationsMatrices_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
      CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,dampingMatrix,err,error,*999)
      updateStiffness=.FALSE.
      updateDamping=dampingMatrix%updateDamping
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a BURGERS equation type of a fluid mechanics equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    NULLIFY(rhsVector)
    CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,1,rhsVector,err,error,*999)
    updateRHS=rhsVector%updateVector

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateResidual.OR.updateRHS)

    IF(update) THEN
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      NULLIFY(dependentField)
      CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
      NULLIFY(materialsField)
      IF(esSpecification(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) &
        & CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        
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
      NULLIFY(materialsInterpPoint)
      IF(esSpecification(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & materialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      ENDIF
      
      aParam=1.0_DP
      bParam=1.0_DP
      cParam=1.0_DP
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        IF(esSpecification(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
            & err,error,*999)
          IF(esSpecification(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
            aParam=materialsInterpPoint%values(1,NO_PART_DERIV)
            bParam=materialsInterpPoint%values(2,NO_PART_DERIV)
            cParam=materialsInterpPoint%values(3,NO_PART_DERIV)
          ELSE
            bParam=materialsInterpPoint%values(1,NO_PART_DERIV)
          ENDIF
        ENDIF
        IF(updateResidual) THEN
          ududX=0.0_DP
          DO columnComponentIdx=1,numberOfColsComponents
            uValue(columnComponentIdx)=dependentInterpPoint%values(columnComponentIdx,NO_PART_DERIV)
            dudX(columnComponentIdx)=0.0_DP
            DO colsXiIdx=1,numberOfXi
              dudXi(columnComponentIdx,colsXiIdx)=dependentInterpPoint%values(columnComponentIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(colsXiIdx))
              dudX(columnComponentIdx)=dudX(columnComponentIdx)+dudXi(columnComponentIdx,colsXiIdx)* &
                & geometricInterpPointMetrics%dXidX(colsXiIdx,columnComponentIdx)
            ENDDO !colsXiIdx
            ududX=ududX+uValue(columnComponentIdx)*dudX(columnComponentIdx)
          ENDDO !componentIdx
        ENDIF
        
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
            IF(updateMatrices) THEN
              DO rowsXiIdx=1,numberOfXi
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx, &
                  & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(rowsXiIdx),gaussPointIdx,rowsdPhidXi(rowsXiIdx),err,error,*999)
              ENDDO !rowsXiIdx
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
                  !Diffusion matrix
                  IF(updateStiffness) THEN
                    DO colsXiIdx=1,numberOfXi
                      CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(colsXiIdx),gaussPointIdx,colsdPhidXi(colsXiIdx), &
                        & err,error,*999)
                    ENDDO !colsXiIdx
                    sum=0.0_DP
                    !Calculate sum
                    DO rowsXiIdx=1,numberOfXi
                      DO colsXiIdx=1,numberOfXi
                        sum=sum+colsdPhidXi(colsXiIdx)*rowsdPhidXi(rowsXiIdx)*geometricInterpPointMetrics%gu(rowsXiIdx,colsXiIdx)
                      ENDDO !colsXiIdx
                    ENDDO !rowsXiIdx
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)- &
                      & bParam*sum*jacobianGaussWeight
                  ENDIF
                  !Mass matrix
                  IF(updateDamping) THEN
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                      & rowsPhi,err,error,*999)
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(colsQuadratureScheme,columnElementParameterIdx,NO_PART_DERIV, &
                      & colsPhi,err,error,*999)
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & aParam*rowsPhi*colsPhi*jacobianGaussWeight
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !Stiffness or Damping
            !Calculate RHS
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
            ENDIF
            !Calculate nonlinear vector
            IF(updateResidual) THEN
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowsQuadratureScheme,rowElementParameterIdx,NO_PART_DERIV, &
                & rowsPhi,err,error,*999)
              nonlinearMatrices%elementResidual%vector(rowElementDOFIdx)= &
                & nonlinearMatrices%elementResidual%vector(rowElementDOFIdx)+ & 
                & cParam*ududX*rowsPhi*jacobianGaussWeight
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
            IF(updateMatrices) THEN
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
                  IF(updateStiffness) THEN                  
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF
                  IF(updateDamping) THEN                  
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                      & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                      & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !update matrix
            IF(updateRHS) THEN
              rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
            IF(updateResidual) THEN
              residualVector%elementVector%vector(rowElementDOFIdx)=residualVector%elementVector%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling

    ENDIF !Update

    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Burgers_FiniteElementResidualEvaluate",err,error)
    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN 1

  END SUBROUTINE Burgers_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

END MODULE BurgersEquationsRoutines
