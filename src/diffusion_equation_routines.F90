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

!>This module handles all diffusion equation routines.
MODULE DiffusionEquationsRoutines

  USE AnalyticAnalysisRoutines
  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BoundaryConditionsRoutines
  USE BoundaryConditionAccessRoutines
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE CoordinateSystemRoutines
  USE DecompositionRoutines
  USE DecompositionAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE MeshAccessRoutines
  USE ProblemAccessRoutines
  USE RegionAccessRoutines
  USE Strings
  USE SolverRoutines
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

  PUBLIC Diffusion_AnalyticFunctionsEvaluate

  PUBLIC Diffusion_BoundaryConditionAnalyticCalculate

  PUBLIC Diffusion_EquationsSetSetup

  PUBLIC Diffusion_EquationsSetSolutionMethodSet
  
  PUBLIC Diffusion_EquationsSetSpecificationSet

  PUBLIC Diffusion_FiniteElementCalculate

  PUBLIC Diffusion_FiniteElementJacobianEvaluate
  
  PUBLIC Diffusion_FiniteElementResidualEvaluate
  
  PUBLIC Diffusion_ProblemSpecificationSet

  PUBLIC Diffusion_ProblemSetup

  PUBLIC Diffusion_PreSolve

  PUBLIC Diffusion_PostSolve

  PUBLIC Diffusion_PostLoop  

!!TODO: should the following two routines really be public???

  PUBLIC Diffusion_PreSolveGetSourceValue
  
  PUBLIC Diffusion_PreSolveStoreCurrentSolution


CONTAINS

  !
  !================================================================================================================================
  !
  
  !>Evaluate the analytic solutions for a diffusion equation
  SUBROUTINE Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
    & globalDerivativeIndex,componentNumber,analyticParameters,materialsParameters,analyticValue,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER, INTENT(IN) :: equationsSet !<The equations set to evaluate
    INTEGER(INTG), INTENT(IN) :: analyticFunctionType !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: x(:) !<x(dimensionIdx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: tangents(:,:) !<tangents(dimensionIdx,xiIdx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: normal(:) !<normal(dimensionIdx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: time !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: variableType !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: globalDerivativeIndex !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: componentNumber !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: analyticParameters(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: materialsParameters(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: analyticValue !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: equationsSubType,esSpecification(3)
    REAL(DP) :: k,phi,A,B,C,D,A1,A2,A3,A4,aParam,bParam,cParam,kParam,lParam,constParam,betaParam,lambdaParam,muParam
    TYPE(VARYING_STRING) :: localError

    !These are parameters for the analytical solution
    !k = 1.0_DP !this is a time decay constant for the exponential term
    !phi = 0.785398163397_DP !pi/4 - this sets the orientation of the solution relative to the axes
    k = 10.0_DP !this is a time decay constant for the exponential term
    phi = 1.0_DP !pi/4 - this sets the orientation of the solution relative to the axes
 
    !Solution parameters for 
    A1 = 0.4_DP
    A2 = 0.3_DP
    A3 = 0.2_DP
    A4 = 0.1_DP

    ENTERS("Diffusion_AnalyticFunctionsEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
     
    equationsSubType=esSpecification(3)
    SELECT CASE(equationsSubType)
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE)
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
            analyticValue=aParam*EXP(4.0_DP*PI**2*kParam*time/lParam**2)*COS(2.0_DP*PI*x(1)/lParam+bParam)+cParam
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
            ANALYTICVALUE=0.0_DP
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
            analyticValue=EXP(-k*time)*SIN((SQRT(k))*(x(1)*COS(phi)+x(2)*SIN(phi)))!Need to specify time, k and phi!
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding 
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
            analyticValue=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            ANALYTICVALUE=0.0_DP !set to zero currently- actual value for diffusion solution needs adding
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
            analyticValue=EXP(A*time)*EXP(B*x(1))*EXP(C*x(2))*EXP(D*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding 
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
            analyticValue=1.0_DP/(betaParam+constParam*EXP(lambdaParam*time+muParam*x(1)))**2
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
            analyticValue=0.0_DP                                 
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
            analyticValue=-2.0_DP/cParam*LOG(betaParam+constParam*EXP(muParam*X(1)-aParam*cParam*time/2.0_DP))
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
            analyticValue=0.0_DP 
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
    CASE(EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE)
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
            analyticValue=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2))
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
            analyticValue=0.0_DP!set to zero currently- actual value for diffusion solution needs adding                                    
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
            analyticValue=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2))
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
            analyticValue=0.0_DP!set to zero currently- actual value for diffusion solution needs adding
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
            analyticValue=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding 
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
            analyticValue=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding
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
            analyticValue=A1*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding 
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
            analyticValue=A2*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding
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
            analyticValue=A3*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP!set to zero currently- actual value for diffusion solution needs adding 
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
            analyticValue=A4*EXP(-1.0_DP*time)*(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
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
            analyticValue=0.0_DP !set to zero currently- actual value for diffusion solution needs adding
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
        localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NumberToVString(equationsSubType,"*",err,error))//" is invalid."
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

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Diffusion_BoundaryConditionAnalyticCalculate(equationsSet,boundaryConditions,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to calculate the boundary conditions for
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions !<A pointer to the boundary conditions to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,componentIdx,derivativeIdx,dimensionIdx,esSpecification(3),globalDerivativeIndex, &
      & imyMatrixIdx,localDofIdx,nodeIdx,numberOfComponents,numberOfDimensions,numberOfNodes,numberOfNodeDerivatives, &
      & numberOfVariables,numberOfVersions,variableIdx,variableType,versionIdx
    INTEGER(INTG), POINTER :: equationsSetParameters(:)
    REAL(DP) :: initialValue,normal(3),tangents(3,3),time,analyticValue,x(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    LOGICAL :: boundaryNode
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(FieldType), POINTER :: analyticField,dependentField,equationsSetField,geometricField,materialsField
    TYPE(FieldVariableType), POINTER :: dependentVariable,geometricVariable
 
    ENTERS("Diffusion_BoundaryConditionAnalyticCalculate",err,error,*999)
    
    IF(.NOT.ASSOCIATED(boundaryConditions)) CALL FlagError("Boundary conditions is not associated.",err,error,*999)
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    CALL EquationsSet_AssertAnalyticIsCreated(equationsSet,err,error,*999)
    CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
    
    NULLIFY(dependentField)
    CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
    NULLIFY(geometricField)
    CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
    NULLIFY(geometricVariable)
    CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
    CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
    NULLIFY(geometricParameters)
    CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
    NULLIFY(analyticField)
    NULLIFY(analyticParameters)    
    CALL EquationsSet_AnalyticFieldExists(equationsSet,analyticField,err,error,*999)
    IF(ASSOCIATED(analyticField)) &
      & CALL Field_ParameterSetDataGet(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,analyticParameters,err,error,*999)
    NULLIFY(materialsField)
    NULLIFY(materialsParameters)
    CALL EquationsSet_MaterialsFieldExists(equationsSet,materialsField,err,error,*999)
    IF(ASSOCIATED(materialsField)) &
      CALL Field_ParameterSetDataGet(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,materialsParameters, &
      & err,error,*999)           
    CALL EquationsSet_AnalyticTimeGet(equationsSet,time,err,error,*999)
    IF(equationsSet%specification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)THEN
      !If a multi-comp model, we will use the equations set field information to assign only the appropriate
      !field variable boundary conditions. Use predetermined mapping from equations set field compartment number
      !to field variable type
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
        CALL FieldVariable_ParameterSetEnsureCreated(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_NumberOfComponentsGet(dependentVariable,numberOfComponents,err,error,*999)
        DO componentIdx=1,numberOfComponents
          CALL FieldVariable_ComponentInterpolationCheck(dependentVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
            & err,error,*999)
          NULLIFY(domain)
          CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
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
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,dimensionIdx,localDofIdx,err,error,*999)
              x(dimensionIdx)=geometricParameters(localDofIdx)
            ENDDO !dimensionIdx
            CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            !Loop over the derivatives
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time,variableType, &
                & globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,analyticValue,err,error,*999)
              !Default to version 1 of each node derivative
              CALL FieldVariable_LocalNodeDOFGet(dependentVariable,1,derivativeIdx,nodeIdx,dimensionIdx,localDofIdx,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateLocalDOF(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDofIdx, &
                & analyticValue,err,error,*999)
              IF(MOD(variableType,FIELD_NUMBER_OF_VARIABLE_SUBTYPES)==FIELD_U_VARIABLE_TYPE) THEN
                IF(boundaryNode) THEN
                  !If we are a boundary node then set the analytic value on the boundary
                  CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDofIdx,BOUNDARY_CONDITION_FIXED, &
                    & analyticValue,err,error,*999)
                ELSE
                  CALL FieldVariable_ParameterSetUpdateLocalDof(dependentVariable,FIELD_VALUES_SET_TYPE,localDofIdx, &
                    & analyticValue,err,error,*999)
                ENDIF
              ENDIF
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !component_idx
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
    ELSE
      !for single physics diffusion problems use standard analytic calculate
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
          CALL FieldVariable_DomainGet(dependentVariable,componentIdx,domain,err,error,*999)
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
              CALL FieldVariable_LocalNodeDOFGet(geometricVariable,1,1,nodeIdx,componentIdx,localDofIdx,err,error,*999)
              x(dimensionIdx)=geometricParameters(localDofIdx)
            ENDDO !dimensionIdx
            CALL DomainNodes_NodeBoundaryNodeGet(domainNodes,nodeIdx,boundaryNode,err,error,*999)
            !Loop over the derivatives
            CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
            DO derivativeIdx=1,numberOfNodeDerivatives
              CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
              CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,0.0_DP, &
                & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,initialValue, &
                & err,error,*999)
              CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,time, &
                & variableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters,analyticValue, &
                & err,error,*999)
              CALL DomainNodes_DerivativeNumberOfVersionsGet(domainNodes,derivativeIdx,nodeIdx,numberOfVersions,err,error,*999)
              DO versionIdx=1,numberOfVersions
                CALL FieldVariable_LocalNodeDOFGet(dependentVariable,versionIdx,derivativeIdx,nodeIdx,componentIdx,localDofIdx, &
                  & err,error,*999)
                CALL Field_ParameterSetUpdateLocalDof(dependentField,variableType,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                  & localDofIdx,analyticValue,err,error,*999)
                IF(variableType==FIELD_U_VARIABLE_TYPE) THEN
                  IF(boundaryNode) THEN
                    !If we are a boundary node then set the analytic value on the boundary
                    CALL BoundaryConditions_SetLocalDOF(boundaryConditions,dependentVariable,localDofIdx, &
                      & BOUNDARY_CONDITION_FIXED,analyticValue,err,error,*999)
                  ELSE
                    !Set the initial condition.
                    CALL FieldVariable_ParameterSetUpdateLocalDof(dependentVariable,FIELD_VALUES_SET_TYPE,localDofIdx, &
                      & initialValue,err,error,*999)
                  ENDIF
                ENDIF
              ENDDO !versionIdx
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        ENDDO !component_idx
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateStart(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(dependentVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
      ENDDO !variableIdx
    ENDIF
    IF(ASSOCIATED(materialsField)) &
      & CALL Field_ParameterSetDataRestore(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,materialsParameters, &
      & err,error,*999)
    IF(ASSOCIATED(analyticField)) &
      & CALL Field_ParameterSetDataRestore(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,analyticParameters, &
      & err,error,*999)            
    CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)            

    EXITS("Diffusion_BoundaryConditionAnalyticCalculate")
    RETURN
999 ERRORS("Diffusion_BoundaryConditionAnalyticCalculate",err,error)
    EXITS("Diffusion_BoundaryConditionAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Diffusion_BoundaryConditionAnalyticCalculate

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a diffusion equation type of an classical field equations set class.
  SUBROUTINE Diffusion_EquationsSetSolutionMethodSet(equationsSet,solutionMethod,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: solutionMethod !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: esSpecification(3)
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Diffusion_EquationsSetSolutionMethodSet",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, & 
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) 
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
        & " is not valid for a diffusion equation type of an classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
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
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_EquationsSetSpecificationSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(ALLOCATED(equationsSet%specification)) CALL FlagError("Equations set specification is already allocated.",err,error,*999)    
    IF(SIZE(specification,1)<3) THEN
      localError="The size of the specified specification array of "// &
        & TRIM(NumberToVString(SIZE(specification,1),"*",err,error))//" is invalid. The size should be >= 3."
      CALL FlagError(localError,err,error,*999)
    END IF
    
    subtype=specification(3)
    SELECT CASE(subtype)
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
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
    ALLOCATE(equationsSet%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
    equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_DIFFUSION_EQUATION_TYPE,subtype]
 
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
  SUBROUTINE Diffusion_EquationsSetSetup(equationsSet,equationsSetSetup,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: equationsSetSetup !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: compartmentIdx,compartmentCount,componentIdx,esSpecification(3),geometricComponentNumber, &
      & geometricMeshComponent,geometricScalingType,lumpingType,materialVariableTypes(FIELD_NUMBER_OF_VARIABLE_TYPES), &
      & myMatrixNumber,numberOfAnalyticComponents,numberOfCompartments,numberOfDependentComponents,numberOfDependentVariables, &
      & numberOfDimensions,numberOfEquationsSetComponents,numberOfEquationsSetVariables, &
      & numberOfMaterialCoefficients(FIELD_NUMBER_OF_VARIABLE_TYPES),numberOfMaterialComponents(FIELD_NUMBER_OF_VARIABLE_TYPES), &
      & numberOfMaterialsCouplingComponents,numberOfMaterialVariables,numberOfSourceComponents,numberOfSourceVariables, &
      & solutionMethod,sparsityType,variableIdx
    INTEGER(INTG), POINTER :: equationsSetFieldData(:)
    INTEGER(INTG), ALLOCATABLE :: variableTypes(:),variableUTypes(:),couplingMatrixStorageType(:),couplingMatrixStructureType(:)
    REAL(DP) :: aParam,bParam,cParam
    TYPE(DecompositionType), POINTER :: geometricDecomposition
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsSetEquationsFieldType), POINTER :: equationsField
    TYPE(EquationsSetMaterialsType), POINTER :: equationsMaterials
    TYPE(EquationsSetSourceType), POINTER :: equationsSource
    TYPE(FieldType), POINTER :: analyticField,dependentField,equationsSetField,geometricField,materialsField,sourceField
    TYPE(RegionType), POINTER :: region
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_EquationsSetSetup",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
      !OK
    CASE DEFAULT
      localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion type of a classical field equations set."
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
        CALL Diffusion_EquationsSetSolutionMethodSet(equationsSet,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
        IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          numberOfEquationsSetVariables = 1
          numberOfEquationsSetComponents = 2
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            !Create the auto created equations set field
            CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsField%equationsSetField, &
              & err,error,*999)
            CALL Field_LabelSet(equationsField%equationsSetField,"Equations Set Field",err,error,*999)
            CALL Field_TypeSetAndLock(equationsField%equationsSetField,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(equationsField%equationsSetField,FIELD_INDEPENDENT_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfVariablesSet(equationsField%equationsSetField,numberOfEquationsSetVariables, &
              & err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsField%equationsSetField,[FIELD_U_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_VariableLabelSet(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,"Equations", &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE, &
              & numberOfEquationsSetComponents,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,numberOfEquationsSetVariables,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfEquationsSetComponents, &
              & err,error,*999)
          ENDIF
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE.OR. &
          & esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE) THEN
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            CALL Field_CreateFinish(equationsField%equationsSetField,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1, &
              & 1_INTG,err,error,*999)
            CALL Field_ComponentValuesInitialise(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2, &
              & 1_INTG,err,error,*999)
          ENDIF
        ENDIF
!!TODO: Check valid setup
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
      !-----------------------------------------------------------------
      ! G e o m e t r i c   f i e l d
      !-----------------------------------------------------------------
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
        !do nothing 
      CASE(EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL Field_ParameterSetEnsureCreated(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_DISPLACEMENT_SET_TYPE, &
            & err,error,*999)
          CALL Field_ParameterSetEnsureCreated(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
            & err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          !Do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a linear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
        SELECT CASE(equationsSetSetup%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          NULLIFY(equationsField)
          CALL EquationsSet_EquationsFieldGet(equationsSet,equationsField,err,error,*999)
          NULLIFY(geometricField)
          CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
          numberOfEquationsSetComponents = 2
          IF(equationsField%equationsSetFieldAutoCreated) THEN
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsField%equationsSetField,geometricDecomposition, &
              & err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsField%equationsSetField,geometricField,err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricComponentNumber,err,error,*999)
            DO componentIdx=1,numberOfEquationsSetComponents
              CALL Field_ComponentMeshComponentSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & geometricComponentNumber,err,error,*999)
              CALL Field_ComponentInterpolationSetAndLock(equationsField%equationsSetField,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !componentIdx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
            CALL Field_ScalingTypeSet(equationsField%equationsSetField,geometricScalingType,err,error,*999)
          ELSE
            !Do nothing
          ENDIF
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          !Do nothing
        CASE DEFAULT
          localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
            & " is invalid for a linear diffusion equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      END SELECT
    CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
      !-----------------------------------------------------------------
      ! D e p e n d e n t   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
!!TODO: setup the coupled field
            CALL FlagError("Not implemented.",err,error,*999)
          ELSE
            !Check the field created by advection-diffusion routines for the coupled problem
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,4,err,error,*999)
            CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE,1,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              componentIdx=1
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_V_VARIABLE_TYPE,componentIdx, &
                & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
              CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,FIELD_DELVDELN_VARIABLE_TYPE, & 
                & componentIdx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
!!TODO: setup field
            CALL FlagError("Not implemented.",err,error,*999)
          ELSE
            !uses number of compartments to check that appropriate number and type of variables have been set on the
            !dependent field
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_DEPENDENT_TYPE,err,error,*999)                 
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
              & err,error,*999)
            numberOfCompartments=equationsSetFieldData(2)
            CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,2*numberOfCompartments,err,error,*999)
            !Create & populate array storing all of the relevant variable types against which to check the field variables
            ALLOCATE(variableTypes(2*numberOfCompartments),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable types.",err,error,*999)
            DO compartmentIdx=1,numberOfCompartments
              variableTypes(2*compartmentIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(compartmentIdx-1))
              variableTypes(2*compartmentIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(compartmentIdx-1))
            ENDDO
            CALL Field_VariableTypesCheck(equationsSetSetup%field,variableTypes,err,error,*999)            
            DO compartmentIdx=1,2*numberOfCompartments
              CALL Field_DimensionCheck(equationsSetSetup%field,variableTypes(compartmentIdx),FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(equationsSetSetup%field,variableTypes(compartmentIdx),FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,variableTypes(compartmentIdx),1,err,error,*999)
            ENDDO !compartmentIdx
            CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              componentIdx=1
              DO compartmentIdx=1,2*numberOfCompartments
                CALL Field_ComponentInterpolationCheck(equationsSetSetup%field,variableTypes(compartmentIdx),componentIdx, &
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
              localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
            IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
          ENDIF
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE, &
          & EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          !Standard case
          numberOfDependentVariables=2
          numberOfDependentComponents=1
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
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
            CALL Field_VariableLabelSet(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
              & err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DimensionSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE, &
              & err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
              & numberOfDependentComponents,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsSet%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
              & numberOfDependentComponents,err,error,*999)
            !Default to the geometric interpolation setup
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            DO componentIdx=1,numberOfDependentComponents
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
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, & 
              & err,error,*999)
            CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfDimensions, &
              & err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_DELUDELN_VARIABLE_TYPE, & 
              & numberOfDimensions,err,error,*999)
            CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
            SELECT CASE(solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              DO componentIdx=1,numberOfDimensions
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
        END SELECT
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSet%dependent%dependentFieldAutoCreated) THEN
          CALL Field_CreateFinish(equationsSet%dependent%dependentField,err,error,*999)
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation"
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
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE)
        !i.e., a(x).\delby{u(x,t)}{t}+div(\sigma(x).grad(u(x,t)))+s(x)=0
        numberOfMaterialVariables=1
        materialVariableTypes(1)=FIELD_U_VARIABLE_TYPE
        numberOfMaterialCoefficients(1)=1 !a
        numberOfMaterialComponents(1)=numberOfMaterialCoefficients(1)+NUMBER_OF_VOIGT(numberOfDimensions) !a + sigma
      CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
        !i.e., a(x).\delby{u(x,t)}{t}+div(\sigma(x).grad(u(x,t)))+b(x)u(x,t)+s(x)=0
        numberOfMaterialVariables=1
        materialVariableTypes(1)=FIELD_U_VARIABLE_TYPE
        numberOfMaterialCoefficients(1)=2 !a + b
        numberOfMaterialComponents(1)=numberOfMaterialCoefficients(1)+NUMBER_OF_VOIGT(numberOfDimensions) !a + b + sigma
      CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE)
        !i.e., a(x).\delby{u(x,t)}{t}+div(\sigma(x).grad(u(x,t)))+b(x)u(x,t)+c(x).u^2(x,t)+s(x)=0
        numberOfMaterialVariables=1
        materialVariableTypes(1)=FIELD_U_VARIABLE_TYPE
        numberOfMaterialCoefficients(1)=3 !a + b + c
        numberOfMaterialComponents(1)=numberOfMaterialCoefficients(1)+NUMBER_OF_VOIGT(numberOfDimensions) !a + b + c + sigma
      CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
        !i.e., a(x).\delby{u(x,t)}{t}+div(\sigma(x).grad(u(x,t)))+b(x)e^[c(x).u(x,t)]+s(x)=0
        numberOfMaterialVariables=1
        materialVariableTypes(1)=FIELD_U_VARIABLE_TYPE
        numberOfMaterialCoefficients(1)=3 !a + b + c
        numberOfMaterialComponents(1)=numberOfMaterialCoefficients(1)+NUMBER_OF_VOIGT(numberOfDimensions) !a + b + c + sigma
      CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
        numberOfMaterialVariables=1
        materialVariableTypes(1)=FIELD_U_VARIABLE_TYPE
        numberOfMaterialCoefficients(1)=2 !??
        numberOfMaterialComponents(1)=numberOfMaterialCoefficients(1)+NUMBER_OF_VOIGT(numberOfDimensions) !?? + sigma
      CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
        numberOfMaterialVariables=2
        materialVariableTypes(1:2)=[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE]
        numberOfMaterialCoefficients(1)=0 !??      
        numberOfMaterialComponents(1)=numberOfDimensions !!OR NUMBER OF VOIGT NUMBER OF DIMENSIONS
        NULLIFY(equationsSetField)
        CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
        NULLIFY(equationsSetFieldData)
        CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
          & err,error,*999)
        numberOfCompartments=equationsSetFieldData(2)
        numberOfMaterialCoefficients(2)=numberOfCompartments !??      
        numberOfMaterialComponents(2)=numberOfCompartments
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(esSpecification(3),"*",err,error))// &
          & " is not valid for a diffusion type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Create the auto created materials field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsMaterials%materialsField,err,error,*999)
          CALL Field_LabelSet(equationsMaterials%materialsField,"Materials Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsMaterials%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsMaterials%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsMaterials%materialsField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsMaterials%materialsField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsMaterials%materialsField,numberOfMaterialVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsMaterials%materialsField, &
            & materialVariableTypes(1:numberOfMaterialVariables),err,error,*999)
          DO variableIdx=1,numberOfMaterialVariables
            CALL Field_VariableLabelSet(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
              & "Materials",err,error,*999)
            CALL Field_DimensionSetAndLock(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsSetAndLock(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
              & numberOfMaterialComponents(variableIdx),err,error,*999)
            CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
            !Default the materials components to the first geometric interpolation setup with constant interpolation
            DO componentIdx=1,numberOfMaterialComponents(variableIdx)
              CALL Field_ComponentInterpolationSet(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
                & componentIdx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL Field_ComponentMeshComponentSet(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
                & componentIdx,geometricMeshComponent,err,error,*999)
            ENDDO !componentIdx
          ENDDO !variableIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsMaterials%materialsField,geometricScalingType,err,error,*999) 
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_MATERIAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,numberOfMaterialVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,materialVariableTypes(1:numberOfMaterialVariables),err,error,*999)
          DO variableIdx=1,numberOfMaterialVariables
            CALL Field_DimensionCheck(equationsSetSetup%field,materialVariableTypes(variableIdx),FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(equationsSetSetup%field,materialVariableTypes(variableIdx),FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,materialVariableTypes(variableIdx), &
              & numberOfMaterialComponents(variableIdx),err,error,*999)
          ENDDO !variableIdx
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsMaterials%materialsFieldAutoCreated) THEN
          !Finish creating the materials field
          CALL Field_CreateFinish(equationsMaterials%materialsField,err,error,*999)
          !Set the default values of the materials values
          DO variableIdx=1,numberOfMaterialVariables
            !Set the coefficient values to 1.0
            DO componentIdx=1,numberOfMaterialCoefficients(variableIdx)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
            !Set the additional conductivity tensor to the identity tensor
            !First do the diagonal
            DO componentIdx=numberOfMaterialCoefficients(variableIdx)+1, &
              & MIN(numberOfMaterialCoefficients(variableIdx)+numberOfDimensions+1,numberOfMaterialComponents(variableIdx))
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
                & FIELD_VALUES_SET_TYPE,componentIdx,1.0_DP,err,error,*999)
            ENDDO !componentIdx
            !Now do the off-diagonal
            DO componentIdx=MIN(numberOfMaterialCoefficients(variableIdx)+numberOfDimensions+2, &
              & numberOfMaterialComponents(variableIdx)),numberOfMaterialComponents(variableIdx)
              CALL Field_ComponentValuesInitialise(equationsMaterials%materialsField,materialVariableTypes(variableIdx), &
                & FIELD_VALUES_SET_TYPE,componentIdx,0.0_DP,err,error,*999)
            ENDDO !componentIdx
          ENDDO !variableIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
      !-----------------------------------------------------------------
      ! S o u r c e   f i e l d
      !-----------------------------------------------------------------
      NULLIFY(equationsSource)
      CALL EquationsSet_SourceGet(equationsSet,equationsSource,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      numberOfSourceVariables=1
      numberOfSourceComponents=1
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          !Create the auto created source field
          CALL Field_CreateStart(equationsSetSetup%fieldUserNumber,region,equationsSource%sourceField,err,error,*999)
          CALL Field_LabelSet(equationsSource%sourceField,"Source Field",err,error,*999)
          CALL Field_TypeSetAndLock(equationsSource%sourceField,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeSetAndLock(equationsSource%sourceField,FIELD_INDEPENDENT_TYPE,err,error,*999)
          NULLIFY(geometricDecomposition)
          CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
          CALL Field_DecompositionSetAndLock(equationsSource%sourceField,geometricDecomposition,err,error,*999)
          CALL Field_GeometricFieldSetAndLock(equationsSource%sourceField,geometricField,err,error,*999)
          CALL Field_NumberOfVariablesSetAndLock(equationsSource%sourceField,numberOfSourceVariables,err,error,*999)
          CALL Field_VariableTypesSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,err,error,*999)
          CALL Field_VariableLabelSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,"Source",err,error,*999)
          CALL Field_DimensionSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE, &
            & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          !Set the number of source components
          CALL Field_NumberOfComponentsSetAndLock(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,numberOfSourceComponents, &
            & err,error,*999)
          !Default the source components to the first geometric interpolation setup with node based interpolation
          CALL Field_ComponentMeshComponentGet(geometricField,FIELD_U_VARIABLE_TYPE,1,geometricMeshComponent,err,error,*999)
          DO componentIdx=1,numberOfSourceComponents
            CALL Field_ComponentMeshComponentSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & geometricMeshComponent,err,error,*999)
            CALL Field_ComponentInterpolationSet(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,componentIdx, &
              & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
          ENDDO !componentIdx
          !Default the field scaling to that of the geometric field
          CALL Field_ScalingTypeGet(geometricField,geometricScalingType,err,error,*999)
          CALL Field_ScalingTypeSet(equationsSource%sourceField,geometricScalingType,err,error,*999)
        ELSE
          !Check the user specified field
          CALL Field_TypeCheck(equationsSetSetup%field,FIELD_GENERAL_TYPE,err,error,*999)
          CALL Field_DependentTypeCheck(equationsSetSetup%field,FIELD_INDEPENDENT_TYPE,err,error,*999)
          CALL Field_NumberOfVariablesCheck(equationsSetSetup%field,numberOfSourceVariables,err,error,*999)
          CALL Field_VariableTypesCheck(equationsSetSetup%field,[FIELD_U_VARIABLE_TYPE],err,error,*999)
          CALL Field_DimensionCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
          CALL Field_DataTypeCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
          CALL Field_NumberOfComponentsCheck(equationsSetSetup%field,FIELD_U_VARIABLE_TYPE,numberOfSourceComponents,err,error,*999)
        ENDIF
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        IF(equationsSource%sourceFieldAutoCreated) THEN
          !Finish creating the source field
          CALL Field_CreateFinish(equationsSource%sourceField,err,error,*999)
          !Set the default values for the source field to 1.0
          DO componentIdx=1,numberOfSourceComponents
            CALL Field_ComponentValuesInitialise(equationsSource%sourceField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
              & componentIdx,1.0_DP,err,error,*999)
          ENDDO !componentIdx
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
      !-----------------------------------------------------------------
      ! A n a l y t i c  T y p e
      !-----------------------------------------------------------------
      NULLIFY(equationsAnalytic)
      CALL EquationsSet_AnalyticGet(equationsSet,equationsAnalytic,err,error,*999)
      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
      CALL Field_NumberOfComponentsGet(geometricField,FIELD_U_VARIABLE_TYPE,numberOfDimensions,err,error,*999)
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        CALL EquationsSet_AssertMaterialsIsFinished(equationsSet,err,error,*999)
        NULLIFY(equationsMaterials)
        CALL EquationsSet_MaterialsGet(equationsSet,equationsMaterials,err,error,*999)
        NULLIFY(equationsSource)
        CALL EquationsSet_SourceExists(equationsSet,equationsSource,err,error,*999)
        NULLIFY(materialsField)
        CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
            IF(numberOfDimensions/=1) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a no source diffusion equation requires that there be 1 geometric dimension."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            IF(ASSOCIATED(equationsSource)) THEN
              localError="Cannot have a source for an analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Check the materials values are constant
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            !Set number of analytic field components
            numberOfAnalyticComponents=4
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1
          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
            !Check that domain is 2D
            IF(numberOfDimensions/=2) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a no source diffusion equation requires that there be 2 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            IF(ASSOCIATED(equationsSource)) THEN
              localError="Cannot have a source for an analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1
          CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
            !Check that domain is 3D
            IF(numberOfDimensions/=3) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a constant source diffusion equation requires that there be 3 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a generalised diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1)
            !Check that domain is 3D
            IF(numberOfDimensions/=3) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a linear source diffusion equation requires that there be 3 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a linear source diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
            !Check that domain is 1D
            IF(numberOfDimensions/=1) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " requires that there be 1 geometric dimension."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Check the materials values are constant
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,3,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            !Check that the a parameter is zero.
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,aParam,err,error,*999)
            IF(ABS(aParam)>ZERO_TOLERANCE) CALL FlagError("The 1st material component must be zero.",err,error,*999)
            !Check that the b parameter is not zero.
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,bParam,err,error,*999)
            IF(bParam<ZERO_TOLERANCE) CALL FlagError("The 2nd material component must be greater than zero.",err,error,*999)
            !Check to ensure we get real solutions
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3,cParam,err,error,*999)
            IF((bParam*cParam)>ZERO_TOLERANCE) &
              & CALL FlagError("The product of the 2nd and 3rd material components must not be positive.",err,error,*999)
            !Set the number of analytic field components
            numberOfAnalyticComponents=1
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a diffusion equation with a quadratic source."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1)
            !Check that domain is 1D
            IF(numberOfDimensions/=1) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " requires that there be 1 geometric dimension."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Check the materials values are constant
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,1,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            CALL Field_ComponentInterpolationCheck(materialsField,FIELD_U_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            CALL Field_ComponentInterpolationCheck(MaterialsField,FIELD_U_VARIABLE_TYPE,3,FIELD_CONSTANT_INTERPOLATION, &
              & err,error,*999)
            !Check that the a parameter is not zero.
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,aParam,err,error,*999)
            IF(ABS(aParam)<ZERO_TOLERANCE) CALL FlagError("The 1st material component must not be zero.",err,error,*999)
            !Check that the c parameter is not zero.
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,3,cParam,err,error,*999)
            IF(ABS(cParam)<ZERO_TOLERANCE) CALL FlagError("The 3rd material component must not be zero.",err,error,*999)
            !Check to ensure we get real solutions
            CALL Field_ParameterSetGetConstant(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,2,bParam,err,error,*999)
            IF((aParam*bParam)>ZERO_TOLERANCE) &
              & CALL FlagError("The product of the 1st and 2nd material components must not be positive.",err,error,*999)
            IF((aParam*cParam)<ZERO_TOLERANCE) &
              & CALL FlagError("The product of the 1st and 3rd material components must not be negative.",err,error,*999)
            !Set the number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_EQUATION_ONE_DIM_1
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a diffusion equation with an exponential source."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
          SELECT CASE(equationsSetSetup%analyticFunctionType)
          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM)
            !Check that domain is 2D
            IF(numberOfDimensions/=2) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a multi-compartment diffusion equation requires that there be 2 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_TWO_DIM
          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM)
            !Check that domain is 3D
            IF(numberOfDimensions/=3) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_MULTI_COMP_DIFFUSION_TWO_COMP_THREE_DIM
          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM)
            !Check that domain is 3D
            IF(numberOfDimensions/=3) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a multi-compartment diffusion equation requires that there be 3 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_MULTI_COMP_DIFFUSION_THREE_COMP_THREE_DIM
          CASE(EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM)
            !Check that domain is 3D
            IF(numberOfDimensions/=3) THEN
              localError="The number of geometric dimensions of "// &
                & TRIM(NumberToVString(numberOfDimensions,"*",err,error))// &
                & " is invalid. The analytic function type of "// &
                & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                & " for a multi-compartment diffusion requires that there be 3 geometric dimensions."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            !Set number of analytic field components
            numberOfAnalyticComponents=0
            !Set analytic function type
            equationsAnalytic%analyticFunctionType=EQUATIONS_SET_MULTI_COMP_DIFFUSION_FOUR_COMP_THREE_DIM
          CASE DEFAULT
            localError="The specified analytic function type of "// &
              & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
              & " is invalid for a multi-compartment diffusion equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The equation set subtype of "// &
            & TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
            & " is invalid for an analytical diffusion equation."
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
            NULLIFY(geometricDecomposition)
            CALL Field_DecompositionGet(geometricField,geometricDecomposition,err,error,*999)
            CALL Field_DecompositionSetAndLock(equationsAnalytic%analyticField,geometricDecomposition,err,error,*999)
            CALL Field_GeometricFieldSetAndLock(equationsAnalytic%analyticField,geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(equationsAnalytic%analyticField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE,err,error,*999)
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
            CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE)
              SELECT CASE(equationsAnalytic%analyticFunctionType)
              CASE(EQUATIONS_SET_DIFFUSION_EQUATION_ONE_DIM_1)
                !Set A
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                !Set B
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
                !Set C
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,3,1.0_DP,err,error,*999)
                !Set L
                CALL Field_ComponentValuesInitialise(equationsAnalytic%analyticField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)                      
              CASE(EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1)
                !Do nothing
              CASE DEFAULT
                localError="The specified analytic function type of "// &
                  & TRIM(NumberToVString(equationsSetSetup%analyticFunctionType,"*",err,error))// &
                  & " is invalid for a generalised diffusion equation."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE)
              !Do nothing
            CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
              !Do nothing
            CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
              !Do nothing
            CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
              !Do nothing
            CASE DEFAULT
              localError="The equation set subtype of "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                & " is invalid for an analytical diffusion equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ENDIF
        ENDIF
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
      !-----------------------------------------------------------------
      ! E q u a t i o n s    t y p e
      !-----------------------------------------------------------------     
      SELECT CASE(equationsSetSetup%actionType)
      CASE(EQUATIONS_SET_SETUP_START_ACTION)
        CALL EquationsSet_AssertDependentIsFinished(equationsSet,err,error,*999)
        NULLIFY(equations)
        CALL Equations_CreateStart(equationsSet,equations,err,error,*999)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
        CASE DEFAULT
          CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
        END SELECT
        CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
      CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
        CALL EquationsSet_SolutionMethodGet(equationsSet,solutionMethod,err,error,*999)
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
          CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
          SELECT CASE(esSpecification(3))
          CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
            & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
            CALL EquationsMappingVector_NumberOfResidualsSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_ResidualNumberOfVariablesSet(vectorMapping,1,1,err,error,*999)
            CALL EquationsMappingVector_ResidualVariableTypesSet(vectorMapping,1,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          CASE(EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_V_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
          CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUBTYPE,EQUATIONS_SET_MULTI_COMP_TRANSPORT_ADVEC_DIFF_SUPG_SUBTYPE)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
            NULLIFY(equationsSetField)
            CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
            CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
              & err,error,*999)
            myMatrixNumber = equationsSetFieldData(1)
            numberOfCompartments = equationsSetFieldData(2)    
            CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,numberOfCompartments-1,err,error,*999)           
            ALLOCATE(variableTypes(2*numberOfCompartments),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable types.",err,error,*999)
            ALLOCATE(variableUTypes(numberOfCompartments-1),STAT=err)
            IF(err/=0) CALL FlagError("Could not allocate variable U types.",err,error,*999)
            DO compartmentIdx=1,numberOfCompartments
              variableTypes(2*compartmentIdx-1)=FIELD_U_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(compartmentIdx-1))
              variableTypes(2*compartmentIdx)=FIELD_DELUDELN_VARIABLE_TYPE+(FIELD_NUMBER_OF_VARIABLE_SUBTYPES*(compartmentIdx-1))
            ENDDO !compartmentIdx
            compartmentCount=0
            DO compartmentIdx=1,numberOfCompartments
              IF(compartmentIdx/=myMatrixNumber)THEN
                compartmentCount=compartmentCount+1
                variableUTypes(compartmentCount)=variableTypes(2*compartmentIdx-1)
              ENDIF
            ENDDO !compartmentIdx
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,variableTypes(2*myMatrixNumber-1),err,error,*999)
            CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,variableUTypes,err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,variableTypes(2*myMatrixNumber),err,error,*999)
            IF(ALLOCATED(variableTypes)) DEALLOCATE(variableTypes)
            IF(ALLOCATED(variableUTypes)) DEALLOCATE(variableUTypes)
          CASE DEFAULT
            CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
            CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
          END SELECT
          NULLIFY(sourceField)
          CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
          IF(ASSOCIATED(sourceField)) THEN
            CALL EquationsMappingVector_NumberOfSourcesSet(vectorMapping,1,err,error,*999)
            CALL EquationsMappingVector_SourcesVariableTypesSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
          ENDIF
          CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
          !Create the equations matrices
          NULLIFY(vectorMatrices)
          CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
          !Set up matrix storage and structure
          CALL Equations_LumpingTypeGet(equations,lumpingType,err,error,*999)
          IF(lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
            !Set up lumping
            CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices,[EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED], &
              & err,error,*999)
            CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
              & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
            CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
              [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
            IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE.OR. &
              & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) THEN
              CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
              SELECT CASE(sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatricesVector_NonLinearStorageTypeSet(vectorMatrices,1,[DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatricesVector_NonLinearStorageTypeSet(vectorMatrices,1,couplingMatrixStorageType,err,error,*999)      
                CALL EquationsMatricesVector_NonLinearStructureTypeSet(vectorMatrices,1,couplingMatrixStructureType,err,error,*999)
              CASE DEFAULT
                localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          ELSE
            CALL Equations_SparsityTypeGet(equations,sparsityType,err,error,*999)
            SELECT CASE(sparsityType)
            CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
              CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices, &
                & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) THEN
                CALL EquationsMatricesVector_NonLinearStorageTypeSet(vectorMatrices,1, &
                  & DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
              ENDIF
            CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
              CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                & err,error,*999)
              CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE.OR. &
                & esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) THEN
                CALL EquationsMatricesVector_NonLinearStorageTypeSet(vectorMatrices,1,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)      
                CALL EquationsMatricesVector_NonLinearStructureTypeSet(vectorMatrices,1,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                  & err,error,*999)
              ELSE IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                ALLOCATE(couplingMatrixStorageType(numberOfCompartments-1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate coupling matrix storage type.",err,error,*999)
                ALLOCATE(couplingMatrixStructureType(numberOfCompartments-1),STAT=err)
                IF(err/=0) CALL FlagError("Could not allocate coupling matrix structure type.",err,error,*999)
                DO compartmentIdx=1,numberOfCompartments-1
                  couplingMatrixStorageType(compartmentIdx)=DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE
                  couplingMatrixStructureType(compartmentIdx)=EQUATIONS_MATRIX_FEM_STRUCTURE
                ENDDO !compartmentIdx
                CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,couplingMatrixStorageType,err,error,*999)      
                CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,couplingMatrixStructureType,err,error,*999)
                IF(ALLOCATED(couplingMatrixStorageType)) DEALLOCATE(couplingMatrixStorageType)
                IF(ALLOCATED(couplingMatrixStructureType)) DEALLOCATE(couplingMatrixStructureType)
              ENDIF
            CASE DEFAULT
              localError="The equations matrices sparsity type of "//TRIM(NumberToVString(sparsityType,"*",err,error))// &
                & " is invalid."
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
          localError="The solution method of "//TRIM(NumberToVString(solutionMethod,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(equationsSetSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(equationsSetSetup%setupType,"*",err,error))// &
        & " is invalid for a linear diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Diffusion_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Diffusion_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_EquationsSetSetup

  !
  !================================================================================================================================
  !
  
  !>Calculates the element stiffness matrices and RHS for a diffusion equation finite element equations set.
  SUBROUTINE Diffusion_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<Th e element number to calculate
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAXIMUM_COUPLING_MATRICES=99
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=3
    INTEGER(INTG) colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,columnXiIdx, &
      & compartmentCount,compartmentIdx,componentIdx,couplingVariableTypes(MAXIMUM_COUPLING_MATRICES), &
      & esSpecification(3),gaussPointIdx,myCompartment,meshCompartment1,meshCompartment2,numberOfColsComponents, &
      & numberOfColumnElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfCoupledColumnElementParameters,numberOfCompartments, &
      & numberOfCouplingComponents,numberOfDimensions,numberOfGauss,numberOfRowElementParameters(MAX_NUMBER_OF_COMPONENTS), &
      & numberOfRowsComponents,numberOfXi,rowComponentIdx,rowElementDOFIdx,rowElementParameterIdx,rowXiIdx,rowsVariableType, &
      & scalingType,xiIdx
    INTEGER(INTG), POINTER :: equationsSetFieldData(:)
    REAL(DP) :: aParam,bParam,cParam,columnComponentPhi(MAX_NUMBER_OF_COMPONENTS),columnPhi, &
      & columndPhidXi(MAX_NUMBER_OF_COMPONENTS),compartmentParam,conductivity(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_COMPONENTS), &
      & couplingParam,gaussWeight,jacobian,jacobianGaussWeight,kParam,rowComponentPhi(MAX_NUMBER_OF_COMPONENTS),rowPhi, &
      & rowdPhidXi(MAX_NUMBER_OF_COMPONENTS),sourceParam,sum
    LOGICAL :: update,updateCoupling,updateCouplingMatrices(MAXIMUM_COUPLING_MATRICES),updateDamping,updateMatrices,updateRHS, &
      & updateSource,updateStiffness
    TYPE(BasisType), POINTER :: coupledColumnBasis,dependentBasis,geometricBasis
    TYPE(BasisPtrType) :: columnBasis(MAX_NUMBER_OF_COMPONENTS),rowBasis(MAX_NUMBER_OF_COMPONENTS)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsMatrixPtrType) :: couplingMatrices(MAXIMUM_COUPLING_MATRICES) 
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,equationsSetField,fibreField,geometricField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: advecDiffDependentCurrentInterpParameters, &
      & advecDiffDependentPreviousInterpParameters,colsInterpParameters,couplingInterpParameters, &
      & diffusionDependentPreviousInterpParameters,fibreInterpParameters,geometricInterpParameters, &
      & rowsInterpParameters,sourceInterpParameters,uMaterialsInterpParameters,vMaterialsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: advecDiffDependentCurrentInterpPoint, &
      & advecDiffDependentPreviousInterpPoint,diffusionDependentPreviousInterpPoint,fibreInterpPoint,geometricInterpPoint, &
      & sourceInterpPoint,uMaterialsInterpPoint,vMaterialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,geometricVariable,rowsVariable
    TYPE(FieldVariablePtrType) :: couplingVariables(MAXIMUM_COUPLING_MATRICES)
    TYPE(QuadratureSchemeType), POINTER :: coupledColumnQuadratureScheme,dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAX_NUMBER_OF_COMPONENTS),rowQuadratureScheme(MAX_NUMBER_OF_COMPONENTS)
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Diffusion_FiniteElementCalculate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
    
    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &      
     & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
     & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
     & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
     & EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Can not calculate finite element stiffness matrices for a nonlinear source.",err,error,*999)
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field equations set class."
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
    NULLIFY(dynamicMapping)
    CALL EquationsMappingVector_DynamicMappingGet(vectorMapping,dynamicMapping,err,error,*999)
    NULLIFY(vectorMatrices)
    CALL EquationsVector_VectorMatricesGet(vectorEquations,vectorMatrices,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    NULLIFY(dampingMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    NULLIFY(sourcesMapping)
    CALL EquationsMappingVector_SourcesMappingExists(vectorMapping,sourcesMapping,err,error,*999)
    NULLIFY(sourceMapping)
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMappingSources_SourceMappingGet(sourcesMapping,1,sourceMapping,err,error,*999)
    ENDIF
    NULLIFY(sourceField)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsSet_SourceFieldGet(equationsSet,sourceField,err,error,*999)
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    NULLIFY(rhsMapping)
    CALL EquationsMappingVector_RHSMappingExists(vectorMapping,rhsMapping,err,error,*999)
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF
    updateCoupling=.FALSE.
    IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
      NULLIFY(equationsSetField)
      NULLIFY(equationsSetFieldData)
      CALL EquationsSet_EquationsSetFieldGet(equationsSet,equationsSetField,err,error,*999)
      CALL Field_ParameterSetDataGet(equationsSetField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,equationsSetFieldData, &
        & err,error,*999)
      myCompartment=equationsSetFieldData(1)
      numberOfCompartments=equationsSetFieldData(2)
      IF(numberOfCompartments>MAXIMUM_COUPLING_MATRICES) THEN
        localError="The number of coupling compartments of "//TRIM(NumberToVString(numberOfCompartments,"*",err,error))// &
          & " is greater than the maximum number of coupling compartments of "// &
          & TRIM(NumberToVString(MAXIMUM_COUPLING_MATRICES,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      NULLIFY(linearMapping)
      CALL EquationsMappingVector_LinearMappingGet(vectorMapping,linearMapping,err,error,*999)
      NULLIFY(linearMatrices)
      CALL EquationsMatricesVector_LinearMatricesGet(vectorMatrices,linearMatrices,err,error,*999)      
      compartmentCount=0
      DO compartmentIdx=1,numberOfCompartments
        IF(compartmentIdx/=myCompartment)THEN
          compartmentCount=compartmentCount+1
          NULLIFY(couplingMatrices(compartmentCount)%ptr)
          CALL EquationsMatricesLinear_EquationsMatrixGet(linearMatrices,compartmentCount,couplingMatrices(compartmentCount)%ptr, &
            & err,error,*999)
          CALL EquationsMatrix_UpdateMatrixGet(couplingMatrices(compartmentCount)%ptr,updateCouplingMatrices(compartmentCount), &
            & err,error,*999)
          updateCoupling=(updateCoupling.OR.updateCouplingMatrices(compartmentCount))
          NULLIFY(couplingVariables(compartmentCount)%ptr)
          CALL EquationsMappingLinear_LinearMatrixVariableGet(linearMapping,compartmentCount, &
            & couplingVariables(compartmentCount)%ptr,err,error,*999)
          CALL EquationsMappingLinear_LinearMatrixVariableTypeGet(linearMapping,compartmentCount, &
            & couplingVariableTypes(compartmentCount),err,error,*999)
        ENDIF
      ENDDO !compartmentIdx
    ENDIF

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateCoupling.OR.updateSource.OR.updateRHS)

    IF(update) THEN
 
      NULLIFY(equationsInterpolation)
      CALL Equations_InterpolationGet(equations,equationsInterpolation,err,error,*999)

      NULLIFY(geometricField)
      CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
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
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
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

      NULLIFY(fibreField)
      NULLIFY(fibreInterpParameters)
      NULLIFY(fibreInterpPoint)
      CALL EquationsSet_FibreFieldExists(equationsSet,fibreField,err,error,*999)
      IF(ASSOCIATED(fibreField)) THEN
        CALL EquationsInterpolation_FibreParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpParameters, &
          & err,error,*999)
        CALL EquationsInterpolation_FibrePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,fibreInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,fibreInterpParameters,err,error,*999)
      ENDIF      
 
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
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      
      NULLIFY(materialsField)
      CALL EquationsSet_MaterialsFieldGet(equationsSet,materialsField,err,error,*999)
      NULLIFY(uMaterialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & uMaterialsInterpParameters,err,error,*999)
      NULLIFY(uMaterialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,uMaterialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,uMaterialsInterpParameters,err,error,*999)
      NULLIFY(vMaterialsInterpParameters)
      NULLIFY(vMaterialsInterpPoint)
      IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
        CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & vMaterialsInterpParameters,err,error,*999)
        CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE,vMaterialsInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,vMaterialsInterpParameters,err,error,*999)
      ENDIF
      
      IF(esSpecificatIon(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
        NULLIFY(advecDiffDependentCurrentInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & advecDiffDependentCurrentInterpParameters,err,error,*999)
        NULLIFY(advecDiffDependentCurrentInterpPoint)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & advecDiffDependentCurrentInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber, &
          & advecDiffDependentCurrentInterpParameters,err,error,*999)
        NULLIFY(advecDiffDependentPreviousInterpParameters)
        CALL EquationsInterpolation_PreviousDependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & advecDiffDependentPreviousInterpParameters,err,error,*999)
        NULLIFY(advecDiffDependentPreviousInterpPoint)
        CALL EquationsInterpolation_PreviousDependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
          & advecDiffDependentPreviousInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_PREVIOUS_VALUES_SET_TYPE,elementNumber, &
          & advecDiffDependentPreviousInterpParameters,err,error,*999)          
        NULLIFY(diffusionDependentPreviousInterpParameters)
        CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & diffusionDependentPreviousInterpParameters,err,error,*999)
        NULLIFY(diffusionDependentPreviousInterpPoint)
        CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_V_VARIABLE_TYPE, &
          & diffusionDependentPreviousInterpPoint,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_PREVIOUS_VALUES_SET_TYPE, elementNumber, &
          & diffusionDependentPreviousInterpParameters,err,error,*999)
      ENDIF
      
      NULLIFY(sourceInterpParameters)
      NULLIFY(sourceInterpPoint)
      NULLIFY(sourceField)
      IF(updateSource) THEN
        CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)      
        IF(ASSOCIATED(sourceField)) THEN
          CALL EquationsInterpolation_SourceParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpParameters, &
            & err,error,*999)
          CALL EquationsInterpolation_SourcePointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,sourceInterpPoint,err,error,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,sourceInterpParameters,err,error,*999)
        ENDIF
      ENDIF
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingDynamic_DynamicVariableGet(dynamicMapping,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999) 
      
      DO rowComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBasis(rowComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis(rowComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBasis(rowComponentIdx)%ptr,numberOfRowElementParameters(rowComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureScheme(rowComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBasis(rowComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureScheme(rowComponentIdx)%ptr,err,error,*999)
      ENDDO !rowComponentIdx
      IF(updateMatrices) THEN
        DO columnComponentIdx=1,numberOfColsComponents
          NULLIFY(columnDomain)
          CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          NULLIFY(columnBasis(columnComponentIdx)%ptr)
          CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis(columnComponentIdx)%ptr,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(columnBasis(columnComponentIdx)%ptr, &
            & numberOfColumnElementParameters(columnComponentIdx),err,error,*999)
          NULLIFY(columnQuadratureScheme(columnComponentIdx)%ptr)
          CALL Basis_QuadratureSchemeGet(columnBasis(columnComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
            & columnQuadratureScheme(columnComponentIdx)%ptr,err,error,*999)
        ENDDO !columnComponentIdx
      ENDIF
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        IF(ASSOCIATED(fibreField)) CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
          & fibreInterpPoint,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,uMaterialsInterpPoint, &
          & err,error,*999)
        IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,vMaterialsInterpPoint, &
            & err,error,*999)
          couplingParam=vMaterialsInterpPoint%values(myCompartment,NO_PART_DERIV)
        ENDIF
        IF(ASSOCIATED(sourceField)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
            & err,error,*999)
          sourceParam=sourceInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        
        aParam=uMaterialsInterpPoint%values(1,NO_PART_DERIV)
        SELECT CASE(esSpecification(3))
        CASE(EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
          & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE)
          bParam=uMaterialsInterpPoint%values(2,NO_PART_DERIV)
          !Calculate conductivity tensor
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & uMaterialsInterpPoint%values(3:3+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)
        CASE DEFAULT
          !Calculate conductivity tensor
          CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
            & uMaterialsInterpPoint%values(2:2+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)  
        END SELECT

        IF(esSpecification(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & advecDiffDependentCurrentInterpPoint,err,error,*999)
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & advecDiffDependentPreviousInterpPoint,err,error,*999)           
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx, &
            & diffusionDependentPreviousInterpPoint,err,error,*999)
          !cParam_1_T=dependentInterpPoint%value(1,NO_PART_DERIV) for the U variable type
          !   This is the value of the solution from the advection-diffusion equation at time T          
          !cParam_1_TPlusOne=
          !   This is the value of the solution from the advection-diffusion equation at time T+deltaT
          !cParam_2_T=dependentInterpPoint%value(1,NO_PART_DERIV) for the V variable type
          !   This is the value of the solution from the diffusion equation at time T
          !cParam=cParam_1_T+cParam_1_TPlusOne+cParam_2_T
          !The value of the rhs term is +0.5*(C_1^{t}+C_1^{t+1}-C_2^{t}) 
          cParam=0.5_DP*advecDiffDependentPreviousInterpPoint%VALUES(1,NO_PART_DERIV)+ &
            & 0.5_DP*advecDiffDependentCurrentInterpPoint%values(1,NO_PART_DERIV)- &
            & diffusionDependentPreviousInterpPoint%values(1,NO_PART_DERIV)          
        ENDIF
        
        !Calculate jacobianGaussWeight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
        !Loop over field components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
              & NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                 DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                    & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  DO xiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                      & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx, &
                      & columndPhidXi(xiIdx),err,error,*999)
                  ENDDO !xiIdx
                  IF(updateStiffness) THEN
                    sum=0.0_DP
                    DO rowXiIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    IF(esSpecification(3)==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE.OR. &
                      & esSpecification(3)==EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE) THEN
                      sum=sum+bParam*rowPhi*columnPhi
                    ELSE IF(esSpecification(3)==EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE) THEN
                      sum=sum*couplingParam
                    ENDIF
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                  IF(updateDamping) THEN
                    sum=aParam*rowPhi*columnPhi
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !updateMatrices
            IF(updateCoupling) THEN
              compartmentCount=0
              DO compartmentIdx=1,numberOfCompartments
                IF(compartmentIdx/=myCompartment) THEN
                  compartmentCount=compartmentCount+1
                  !need to test for the case where compartmentIdx==mycompartment
                  !the coupling terms then needs to be added into the stiffness matrix
                  IF(updateCouplingMatrices(compartmentCount)) THEN
                    compartmentParam=vMaterialsInterpPoint%values(compartmentIdx,NO_PART_DERIV)

!! TODO: Cache the quadrature schemes etc.
                    
                    !Loop over element columns
                    columnElementDOFIdx=0
                    NULLIFY(couplingVariables(compartmentCount)%ptr)
                    CALL FieldVariable_NumberOfComponentsGet(couplingVariables(compartmentCount)%ptr,numberOfCouplingComponents, &
                      & err,error,*999)
                    DO columnComponentIdx=1,numberOfCouplingComponents
                      NULLIFY(columnDomain)
                      CALL FieldVariable_ComponentDomainGet(couplingVariables(compartmentCount)%ptr,columnComponentIdx, &
                        & columnDomain,err,error,*999)
                      NULLIFY(columnDomainTopology)
                      CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                      NULLIFY(columnDomainElements)
                      CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                      NULLIFY(coupledColumnBasis)
                      CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,coupledColumnBasis,err,error,*999)
                      NULLIFY(coupledColumnQuadratureScheme)
                      CALL Basis_QuadratureSchemeGet(coupledColumnBasis,BASIS_DEFAULT_QUADRATURE_SCHEME, &
                        & coupledColumnQuadratureScheme,err,error,*999)
                      CALL Basis_NumberOfElementParametersGet(coupledColumnBasis,numberOfCoupledColumnElementParameters, &
                        & err,error,*999)
                      DO columnElementParameterIdx=1,numberOfCoupledColumnElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        CALL BasisQuadratureScheme_GaussBasisFunctionGet(coupledColumnQuadratureScheme,columnElementParameterIdx, &
                          & NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)                        
                        sum=compartmentParam*rowPhi*columnPhi                        
                        couplingMatrices(compartmentCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & couplingMatrices(compartmentCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + & 
                          & sum*jacobianGaussWeight
                      ENDDO !columnElementParameterIdx
                    ENDDO !columnComponentIdx
                  ENDIF
                ENDIF
              ENDDO !compartmentIdx
            ENDIF
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                & sourceParam*rowPhi*jacobianGaussWeight               
            ENDIF
            IF(updateRHS) THEN
              IF(esSpecification(3)==EQUATIONS_SET_COUPLED_SOURCE_DIFFUSION_ADVEC_DIFF_SUBTYPE) THEN
                rhsVector%elementVector%vector(rowElementDOFIdx)=rhsVector%elementVector%vector(rowElementDOFIdx)- &
                  & cParam*rowPhi*jacobianGaussWeight
              ELSE
                rhsVector%elementVector%vector(rowElementDOFIdx)=0.0_DP
              ENDIF
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
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
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
            IF(updateCoupling) THEN
              compartmentCount=0
              DO compartmentIdx=1,numberOfCompartments
                IF(compartmentIdx/=myCompartment) THEN
                  compartmentCount=compartmentCount+1
                  IF(updateCouplingMatrices(compartmentCount)) THEN

!!TODO : Cache number of parameters etc.
                    
                    !Loop over element columns
                    columnElementDOFIdx=0
                    CALL FieldVariable_NumberOfComponentsGet(couplingVariables(compartmentCount)%ptr,numberOfCouplingComponents, &
                      & err,error,*999)
                    NULLIFY(couplingInterpParameters)
                    CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation, &
                      & couplingVariableTypes(compartmentCount),couplingInterpParameters,err,error,*999)
                    CALL Field_InterpolationParametersScaleFactorsElementGet(elementNumber,couplingInterpParameters,err,error,*999)
                    DO columnComponentIdx=1,numberOfCouplingComponents
                      NULLIFY(columnDomain)
                      CALL FieldVariable_ComponentDomainGet(couplingVariables(compartmentCount)%ptr,columnComponentIdx, &
                        & columnDomain,err,error,*999)
                      NULLIFY(columnDomainTopology)
                      CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
                      NULLIFY(columnDomainElements)
                      CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
                      NULLIFY(coupledColumnBasis)
                      CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,coupledColumnBasis,err,error,*999)
                      CALL Basis_NumberOfElementParametersGet(coupledColumnBasis,numberOfCoupledColumnElementParameters, &
                        & err,error,*999)
                      DO columnElementParameterIdx=1,numberOfCoupledColumnElementParameters
                        columnElementDOFIdx=columnElementDOFIdx+1
                        couplingMatrices(compartmentCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) = &
                          & couplingMatrices(compartmentCount)%ptr%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx) + & 
                          & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                          & couplingInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
                      ENDDO !columnElementParameterGet
                    ENDDO !columnComponentIdx
                  ENDIF
                ENDIF
              ENDDO !compartmentIdx
            ENDIF !updateCoupling
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
         
    EXITS("Diffusion_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Diffusion_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a diffusion equation finite element equations set.
  SUBROUTINE Diffusion_FiniteElementJacobianEvaluate(equationsSet,elementNumber,err,error,*)
    
    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element Jacobian evaluation on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=3
    INTEGER(INTG) colsVariableType,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx,esSpecification(3), &
      & gaussPointIdx,numberOfColsComponents,numberOfColumnElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfDimensions, &
      & numberOfGauss,numberOfRowElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfRowsComponents,numberOfXi,rowComponentIdx, &
      & rowElementDOFIdx,rowElementParameterIdx,rowsVariableType,scalingType
    REAL(DP) :: bParam,cParam,columnPhi,gaussWeight,jacobian,jacobianGaussWeight,rowPhi,sum,uValue,VALUE
    LOGICAL :: updateJacobian
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(BasisPtrType) :: columnBasis(MAX_NUMBER_OF_COMPONENTS),rowBasis(MAX_NUMBER_OF_COMPONENTS)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,geometricField,materialsField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,dependentInterpParameters,geometricInterpParameters, &
      & materialsInterpParameters,rowsInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,geometricInterpPoint,materialsInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,geometricVariable,rowsVariable
    TYPE(JacobianMatrixType), POINTER :: jacobianMatrix
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAX_NUMBER_OF_COMPONENTS),rowQuadratureScheme(MAX_NUMBER_OF_COMPONENTS)
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("Diffusion_FiniteElementJacobianEvaluate",err,error,*999)

    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
      CALL FlagError("Can not evaluate a Jacobian for a linear diffusion equation.",err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)      
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field equations set class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
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
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
      
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
      
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      !Cache row and column bases and quadrature schemes to avoid repeated calculations
      IF(numberOfRowsComponents>MAX_NUMBER_OF_COMPONENTS) THEN
        localError="The number of rows components of "//TRIM(NumberToVString(numberOfRowsComponents,"*",err,error))// &
          & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAX_NUMBER_OF_COMPONENTS,"*", &
          & err,error))//". Increase MAX_NUMBER_OF_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO rowComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBasis(rowComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis(rowComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBasis(rowComponentIdx)%ptr,numberOfRowElementParameters(rowComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureScheme(rowComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBasis(rowComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureScheme(rowComponentIdx)%ptr,err,error,*999)
      ENDDO !rowComponentIdx
      IF(numberOfColsComponents>MAX_NUMBER_OF_COMPONENTS) THEN
        localError="The number of columns components of "//TRIM(NumberToVString(numberOfColsComponents,"*",err,error))// &
          & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAX_NUMBER_OF_COMPONENTS,"*", &
          & err,error))//". Increase MAX_NUMBER_OF_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO columnComponentIdx=1,numberOfColsComponents
        NULLIFY(columnDomain)
        CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
        NULLIFY(columnDomainTopology)
        CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
        NULLIFY(columnDomainElements)
        CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
        NULLIFY(columnBasis(columnComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis(columnComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(columnBasis(columnComponentIdx)%ptr, &
          & numberOfColumnElementParameters(columnComponentIdx),err,error,*999)
        NULLIFY(columnQuadratureScheme(columnComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(columnBasis(columnComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & columnQuadratureScheme(columnComponentIdx)%ptr,err,error,*999)
      ENDDO !columnComponentIdx
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics,err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        
        uValue=dependentInterpPoint%values(1,NO_PART_DERIV)
        bParam=materialsInterpPoint%values(2,NO_PART_DERIV)
        IF(esSpecification(3)==EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE) &
          & cParam=materialsInterpPoint%values(3,NO_PART_DERIV)

        !Calculate jacobianGaussWeight.
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
         
        !Loop over field components
        rowElementDOFIdx=0
        DO rowComponentIdx=1,numberOfRowsComponents
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
              & NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
            !Loop over element columns
            columnElementDOFIdx=0
            DO columnComponentIdx=1,numberOfColsComponents
              DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                columnElementDOFIdx=columnElementDOFIdx+1
                CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                  & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                SELECT CASE(esSpecification(3))
                CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE)
                  sum=-2.0_DP*bParam*rowPhi*columnPhi*uValue
                CASE(EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
                  sum=-bParam*cParam*rowPhi*columnPhi*EXP(cParam*uValue)                  
                CASE DEFAULT
                  localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
                    & " is not valid for a diffusion equation type of a classical field equations set class."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                
                jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                  & sum*jacobianGaussWeight
                
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
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            !Loop over element columns
            columnElementDOFIdx=0
            DO columnComponentIdx=1,numberOfColsComponents
              DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                columnElementDOFIdx=columnElementDOFIdx+1
                jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                  & jacobianMatrix%elementJacobian%matrix(rowElementDOFIdx,columnElementDOFIdx)* &
                  & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)* &
                  & colsInterpParameters%scaleFactors(columnElementParameterIdx,columnComponentIdx)
              ENDDO !columnElementParameterIdx
            ENDDO !columnComponentIdx
          ENDDO !rowElementParameterIdx
        ENDDO !rowComponentIdx
      ENDIF !scaling
    ENDIF !update
           
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
  SUBROUTINE Diffusion_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: elementNumber !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG), PARAMETER :: MAX_NUMBER_OF_COMPONENTS=3
    INTEGER(INTG) :: colsVariableType,componentIdx,columnComponentIdx,columnElementDOFIdx,columnElementParameterIdx, &
      & columnXiIdx,esSpecification(3),gaussPointIdx,numberOfColsComponents, &
      & numberOfColumnElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfDimensions,numberOfGauss,numberOfRowsComponents, &
      & numberOfRowElementParameters(MAX_NUMBER_OF_COMPONENTS),numberOfXi,rowComponentIdx,rowElementDOFIdx, &
      & rowElementParameterIdx,rowsVariableType,rowXiIdx,scalingType,xiIdx
    REAL(DP) :: aParam,bParam,cParam,columnPhi,columndPhidXi(MAX_NUMBER_OF_COMPONENTS), &
      & conductivity(MAX_NUMBER_OF_COMPONENTS,MAX_NUMBER_OF_COMPONENTS),gaussWeight,jacobian,jacobianGaussWeight, &
      & kParam,rowPhi,rowdPhidXi(MAX_NUMBER_OF_COMPONENTS),sourceValue,sum,rowComponentPhi(MAX_NUMBER_OF_COMPONENTS), &
      & columnComponentPhi(3),uValue
    LOGICAL :: update,updateDamping,updateMatrices,updateResidual,updateRHS,updateSource,updateStiffness
    TYPE(BasisType), POINTER :: dependentBasis,geometricBasis
    TYPE(BasisPtrType) :: columnBasis(MAX_NUMBER_OF_COMPONENTS),rowBasis(MAX_NUMBER_OF_COMPONENTS)
    TYPE(DecompositionType), POINTER :: dependentDecomposition,geometricDecomposition
    TYPE(DomainType), POINTER :: columnDomain,dependentDomain,geometricDomain,rowDomain
    TYPE(DomainElementsType), POINTER :: columnDomainElements,dependentDomainElements,geometricDomainElements,rowDomainElements
    TYPE(DomainTopologyType), POINTER :: columnDomainTopology,dependentDomainTopology,geometricDomainTopology,rowDomainTopology   
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsInterpolationType), POINTER :: equationsInterpolation
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingLHSType), POINTER :: lhsMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingResidualType), POINTER :: residualMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMappingSourcesType), POINTER :: sourcesMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesResidualType), POINTER :: residualVector
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatricesSourcesType), POINTER :: sourceVectors
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatrixType), POINTER :: stiffnessMatrix,dampingMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: dependentField,fibreField,geometricField,materialsField,sourceField
    TYPE(FieldInterpolationParametersType), POINTER :: colsInterpParameters,dependentInterpParameters,fibreInterpParameters, &
      & geometricInterpParameters,materialsInterpParameters,rowsInterpParameters,sourceInterpParameters
    TYPE(FieldInterpolatedPointType), POINTER :: dependentInterpPoint,fibreInterpPoint,geometricInterpPoint,materialsInterpPoint, &
      & sourceInterpPoint
    TYPE(FieldInterpolatedPointMetricsType), POINTER :: geometricInterpPointMetrics
    TYPE(FieldVariableType), POINTER :: colsVariable,dependentVariable,geometricVariable,rowsVariable
    TYPE(QuadratureSchemeType), POINTER :: dependentQuadratureScheme,geometricQuadratureScheme
    TYPE(QuadratureSchemePtrType) :: columnQuadratureScheme(MAX_NUMBER_OF_COMPONENTS),rowQuadratureScheme(MAX_NUMBER_OF_COMPONENTS)
    TYPE(VARYING_STRING) :: localError
     
    ENTERS("Diffusion_FiniteElementResidualEvaluate",err,error,*999)
    
    CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)

    SELECT CASE(esSpecification(3))
    CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_MULTI_COMP_TRANSPORT_DIFFUSION_SUBTYPE)
      CALL FlagError("Can not evaluate a residual for a linear diffusion equation.",err,error,*999)
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
      !OK
    CASE(EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
      & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)      
    CASE DEFAULT
      localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field equations set class."
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
    CALL EquationsMatricesResidual_UpdateVectorGet(residualVector,updateResidual,err,error,*999)
    NULLIFY(dynamicMatrices)
    CALL EquationsMatricesVector_DynamicMatricesGet(vectorMatrices,dynamicMatrices,err,error,*999)
    NULLIFY(stiffnessMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,1,stiffnessMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(stiffnessMatrix,updateStiffness,err,error,*999)
    NULLIFY(dampingMatrix)
    CALL EquationsMatricesDynamic_EquationsMatrixGet(dynamicMatrices,2,dampingMatrix,err,error,*999)
    CALL EquationsMatrix_UpdateMatrixGet(dampingMatrix,updateDamping,err,error,*999)
    NULLIFY(sourceVectors)
    NULLIFY(sourceVector)
    updateSource=.FALSE.
    IF(ASSOCIATED(sourcesMapping)) THEN
      CALL EquationsMatricesVector_SourceVectorsGet(vectorMatrices,sourceVectors,err,error,*999)
      CALL EquationsMatricesSources_SourceVectorGet(sourceVectors,1,sourceVector,err,error,*999)
      CALL EquationsMatricesSource_UpdateVectorGet(sourceVector,updateSource,err,error,*999)
    ENDIF
    NULLIFY(rhsVector)
    updateRHS=.FALSE.
    IF(ASSOCIATED(rhsMapping)) THEN
      CALL EquationsMatricesVector_RHSVectorGet(vectorMatrices,rhsVector,err,error,*999)
      CALL EquationsMatricesRHS_UpdateVectorGet(rhsVector,updateRHS,err,error,*999)
    ENDIF

    updateMatrices=(updateStiffness.OR.updateDamping)
    update=(updateMatrices.OR.updateSource.OR.updateResidual.OR.updateRHS)

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
      NULLIFY(dependentDomain)
      CALL Decomposition_DomainGet(dependentDecomposition,0,dependentDomain,err,error,*999)
      NULLIFY(dependentDomainTopology)
      CALL Domain_DomainTopologyGet(dependentDomain,dependentDomainTopology,err,error,*999)
      NULLIFY(dependentDomainElements)
      CALL DomainTopology_DomainElementsGet(dependentDomainTopology,dependentDomainElements,err,error,*999)
      NULLIFY(dependentBasis)
      CALL DomainElements_ElementBasisGet(dependentDomainElements,elementNumber,dependentBasis,err,error,*999)
    
      NULLIFY(rowsVariable)
      CALL EquationsMappingLHS_LHSVariableGet(lhsMapping,rowsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(rowsVariable,rowsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(rowsVariable,numberOfRowsComponents,err,error,*999)
      
      NULLIFY(colsVariable)
      CALL EquationsMappingResidual_VariableGet(residualMapping,1,colsVariable,err,error,*999)
      CALL FieldVariable_VariableTypeGet(colsVariable,colsVariableType,err,error,*999)
      CALL FieldVariable_NumberOfComponentsGet(colsVariable,numberOfColsComponents,err,error,*999)
      
      NULLIFY(geometricQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(geometricBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,geometricQuadratureScheme,err,error,*999)
      NULLIFY(dependentQuadratureScheme)
      CALL Basis_QuadratureSchemeGet(dependentBasis,BASIS_DEFAULT_QUADRATURE_SCHEME,dependentQuadratureScheme,err,error,*999)
      CALL BasisQuadratureScheme_NumberOfGaussGet(dependentQuadratureScheme,numberOfGauss,err,error,*999)
      
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
      
      NULLIFY(dependentInterpParameters)
      CALL EquationsInterpolation_DependentParametersGet(equationsInterpolation,colsVariableType, &
        & dependentInterpParameters,err,error,*999)
      NULLIFY(dependentInterpPoint)
      CALL EquationsInterpolation_DependentPointGet(equationsInterpolation,colsVariableType,dependentInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,dependentInterpParameters,err,error,*999)
    
      NULLIFY(materialsInterpParameters)
      CALL EquationsInterpolation_MaterialsParametersGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE, &
        & materialsInterpParameters,err,error,*999)
      NULLIFY(materialsInterpPoint)
      CALL EquationsInterpolation_MaterialsPointGet(equationsInterpolation,FIELD_U_VARIABLE_TYPE,materialsInterpPoint, &
        & err,error,*999)
      CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,elementNumber,materialsInterpParameters,err,error,*999)
      
      NULLIFY(sourceInterpParameters)
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
      
      !Cache row and column bases and quadrature schemes to avoid repeated calculations
      IF(numberOfRowsComponents>MAX_NUMBER_OF_COMPONENTS) THEN
        localError="The number of rows components of "//TRIM(NumberToVString(numberOfRowsComponents,"*",err,error))// &
          & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAX_NUMBER_OF_COMPONENTS,"*", &
          & err,error))//". Increase MAX_NUMBER_OF_COMPONENTS."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      DO rowComponentIdx=1,numberOfRowsComponents
        NULLIFY(rowDomain)
        CALL FieldVariable_ComponentDomainGet(rowsVariable,rowComponentIdx,rowDomain,err,error,*999)
        NULLIFY(rowDomainTopology)
        CALL Domain_DomainTopologyGet(rowDomain,rowDomainTopology,err,error,*999)
        NULLIFY(rowDomainElements)
        CALL DomainTopology_DomainElementsGet(rowDomainTopology,rowDomainElements,err,error,*999)
        NULLIFY(rowBasis(rowComponentIdx)%ptr)
        CALL DomainElements_ElementBasisGet(rowDomainElements,elementNumber,rowBasis(rowComponentIdx)%ptr,err,error,*999)
        CALL Basis_NumberOfElementParametersGet(rowBasis(rowComponentIdx)%ptr,numberOfRowElementParameters(rowComponentIdx), &
          & err,error,*999)
        NULLIFY(rowQuadratureScheme(rowComponentIdx)%ptr)
        CALL Basis_QuadratureSchemeGet(rowBasis(rowComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
          & rowQuadratureScheme(rowComponentIdx)%ptr,err,error,*999)
      ENDDO !rowComponentIdx
      IF(updateMatrices) THEN
        IF(numberOfColsComponents>MAX_NUMBER_OF_COMPONENTS) THEN
          localError="The number of columns components of "//TRIM(NumberToVString(numberOfColsComponents,"*",err,error))// &
            & " is greater than the maximum number of components of "//TRIM(NumberToVString(MAX_NUMBER_OF_COMPONENTS,"*", &
            & err,error))//". Increase MAX_NUMBER_OF_COMPONENTS."
          CALL FlagError(localError,err,error,*999)
        ENDIF
        DO columnComponentIdx=1,numberOfColsComponents
          NULLIFY(columnDomain)
          CALL FieldVariable_ComponentDomainGet(colsVariable,columnComponentIdx,columnDomain,err,error,*999)
          NULLIFY(columnDomainTopology)
          CALL Domain_DomainTopologyGet(columnDomain,columnDomainTopology,err,error,*999)
          NULLIFY(columnDomainElements)
          CALL DomainTopology_DomainElementsGet(columnDomainTopology,columnDomainElements,err,error,*999)
          NULLIFY(columnBasis(columnComponentIdx)%ptr)
          CALL DomainElements_ElementBasisGet(columnDomainElements,elementNumber,columnBasis(columnComponentIdx)%ptr,err,error,*999)
          CALL Basis_NumberOfElementParametersGet(columnBasis(columnComponentIdx)%ptr, &
            & numberOfColumnElementParameters(columnComponentIdx),err,error,*999)
          NULLIFY(columnQuadratureScheme(columnComponentIdx)%ptr)
          CALL Basis_QuadratureSchemeGet(columnBasis(columnComponentIdx)%ptr,BASIS_DEFAULT_QUADRATURE_SCHEME, &
            & columnQuadratureScheme(columnComponentIdx)%ptr,err,error,*999)
        ENDDO !columnComponentIdx
      ENDIF
      
      !Loop over gauss points
      DO gaussPointIdx=1,numberOfGauss
        
        CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,geometricInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolatedPointMetricsCalculate(numberOfXi,geometricInterpPointMetrics, &
          & err,error,*999)
        IF(ASSOCIATED(fibreField)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,fibreInterpPoint, &
            & err,error,*999)
        ENDIF
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,dependentInterpPoint, &
          & err,error,*999)
        CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,materialsInterpPoint, &
          & err,error,*999)
        IF(ASSOCIATED(sourceField)) THEN
          CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gaussPointIdx,sourceInterpPoint, &
            & err,error,*999)
          sourceValue=sourceInterpPoint%values(1,NO_PART_DERIV)
        ENDIF
        
        uValue=dependentInterpPoint%values(1,NO_PART_DERIV)
        aParam=materialsInterpPoint%values(1,NO_PART_DERIV)
        bParam=materialsInterpPoint%values(2,NO_PART_DERIV)
        cParam=materialsInterpPoint%values(3,NO_PART_DERIV)
        
        CALL CoordinateSystem_MaterialTransformSymTensor2(geometricInterpPointMetrics,fibreInterpPoint, &
          & materialsInterpPoint%values(4:4+NUMBER_OF_VOIGT(numberOfDimensions),NO_PART_DERIV),conductivity,err,error,*999)
        
        !Calculate jacobianGaussWeight.
!!TODO: Think about symmetric problems. 
        CALL FieldInterpolatedPointMetrics_JacobianGet(geometricInterpPointMetrics,jacobian,err,error,*999)
        CALL BasisQuadratureScheme_GaussWeightGet(dependentQuadratureScheme,gaussPointIdx,gaussWeight,err,error,*999)
        jacobianGaussWeight=jacobian*gaussWeight
        
        !Loop over field components
        rowElementDOFIdx=0          
        DO rowComponentIdx=1,numberOfRowsComponents
          !Loop over element rows
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1
            CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
              & NO_PART_DERIV,gaussPointIdx,rowPhi,err,error,*999)
            DO xiIdx=1,numberOfXi
              CALL BasisQuadratureScheme_GaussBasisFunctionGet(rowQuadratureScheme(rowComponentIdx)%ptr,rowElementParameterIdx, &
                & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx,rowdPhidXi(xiIdx),err,error,*999)
            ENDDO !xiIdx
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
                  columnElementDOFIdx=columnElementDOFIdx+1
                  CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                    & columnElementParameterIdx,NO_PART_DERIV,gaussPointIdx,columnPhi,err,error,*999)
                  DO xiIdx=1,numberOfXi
                    CALL BasisQuadratureScheme_GaussBasisFunctionGet(columnQuadratureScheme(columnComponentIdx)%ptr, &
                      & columnElementParameterIdx,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(xiIdx),gaussPointIdx, &
                      & columndPhidXi(xiIdx),err,error,*999)
                  ENDDO !xiIdx
                  IF(updateStiffness) THEN
                    sum=0.0_DP
                    DO rowXiIdx=1,numberOfXi
                      DO columnXiIdx=1,numberOfXi
                        DO xiIdx=1,numberOfXi
                          sum=sum+conductivity(rowXiIdx,xiIdx)*rowdPhidXi(rowXiIdx)*columndPhidXi(columnXiIdx)* &
                            & geometricInterpPointMetrics%gu(rowXiIdx,xiIdx)
                        ENDDO !xiIdx
                      ENDDO !columnXiIdx
                    ENDDO !rowXiIdx
                    IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) sum=sum+bParam*rowPhi*columnPhi
                    stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & stiffnessMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                  IF(updateDamping) THEN
                    sum=aParam*rowPhi*columnPhi
                    dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)= &
                      & dampingMatrix%elementMatrix%matrix(rowElementDOFIdx,columnElementDOFIdx)+ &
                      & sum*jacobianGaussWeight
                  ENDIF
                ENDDO !columnElementParameterIdx
              ENDDO !columnComponentIdx
            ENDIF !updateMatrices
            IF(updateResidual) THEN
              IF(esSpecification(3)==EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE) THEN
                sum=cParam*uValue**2*rowPhi
              ELSE
!!TODO: Handle floating point exceptions better
                IF((cParam*uValue)>20000.0_DP) THEN
                  localError="The value of "//TRIM(NumberToVString(cParam*uValue,"*",err,error))// &
                    & " is out of range for an exponential function."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                sum=bParam*EXP(cParam*uValue)*rowPhi
              ENDIF
              residualVector%elementResidual%vector(rowElementDOFIdx)=residualVector%elementResidual%vector(rowElementDOFIdx)+ &
                & sum*jacobianGaussWeight               
            ENDIF
            IF(updateSource) THEN
              sourceVector%elementVector%vector(rowElementDOFIdx)=sourceVector%elementVector%vector(rowElementDOFIdx)+ &
                & sourceValue*rowPhi*jacobianGaussWeight               
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
          DO rowElementParameterIdx=1,numberOfRowElementParameters(rowComponentIdx)
            rowElementDOFIdx=rowElementDOFIdx+1                    
            IF(updateMatrices) THEN
              !Loop over element columns
              columnElementDOFIdx=0
              DO columnComponentIdx=1,numberOfColsComponents
                DO columnElementParameterIdx=1,numberOfColumnElementParameters(columnComponentIdx)
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
            ENDIF !update matrices
            IF(updateResidual) THEN
              residualVector%elementResidual%vector(rowElementDOFIdx)=residualVector%elementResidual%vector(rowElementDOFIdx)* &
                & rowsInterpParameters%scaleFactors(rowElementParameterIdx,rowComponentIdx)
            ENDIF
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
                
    EXITS("Diffusion_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Diffusion_FiniteElementResidualEvaluate",err,error)
    EXITS("Diffusion_FiniteElementResidualEvaluate")
    RETURN 1
    
  END SUBROUTINE Diffusion_FiniteElementResidualEvaluate
 
  !
  !================================================================================================================================
  !

  !>Performs pre-solve operations for the a diffusion problem.
  SUBROUTINE Diffusion_PreSolve(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    INTEGER(INTG) :: pSpecification(3),outputType
    LOGICAL :: updateBoundaryConditions,updateMaterials
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem    
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolve",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
   
    updateMaterials = .FALSE.    
    updateBoundaryConditions = .TRUE.

    !IF(updateMaterials) THEN
    !  CALL Diffusion_PreSolveUpdateMaterials(solver,err,error,*999)
    !ENDIF
    
    !IF(updateBoundaryConditions) THEN
    !  CALL Diffusion_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    !ENDIF

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      !Do nothing ???
      CALL Diffusion_PreSolveUpdateAnalyticValues(solver,err,error,*999)
    CASE(PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      IF(outputType>=SOLVER_PROGRESS_OUTPUT) CALL WriteString(GENERAL_OUTPUT_TYPE,"ALE diffusion pre solve... ",err,error,*999)
      IF(.NOT.solver%dynamicSolver%ale) CALL FlagError("Mesh motion calculation not successful for ALE problem.",err,error,*999)
      !First update mesh and calculate boundary velocity values
      CALL Diffusion_PreSolveALEUpdateMesh(solver,err,error,*999)
      !Then apply both normal and moving mesh boundary conditions
      !CALL Diffusion_PreSolveUpdateBoundaryConditions(solver,err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  SUBROUTINE Diffusion_PreSolveUpdateBoundaryConditions(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("Diffusion_PreSolveUpdateBoundaryConditions",err,error,*999)

    !This routine previously set analytic BCs, but this has been moved. Needs rewriting to set
    CALL FlagError("Not implemented.",err,error,*999)
    
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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,boundaryConditionType,componentIdx,derivativeIdx,dimensionIdx, &
      & dynamicVariableType,equationsSetIdx,globalDerivativeIndex,globalDofIdx,localDofIdx,nodeIdx,numberOfComponents, &
      & numberOfDimensions,numberOfEquationsSets,numberOfNodes,numberOfNodeDerivatives,pSpecification(3)
    REAL(DP) :: A1,currentTime,D1,timeIncrement,normal(3),tangents(3,3),VALUE,X(3)
    REAL(DP), POINTER :: analyticParameters(:),geometricParameters(:),materialsParameters(:)
    TYPE(BoundaryConditionsType), POINTER :: boundaryConditions
    TYPE(BoundaryConditionsVariableType), POINTER :: boundaryConditionsVariable
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(DomainType), POINTER :: domain
    TYPE(DomainMappingType), POINTER :: domainMapping
    TYPE(DomainNodesType), POINTER :: domainNodes
    TYPE(DomainTopologyType), POINTER :: domainTopology
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: analyticField,dependentField,geometricField,materialsField,sourceField
    TYPE(FieldVariableType), POINTER :: dynamicVariable,geometricVariable,sourceVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolveUpdateAnalyticValues",err,error,*999)

    A1=0.4_DP
    D1=1.0_DP

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      !loop over all the equation sets and set the appropriate field variable type BCs and
      !the source field associated with each equation set
      CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
      DO equationsSetIdx=1,numberOfEquationsSets
        NULLIFY(equationsSet)
        CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
        NULLIFY(equationsAnalytic)
        CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
        IF(ASSOCIATED(equationsAnalytic)) THEN
          CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
          CALL EquationsSet_AnalyticTimeSet(equationsSet,currentTime,err,error,*999)
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
          CALL FieldVariable_VariableTypeGet(dynamicVariable,dynamicVariableType,err,error,*999)
          NULLIFY(domainMapping)
          CALL FieldVariable_DomainMappingGet(dynamicVariable,domainMapping,err,error,*999)
          NULLIFY(boundaryConditions)
          CALL SolverEquations_BoundaryConditionsGet(solverEquations,boundaryConditions,err,error,*999)
          NULLIFY(boundaryConditionsVariable)
          CALL BoundaryConditions_VariableGet(boundaryConditions,dynamicVariable,boundaryConditionsVariable,err,error,*999)
          CALL FieldVariable_NumberOfComponentsGet(dynamicVariable,numberOfComponents,err,error,*999)
          DO componentIdx=1,numberOfComponents
            CALL FieldVariable_ComponentInterpolationCheck(dynamicVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
              & err,error,*999)
            NULLIFY(domain)
            CALL FieldVariable_DomainGet(dynamicVariable,componentIdx,domain,err,error,*999)
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
                x(dimensionIdx)=geometricParameters(localDofIdx)
              ENDDO !dimensionIdx
              !Loop over the derivatives
              CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
              DO derivativeIdx=1,numberOfNodeDerivatives
                CALL DomainNodes_DerivativeGlobalIndexGet(domainNodes,derivativeIdx,nodeIdx,globalDerivativeIndex,err,error,*999)
                CALL Diffusion_AnalyticFunctionsEvaluate(equationsSet,analyticFunctionType,x,tangents,normal,currentTime, &
                  & dynamicVariableType,globalDerivativeIndex,componentIdx,analyticParameters,materialsParameters, &
                  & VALUE,err,error,*999)
                !Default to version 1 of each node derivative
                CALL FieldVariable_LocalNodeDOFGet(dynamicVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx,err,error,*999)
                CALL DomainMapping_LocalToGlobalGet(domainMapping,localDOFIdx,globalDOFIdx,err,error,*999)
                CALL FieldVariable_ParameterSetUpdateLocalDOF(dynamicVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,localDofIdx, &
                  & VALUE,err,error,*999)
                boundaryConditionType=boundaryConditionsVariable%DOFTypes(globalDofIdx)
                IF(boundaryConditionType==BOUNDARY_CONDITION_FIXED) THEN
                  CALL FieldVariable_ParameterSetUpdateLocalDOF(dynamicVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,VALUE, &
                    & err,error,*999)
                ENDIF
              ENDDO !derivativeIdx
            ENDDO !nodeIdx
          ENDDO !componentIdx
          CALL FieldVariable_ParameterSetUpdateStart(dynamicVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateStart(dynamicVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dynamicVariable,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
          CALL FieldVariable_ParameterSetUpdateFinish(dynamicVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
          IF(problem%specification(3)==PROBLEM_LINEAR_DIFFUSION_SUBTYPE) THEN
            !>Set the source field to a specified analytical function
            NULLIFY(sourceField)
            CALL EquationsSet_SourceFieldExists(equationsSet,sourceField,err,error,*999)
            IF(ASSOCIATED(sourceField)) THEN
              NULLIFY(sourceVariable)
              CALL Field_VariableGet(sourceField,FIELD_U_VARIABLE_TYPE,sourceVariable,err,error,*999)
              CALL FieldVariable_NumberOfComponentsGet(sourceVariable,numberOfComponents,err,error,*999)
              DO componentIdx=1,numberOfComponents
                CALL FieldVariable_ComponentInterpolationCheck(sourceVariable,componentIdx,FIELD_NODE_BASED_INTERPOLATION, &
                  & err,error,*999)
                NULLIFY(domain)
                CALL FieldVariable_DomainGet(dynamicVariable,componentIdx,domain,err,error,*999)
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
                    x(dimensionIdx)=geometricParameters(localDofIdx)
                  ENDDO !dimensionIdx
                  !Loop over the derivatives
                  CALL DomainNodes_NodeNumberOfDerivativesGet(domainNodes,nodeIdx,numberOfNodeDerivatives,err,error,*999)
                  DO derivativeIdx=1,numberOfNodeDerivatives
                    SELECT CASE(analyticFunctionType)
                    CASE(EQUATIONS_SET_DIFFUSION_EQUATION_THREE_DIM_1)
                      VALUE=-1*A1*EXP(-1*currentTime)*(X(1)*X(1)+X(2)*X(2)+X(3)*X(3)+6)
                    CASE DEFAULT
                      localError="The analytic function type of "//TRIM(NumberToVString(analyticFunctionType,"*", err,error))// &
                        & " is invalid."
                      CALL FlagError(localError,err,error,*999)
                    END SELECT
                    !Default to version 1 of each node derivative
                    CALL FieldVariable_LocalNodeDOFGet(sourceVariable,1,derivativeIdx,nodeIdx,componentIdx,localDOFIdx, &
                      & err,error,*999)
                    CALL FieldVariable_ParameterSetUpdateLocalDof(sourceVariable,FIELD_VALUES_SET_TYPE,localDOFIdx,VALUE, &
                      & err,error,*999)
                  ENDDO !derivativeIdx
                ENDDO !nodeIdx
              ENDDO !componentIdx
              CALL FieldVariable_ParameterSetUpdateStart(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
              CALL FieldVariable_ParameterSetUpdateFinish(sourceVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
            ENDIF
          ENDIF
          IF(ASSOCIATED(materialsField)) & 
            & CALL Field_ParameterSetDataRestore(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & materialsParameters,err,error,*999)
          IF(ASSOCIATED(analyticField)) &
            & CALL Field_ParameterSetDataRestore(analyticField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & analyticParameters,err,error,*999)
          CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_VALUES_SET_TYPE,geometricParameters,err,error,*999)
        ENDIF !Analytic field
      ENDDO !equationsSetIdx
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  SUBROUTINE Diffusion_PreSolveALEUpdateMesh(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dofNumber,esSpecification(3),inputOption,inputType,numberOfDOFsToPrint,numberOfDimensions, &
      & outputType,pSpecification(3),totalNumberOfDofs
    REAL(DP) :: alpha,currentTime,timeIncrement
    REAL(DP), POINTER :: inputData1(:),meshDisplacementValues(:)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: geometricField
    TYPE(FieldVariableType), POINTER :: geometricVariable
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverALEDiffusion
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolveALEUpdateMesh",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE)
      CALL ControlLoop_CurrentTimesGet(controlLoop,currentTime,timeIncrement,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)     
      NULLIFY(solverEquations)
      CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
      NULLIFY(solverMapping)
      CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
      NULLIFY(equationsSet)
      CALL SolverMapping_EquationsSetGet(solverMapping,1,equationsSet,err,error,*999)
      CALL EquationsSet_SpecificationGet(equationsSet,3,esSpecification,err,error,*999)
      SELECT CASE(esSpecification(3))
      CASE(EQUATIONS_SET_GENERALISED_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_DIFFUSION_SUBTYPE)
        ! do nothing ???
      CASE(EQUATIONS_SET_GENERALISED_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_LINEAR_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_QUADRATIC_SOURCE_ALE_DIFFUSION_SUBTYPE, &
        & EQUATIONS_SET_EXPONENTIAL_SOURCE_ALE_DIFFUSION_SUBTYPE)
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Diffusion update mesh ... ",err,error,*999)
        NULLIFY(geometricField)
        CALL EquationsSet_GeometricFieldGet(equationsSet,geometricField,err,error,*999)
        NULLIFY(geometricVariable)
        CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,geometricVariable,err,error,*999)
        !--- First, read mesh displacement values from file
        CALL FieldVariable_NumberOfComponentsGet(geometricVariable,numberOfDimensions,err,error,*999)
          
        inputType=42
        inputOption=2
        NULLIFY(inputData1)
        !CALL Field_ParameterSetDataGet(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,inputData1,err,error,*999)
        !CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputData1,numberOfDimensions,inputType,inputOption,currentTime)
        
        NULLIFY(meshDisplacementValues)
        CALL FieldVariable_ParameterSetDataGet(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,meshDisplacementValues, &
          & err,error,*999)
        IF(diagnostics1) THEN
          numberOfDOFsToPrint = SIZE(meshDisplacementValues,1)
          CALL WRITE_STRING_VECTOR(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDOFsToPrint,numberOfDOFsToPrint,numberOfDOFsToPrint,&
            & meshDisplacementValues,'(" meshDisplacementValues = ",3(X,E13.6))','3(3(X,E13.6))', &
            & err,error,*999)
        ENDIF

        !CALL FLUID_MECHANICS_IO_READ_DATA(SOLVER_LINEAR_TYPE,inputData1,numberOfDimensions,inputType,inputOption,currentTime)

        CALL FieldVariable_TotalNumberOfDOFsGet(geometricVariable,totalNumberOfDOFs,err,error,*999)

        !--- Second, update geometric field
        DO dofNumber=1,totalNumberOfDofs
          CALL FieldVariable_ParameterSetAddLocalDOF(geometricVariable,FIELD_VALUES_SET_TYPE,dofNumber, &
            & meshDisplacementValues(dofNumber),err,error,*999)
        ENDDO !dofNumber
        
        CALL FieldVariable_ParameterSetUpdateStart(geometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)
        CALL FieldVariable_ParameterSetUpdateFinish(geometricVariable,FIELD_VALUES_SET_TYPE,err,error,*999)

        !--- Third, use displacement values to calculate velocity values
        alpha=1.0_DP/timeIncrement
        CALL FieldVariable_ParameterSetsCopy(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,FIELD_MESH_VELOCITY_SET_TYPE, &
          & alpha,err,error,*999)
        
        CALL FieldVariable_ParameterSetDataRestore(geometricVariable,FIELD_MESH_DISPLACEMENT_SET_TYPE,meshDisplacementValues, &
          & err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(esSpecification(3),"*",err,error))// &
          & " is not valid for a diffusion equation type of a classical field problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Diffusion_PreSolveALEUpdateMesh")
    RETURN
999 ERRORSEXITS("Diffusion_PreSolveALEUpdateMesh",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveALEUpdateMesh
  
  !   
  !================================================================================================================================
  !
  
  SUBROUTINE Diffusion_PreSolveStoreCurrentSolution(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,numberOfComponentsDependentDiffusionOne,outputType,pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSetDiffusionOne
    TYPE(FieldType), POINTER :: dependentFieldDiffusionOne
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverDiffusionOne
    TYPE(SolverEquationsType), POINTER :: solverEquationsDiffusionOne
    TYPE(SolverMappingType), POINTER :: solverMappingDiffusionOne
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("Diffusion_PreSolveStoreCurrentSolution",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_COUPLED_DIFFUSION_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      IF(solverNumber==1) THEN
        !--- Get the dependent field of the diffusion-one equations
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) & 
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Store value diffusion-one dependent field at time, t ... ",err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverDiffusionOne)
        CALL Solvers_SolverGet(solvers,1,solverDiffusionOne,err,error,*999)
        NULLIFY(solverEquationsDiffusionOne)
        CALL Solver_SolverEquationsGet(solverDiffusionOne,solverEquationsDiffusionOne,err,error,*999)
        NULLIFY(solverMappingDiffusionOne)
        CALL SolverEquations_SolverMappingGet(solverEquationsDiffusionOne,solverMappingDiffusionOne,err,error,*999)
        NULLIFY(equationsSetDiffusionOne)
        CALL SolverMapping_EquationsSetGet(solverMappingDiffusionOne,1,equationsSetDiffusionOne,err,error,*999)
        NULLIFY(dependentFieldDiffusionOne)
        CALL EquationsSet_DependentFieldGet(equationsSetDiffusionOne,dependentFieldDiffusionOne,err,error,*999)
        CALL Field_NumberOfComponentsGet(dependentFieldDiffusionOne,FIELD_U_VARIABLE_TYPE, &
          & numberOfComponentsDependentDiffusionOne,err,error,*999)

        !--- Copy the current time value parameters set from diffusion-one's dependent field 
        DO componentIdx=1,numberOfComponentsDependentDiffusionOne
          CALL Field_ParametersToFieldParametersCopy(dependentFieldDiffusionOne, & 
            & FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,dependentFieldDiffusionOne, & 
            & FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,componentIdx,err,error,*999)
        ENDDO !componentIdx        

        !IF(diagnostics3) THEN
        !  NULLIFY(dummyValues2)
        !  CALL Field_ParameterSetDataGet(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !  numberOfDOFsToPrint = SIZE(dummyValues2,1)
        !  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDOFsToPrint,numberOfDOFsToPrint,numberOfDOFsToPrint, &
        !    & dummyValues2,'(" dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))', &
        !    & '4(4(X,E13.6))',err,error,*999)
        !  CALL Field_ParameterSetDataRestore(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !ENDIF
        
      ENDIF
    CASE(PROBLEM_COUPLED_DIFFUSION_ADVEC_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      IF(solverNumber==2) THEN
        !--- Get the dependent field of the diffusion equations
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE, &
          & "Store value of diffusion solution  (dependent field - V variable_type) at time, t ... ",err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverDiffusionOne)
        CALL Solvers_SolverGet(solvers,2,solverDiffusionOne,err,error,*999)
        NULLIFY(solverEquationsDiffusionOne)
        CALL Solver_SolverEquationsGet(solverDiffusionOne,solverEquationsDiffusionOne,err,error,*999)
        NULLIFY(solverMappingDiffusionOne)
        CALL SolverEquations_SolverMappingGet(solverEquationsDiffusionOne,solverMappingDiffusionOne,err,error,*999)
        NULLIFY(equationsSetDiffusionOne)
        CALL SolverMapping_EquationsSetGet(solverMappingDiffusionOne,1,equationsSetDiffusionOne,err,error,*999)
        NULLIFY(dependentFieldDiffusionOne)
        CALL EquationsSet_DependentFieldGet(equationsSetDiffusionOne,dependentFieldDiffusionOne,err,error,*999)
        CALL Field_NumberOfComponentsGet(dependentFieldDiffusionOne,FIELD_V_VARIABLE_TYPE, &
          & numberOfComponentsDependentDiffusionOne,err,error,*999)

        !--- Copy the current time value parameters set from diffusion-one's dependent field 
        DO componentIdx=1,numberOfComponentsDependentDiffusionOne
          CALL Field_ParametersToFieldParametersCopy(dependentFieldDiffusionOne, & 
            & FIELD_V_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,dependentFieldDiffusionOne, & 
            & FIELD_V_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,componentIdx,err,error,*999)
        ENDDO !componentIdx

        !IF(diagnostics3) THEN
        !  NULLIFY(dummyValues2)
        !  CALL Field_ParameterSetDataGet(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !  numberOfDOFsToPrint = SIZE(dummyValues2,1)
        !  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDOFsToPrint,numberOfDOFsToPrint,numberOfDOFsToPrint, &
        !    & dummyValues2,'(" dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))', &
        !    & '4(4(X,E13.6))',err,error,*999)
        !  CALL Field_ParameterSetDataRestore(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !ENDIF
        
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

    EXITS("Diffusion_PreSolveStoreCurrentSolution")
    RETURN
999 ERRORS("Diffusion_PreSolveStoreCurrentSolution",err,error)
    EXITS("Diffusion_PreSolveStoreCurrentSolution")
    RETURN 1
    
  END SUBROUTINE Diffusion_PreSolveStoreCurrentSolution
  
  !   
  !================================================================================================================================
  !
  SUBROUTINE Diffusion_PreSolveGetSourceValue(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: componentIdx,numberOfComponentsDependentDiffusionTwo,numberOfComponentsSourceDiffusionOne, &
      & outputType,pSpecification(3),solverNumber
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSetDiffusionOne,equationsSetDiffusionTwo
    TYPE(FieldType), POINTER :: dependentFieldDiffusionTwo,sourceFieldDiffusionOne
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverType), POINTER :: solverDiffusionOne,solverDiffusionTwo
    TYPE(SolverSType), POINTER :: solvers
    TYPE(SolverEquationsType), POINTER :: solverEquationsDiffusionOne,solverEquationsDiffusionTwo
    TYPE(SolverMappingType), POINTER :: solverMappingDiffusionOne,solverMappingDiffusionTwo
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PreSolveGetSourceValue",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)

    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_COUPLED_DIFFUSION_DIFFUSION_SUBTYPE)
      CALL Solver_GlobalNumberGet(solver,solverNumber,err,error,*999)
      CALL Solver_OutputTypeGet(solver,outputType,err,error,*999)
      IF(solverNumber==1) THEN
        !--- Get the dependent field of the diffusion_two equations
        IF(outputType>=SOLVER_PROGRESS_OUTPUT) &
          & CALL WriteString(GENERAL_OUTPUT_TYPE,"Update diffusion-one source field ... ",err,error,*999)
        NULLIFY(solvers)
        CALL Solver_SolversGet(solver,solvers,err,error,*999)
        NULLIFY(solverDiffusionTwo)
        CALL Solvers_SolverGet(solvers,2,solverDiffusionTwo,err,error,*999)
        NULLIFY(solverEquationsDiffusionTwo)
        CALL Solver_SolverEquationsGet(solverDiffusionTwo,solverEquationsDiffusionTwo,err,error,*999)
        NULLIFY(solverMappingDiffusionTwo)
        CALL SolverEquations_SolverMappingGet(solverEquationsDiffusionTwo,solverMappingDiffusionTwo,err,error,*999)
        NULLIFY(equationsSetDiffusionTwo)
        CALL SolverMapping_EquationsSetGet(solverMappingDiffusionTwo,1,equationsSetDiffusionTwo,err,error,*999)
        NULLIFY(dependentFieldDiffusionTwo)
        CALL EquationsSet_DependentFieldGet(equationsSetDiffusionTwo,dependentFieldDiffusionTwo,err,error,*999)
        CALL Field_NumberOfComponentsGet(dependentFieldDiffusionTwo,FIELD_U_VARIABLE_TYPE, &
          & numberOfComponentsDependentDiffusionTwo,err,error,*999)

        !--- Get the source field for the diffusion_one equations
        CALL Solvers_SolverGet(solvers,1,solverDiffusionOne,err,error,*999)
        NULLIFY(solverEquationsDiffusionOne)
        CALL Solver_SolverEquationsGet(solverDiffusionOne,solverEquationsDiffusionOne,err,error,*999)
        NULLIFY(solverMappingDiffusionOne)
        CALL SolverEquations_SolverMappingGet(solverEquationsDiffusionOne,solverMappingDiffusionOne,err,error,*999)
        NULLIFY(equationsSetDiffusionOne)
        CALL SolverMapping_EquationsSetGet(solverMappingDiffusionOne,1,equationsSetDiffusionOne,err,error,*999)
        NULLIFY(sourceFieldDiffusionOne)
        CALL EquationsSet_SourceFieldGet(equationsSetDiffusionOne,sourceFieldDiffusionOne,err,error,*999)
        CALL Field_NumberOfComponentsGet(sourceFieldDiffusionOne,FIELD_U_VARIABLE_TYPE, &
          & numberOfComponentsSourceDiffusionOne,err,error,*999)
        
        !--- Copy the result from diffusion-two's dependent field to diffusion-one's source field
        IF(numberOfComponentsSourceDiffusionOne/=numberOfComponentsDependentDiffusionTwo) THEN
          localError="Number of components of diffusion-two dependent field "// &
            & "is not consistent with diffusion-one-equation source field."
          CALL FlagError(localError,err,error,*999)
        END IF
        
        DO componentIdx=1,numberOfComponentsSourceDiffusionOne
          CALL Field_ParametersToFieldParametersCopy(dependentFieldDiffusionTwo,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
            & componentIdx,sourceFieldDiffusionOne,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,componentIdx,err,error,*999)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FIELDMESHDISPLACEMENTTYPE needs to be changed to appropriate type for this problem
        ENDDO !componentIdx

        !IF(diagnostics3) THEN
        !  NULLIFY(dummyValues2)
        !  CALL Field_ParameterSetDataGet(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !  numberOfDOFsToPrint = SIZE(dummyValues2,1)
        !  CALL WriteStringVector(DIAGNOSTIC_OUTPUT_TYPE,1,1,numberOfDOFsToPrint,numberOfDOFsToPrint,numberOfDOFsToPrint, &
        !    & dummyValues2,'(" dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE = ",4(X,E13.6))', &
        !    & '4(4(X,E13.6))',err,error,*999)
        !  CALL Field_ParameterSetDataRestore(dependentFieldFiniteElasticity,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
        !    & dummyValues2,err,error,*999)
        !ENDIF
        
      ENDIF
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT

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
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(ProblemType), POINTER :: problem
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PostSolve",err,error,*999)
    
    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
   
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE)
      !CALL Diffusion_PostSolveOuputData(solver,err,error,*999)
    CASE(PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
    CASE(PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
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
  SUBROUTINE Diffusion_PostSolveOutputData(solver,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: solver !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: analyticFunctionType,currentIteration,equationsSetIdx,inputIteration,numberOfEquationsSets, &
      & outputIteration,pSpecification(3)
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    CHARACTER(14) :: file,outputFile
    TYPE(ControlLoopType), POINTER :: controlLoop
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(EquationsSetAnalyticType), POINTER :: equationsAnalytic
    TYPE(FieldType), POINTER :: dependentField
    TYPE(ProblemType), POINTER :: problem
    TYPE(SolverEquationsType), POINTER :: solverEquations 
    TYPE(SolverMappingType), POINTER :: solverMapping 
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_PostSolveOutputData",err,error,*999)

    NULLIFY(controlLoop)
    CALL Solver_ControlLoopGet(solver,controlLoop,err,error,*999)
    NULLIFY(problem)
    CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    
    SELECT CASE(pSpecification(3))
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE)
      CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
        & outputIteration,inputIteration,err,error,*999)
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
            NULLIFY(equationsAnalytic)
            CALL EquationsSet_AnalyticExists(equationsSet,equationsAnalytic,err,error,*999)
            IF(ASSOCIATED(equationsAnalytic)) THEN
              CALL EquationsSet_AnalyticFunctionTypeGet(equationsSet,analyticFunctionType,err,error,*999)
              IF(analyticFunctionType==EQUATIONS_SET_DIFFUSION_EQUATION_TWO_DIM_1 .OR. &
                & analyticFunctionType==EQUATIONS_SET_LINEAR_SOURCE_DIFFUSION_EQUATION_THREE_DIM_1) THEN
                NULLIFY(dependentField)
                CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
                CALL AnalyticAnalysis_Output(dependentField,file,err,error,*999)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO !equationsSetIdx
    CASE(PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE)
      ! do nothing ???
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="Problem subtype "//TRIM(NumberToVString(pSpecification(3),"*",err,error))// &
        & " is not valid for a diffusion equation type of a classical field problem class."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    
    EXITS("Diffusion_PostSolveOutputData")
    RETURN
999 ERRORSEXITS("Diffusion_PostSolveOutputData",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PostSolveOutputData

  !
  !================================================================================================================================
  !

  !>Sets up the diffusion equations.
  SUBROUTINE Diffusion_ProblemSetup(problem,problemSetup,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: problemSetup !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype,pSpecification(3)
    TYPE(ControlLoopType), POINTER :: controlLoop,controlLoopRoot
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Diffusion_ProblemSetup",err,error,*999)

    CALL Problem_SpecificationGet(problem,3,pSpecification,err,error,*999)
    problemSubtype=pSpecification(3)

    SELECT CASE(problemSubtype)
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE)
      !ok
    CASE DEFAULT
      localError="The problem subtype of "//TRIM(NumberToVString(problemSubtype,"*",err,error))// &
        & " does not equal a linear diffusion equation subtype."
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
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_CONTROL_TYPE)
      SELECT CASE(problemSetup%actionType)
      CASE(PROBLEM_SETUP_START_ACTION)
        !Set up a time control loop
        NULLIFY(controlLoop)
        CALL ControlLoop_CreateStart(problem,controlLoop,err,error,*999)
        CALL ControlLoop_TypeSet(controlLoop,CONTROL_TIME_LOOP_TYPE,err,error,*999)
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
          & " is invalid for a linear diffusion equation."
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
        !Set the solver to be a first order dynamic solver
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        CALL Solver_TypeSet(solver,SOLVER_DYNAMIC_TYPE,err,error,*999)
        CALL Solver_LabelSet(solver,"Dynamic solver",err,error,*999)
        CALL Solver_DynamicOrderSet(solver,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
        IF(problemSubtype==PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE.OR. &
          & problemSubtype==PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE) THEN
          CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
        ELSE
          CALL Solver_DynamicLinearityTypeSet(solver,SOLVER_DYNAMIC_LINEAR,err,error,*999)
        ENDIF
        !Set solver defaults
        CALL Solver_DynamicDegreeSet(solver,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
        CALL Solver_DynamicSchemeSet(solver,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
        CALL Solver_LibraryTypeSet(solver,SOLVER_CMISS_LIBRARY,err,error,*999)
      CASE(PROBLEM_SETUP_FINISH_ACTION)
        !Get the solvers
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        !Finish the solvers creation
        CALL Solvers_CreateFinish(solvers,err,error,*999)
      CASE DEFAULT
        localError="The action type of "//TRIM(NumberToVString(problemSetup%actionType,"*",err,error))// &
          & " for a setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
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
        IF(problemSubtype==PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE.OR. &
          & problemSubtype==PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE) THEN
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_NONLINEAR,err,error,*999)
        ELSE          
          CALL SolverEquations_LinearityTypeSet(solverEquations,SOLVER_EQUATIONS_LINEAR,err,error,*999)
        ENDIF
        CALL SolverEquations_TimeDependenceTypeSet(solverEquations,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
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
          & " is invalid for a linear diffusion equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE DEFAULT
      localError="The setup type of "//TRIM(NumberToVString(problemSetup%setupType,"*",err,error))// &
        & " is invalid for a linear diffusion equation."
      CALL FlagError(localError,err,error,*999)
    END SELECT
       
    EXITS("Diffusion_ProblemSetup")
    RETURN
999 ERRORSEXITS("Diffusion_ProblemSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_ProblemSetup
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a diffusion equation problem.
  SUBROUTINE Diffusion_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: problemSubtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Diffusion_ProblemSpecificationSet",err,error,*999)

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
    CASE(PROBLEM_LINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_DIFFUSION_SUBTYPE, &
      & PROBLEM_LINEAR_ALE_DIFFUSION_SUBTYPE, &
      & PROBLEM_NONLINEAR_ALE_DIFFUSION_SUBTYPE)
!!TODO: WHY DO WE HAVE ALE PROBLEMS AND EQUATIONS IN DIFFUSION PROBLEMS??? THERE IS NO VELOCITY????
      !All ok
    CASE DEFAULT
      localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
        & " is not valid for a diffusion type of a classical field problem."
      CALL FlagError(localError,err,error,*999)
    END SELECT
    ALLOCATE(problem%specification(3),stat=err)
    IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
    problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_DIFFUSION_EQUATION_TYPE,problemSubtype]

    EXITS("Diffusion_ProblemSpecificationSet")
    RETURN
999 ERRORS("Diffusion_ProblemSpecificationSet",err,error)
    EXITS("Diffusion_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Diffusion_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE Diffusion_PostLoop(controlLoop,err,error,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: controlLoop !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: currentIteration,equationsSetIdx,inputIteration,loopType,numberOfEquationsSets,outputIteration,outputType, &
      & regionUserNumber
    REAL(DP) :: currentTime,startTime,stopTime,timeIncrement
    TYPE(ControlLoopType), POINTER :: parentLoop
    TYPE(ControlLoopTimeType), POINTER :: timeLoop,timeLoopParent
    TYPE(EquationsSetType), POINTER :: equationsSet
    TYPE(FieldType), POINTER :: dependentField
    TYPE(ProblemType), POINTER :: problem
    TYPE(RegionType), POINTER :: dependentRegion   
    TYPE(SolverType), POINTER :: solver
    TYPE(SolverEquationsType), POINTER :: solverEquations
    TYPE(SolverMappingType), POINTER :: solverMapping
    TYPE(SolversType), POINTER :: solvers
    TYPE(VARYING_STRING) :: filename,localError,method

    ENTERS("Diffusion_PostLoop",err,error,*999)

    CALL ControlLoop_OutputTypeGet(controlLoop,outputType,err,error,*999)
    IF(outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
!!TODO: 1) WHY IS OUTPUT PRELOOP RATHER THAN POST LOOP 2) WHY IS OUTPUT TIME USED TO DETERMINE DEPENDENT OUTPUT???
      CALL ControlLoop_TypeGet(controlLoop,loopType,err,error,*999)
      SELECT CASE(loopType)
      CASE(CONTROL_SIMPLE_TYPE)
        !do nothing
      CASE(CONTROL_FIXED_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_TIME_LOOP_TYPE)
        CALL ControlLoop_CurrentTimeInformationGet(controlLoop,currentTime,timeIncrement,startTime,stopTime,currentIteration, &
          & outputIteration,inputIteration,err,error,*999)
        !Export the dependent field for this time step
        NULLIFY(problem)
        CALL ControlLoop_ProblemGet(controlLoop,problem,err,error,*999)
        NULLIFY(solvers)
        CALL ControlLoop_SolversGet(controlLoop,solvers,err,error,*999)
        NULLIFY(solver)
        CALL Solvers_SolverGet(solvers,1,solver,err,error,*999)
        NULLIFY(solverEquations)
        CALL Solver_SolverEquationsGet(solver,solverEquations,err,error,*999)
        NULLIFY(solverMapping)
        CALL SolverEquations_SolverMappingGet(solverEquations,solverMapping,err,error,*999)
        !Loop over the equations sets associated with the solver
        CALL SolverMapping_NumberOfEquationsSetsGet(solverMapping,numberOfEquationsSets,err,error,*999)
        DO equationsSetIdx=1,numberOfEquationsSets
          NULLIFY(equationsSet)
          CALL SolverMapping_EquationsSetGet(solverMapping,equationsSetIdx,equationsSet,err,error,*999)
          NULLIFY(dependentField)
          CALL EquationsSet_DependentFieldGet(equationsSet,dependentField,err,error,*999)
          NULLIFY(dependentRegion)
          CALL Field_RegionGet(dependentField,dependentRegion,err,error,*999)
          CALL Region_UserNumberGet(dependentRegion,regionUserNumber,err,error,*999)
          filename="Time_"//TRIM(NumberToVString(regionUserNumber,"*",err,error))// &
            & "_"//TRIM(NumberToVString(currentIteration,"*",err,error))
          method="FORTRAN"
          IF(outputIteration/=0) THEN
            IF(MOD(currentIteration,outputIteration)==0) THEN
              CALL FIELD_IO_NODES_EXPORT(dependentRegion%fields,filename,method,err,error,*999)
            ENDIF
          ENDIF
        ENDDO !equationsSetIdx        
      CASE(CONTROL_WHILE_LOOP_TYPE)
        !do nothing
      CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
        !do nothing
      CASE DEFAULT
        localError="The control loop type of "//TRIM(NumberToVString(loopType,"*",err,error))//" is invalid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ENDIF

    EXITS("Diffusion_PostLoop")
    RETURN
999 ERRORSEXITS("Diffusion_PostLoop",err,error)
    RETURN 1
    
  END SUBROUTINE Diffusion_PostLoop

  !
  !================================================================================================================================
  !

END MODULE DiffusionEquationsRoutines 
