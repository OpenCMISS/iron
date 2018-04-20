!> \file
!> \author David Ladd
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

!>This module handles all Burgers equation routines.
MODULE BURGERS_EQUATION_ROUTINES

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
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FIELD_IO_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
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

  PUBLIC BURGERS_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC Burgers_EquationsSetSolutionMethodSet

  PUBLIC Burgers_EquationsSetSpecificationSet

  PUBLIC Burgers_FiniteElementJacobianEvaluate

  PUBLIC Burgers_FiniteElementResidualEvaluate

  PUBLIC Burgers_ProblemSpecificationSet

  PUBLIC BURGERS_EQUATION_PROBLEM_SETUP

  PUBLIC BURGERS_EQUATION_PRE_SOLVE,BURGERS_EQUATION_POST_SOLVE

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  !Calculates a one-dimensional dynamic solution to the burgers equation
  SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,version_idx
    REAL(DP) :: VALUE,X(3),INITIAL_VALUE
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:),GEOMETRIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX,ANALYTIC_FUNCTION_TYPE
    !THESE ARE TEMPORARY VARIABLES - they need to be replace by constant field values and the current simulation time
    REAL(DP) :: TIME,NORMAL(3),TANGENTS(3,3)
    !CURRENT_TIME = 1.2_DP

    ENTERS("Burgers_BoundaryConditionsAnalyticCalculate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN
            ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            NULLIFY(GEOMETRIC_PARAMETERS)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)
            NULLIFY(ANALYTIC_VARIABLE)
            NULLIFY(ANALYTIC_PARAMETERS)
            IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
              CALL Field_VariableGet(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & ANALYTIC_PARAMETERS,err,error,*999)
            ENDIF
            NULLIFY(materialsField)
            NULLIFY(MATERIALS_VARIABLE)
            NULLIFY(MATERIALS_PARAMETERS)
            IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
              materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
              CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & MATERIALS_PARAMETERS,err,error,*999)
            ENDIF
            ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
            TIME=EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                variable_type=dependentField%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL Field_ParameterSetEnsureCreated(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                      DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                      IF(ASSOCIATED(DOMAIN)) THEN
                        IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
                              DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                !Default to version 1 of each node derivative
                                local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                  & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                              ENDDO !dim_idx
                              !Loop over the derivatives
                              DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX
                                CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE, &
                                  & X,TANGENTS,NORMAL,0.0_DP,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                  & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,INITIAL_VALUE,err,error,*999)
                                CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE, &
                                  & X,TANGENTS,NORMAL,TIME,variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                  & ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,err,error,*999)
                                DO version_idx=1,DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%numberOfVersions
                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(version_idx)
                                  CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                  IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                    IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                      !If we are a boundary node then set the analytic value on the boundary
                                      CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                        & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                    ELSE
                                      !Set the initial condition.
                                      CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                        & FIELD_VALUES_SET_TYPE,local_ny,INITIAL_VALUE,err,error,*999)
                                    ENDIF
                                  ENDIF
                                ENDDO !version_idx
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ELSE
                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Domain topology is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Domain is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Boundary conditions are not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Burgers_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORSEXITS("Burgers_BoundaryConditionsAnalyticCalculate",err,error)
    RETURN 1

  END SUBROUTINE Burgers_BoundaryConditionsAnalyticCalculate


  !
  !================================================================================================================================
  !
  !>Evaluate the analytic solutions for a Burgers equation
  SUBROUTINE Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,X, &
    & TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
    & VALUE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: EQUATIONS_SET !<The equations set to evaluate
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: X(:) !<X(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: TANGENTS(:,:) !<TANGENTS(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: NORMAL(:) !<NORMAL(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: TIME !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: GLOBAL_DERIVATIVE !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: ANALYTIC_PARAMETERS(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: MATERIALS_PARAMETERS(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analytic function value.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,D_PARAM,E_PARAM,X0_PARAM
    INTEGER(INTG) :: EQUATIONS_SUBTYPE
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_AnalyticFunctionsEvaluate",err,error,*999)

    IF(.NOT.ASSOCIATED(EQUATIONS_SET)) THEN
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ELSE
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      ELSE
        EQUATIONS_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      END IF
    END IF
    SELECT CASE(EQUATIONS_SUBTYPE)
    CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
        !For del[u]/del[t] + u.(del[u]/del[x]) = nu.(del^2[u]/del[x]^2)
        !u(x,t)=1-tanh(x-x_0-t)/(2.nu))   with BCs,
        !u(0,t) = 2, u_{n} = 2.u_{n-1} - u_{n-2}
        !see http://www.cfd-online.com/Wiki/Burgers_equation
        !OpenCMISS has del[u]/del[t] + K.(del^2[u]/del[x]^2) + u.(del[u]/del[x]) = 0,
        !u(x,t)= 1 - tanh(x-x_0 - t)/(2.K)
        B_PARAM=MATERIALS_PARAMETERS(1)  !nu
        X0_PARAM=ANALYTIC_PARAMETERS(1)  !x_0
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=1.0_DP - TANH((X(1)-X0_PARAM-TIME)/(2.0_DP*B_PARAM))
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            CALL FlagError("Not implemented.",err,error,*999)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
          & " is invalid for a Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
      !a.del u/del t + b.del^2 u/del x^2 + c.u.del u/del x = 0
      A_PARAM=MATERIALS_PARAMETERS(1)
      B_PARAM=MATERIALS_PARAMETERS(2)
      C_PARAM=MATERIALS_PARAMETERS(3)
      SELECT CASE(ANALYTIC_FUNCTION_TYPE)
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1)
        !Analytic solution is u(x,t)=(D+a.x)/(E+c.t)
        D_PARAM = ANALYTIC_PARAMETERS(1)
        E_PARAM = ANALYTIC_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=(D_PARAM+A_PARAM*X(1))/(E_PARAM+C_PARAM*TIME)
          CASE(GLOBAL_DERIV_S1)
            VALUE=D_PARAM/(E_PARAM+C_PARAM*TIME)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            VALUE=0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
        !Analytic_solution=a.D+2.b/c(x-c.D.t+E)
        D_PARAM = ANALYTIC_PARAMETERS(1)
        E_PARAM = ANALYTIC_PARAMETERS(2)
        SELECT CASE(VARIABLE_TYPE)
        CASE(FIELD_U_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=A_PARAM*D_PARAM+2.0_DP*B_PARAM/(C_PARAM*(X(1)-C_PARAM*D_PARAM*TIME+E_PARAM))
          CASE(GLOBAL_DERIV_S1)
            VALUE=-2.0_DP*B_PARAM/(C_PARAM*(X(1)-C_PARAM*D_PARAM*TIME+E_PARAM)**2)
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(FIELD_DELUDELN_VARIABLE_TYPE)
          SELECT CASE(GLOBAL_DERIVATIVE)
          CASE(NO_GLOBAL_DERIV)
            VALUE=0.0_DP
          CASE(GLOBAL_DERIV_S1)
            VALUE=0.0_DP
          CASE DEFAULT
            localError="The global derivative index of "//TRIM(NUMBER_TO_VSTRING(GLOBAL_DERIVATIVE,"*",err,error))// &
              & " is invalid."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The variable type of "//TRIM(NUMBER_TO_VSTRING(VARIABLE_TYPE,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The analytic function type of "// &
          & TRIM(NUMBER_TO_VSTRING(ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
          & " is invalid for a generalised Burgers equation."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE DEFAULT
      localError="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SUBTYPE,"*",err,error))// &
        & " is invalid."
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
  SUBROUTINE Burgers_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_EquationsSetSolutionMethodSet",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
        & EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
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
          localError="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Burgers equation type of an classical field equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Burgers_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers equation set.", &
          & err,error,*999)
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
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_BURGERS_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

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
  SUBROUTINE BURGERS_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a BURGERS equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
      & NUMBER_OF_DEPENDENT_COMPONENTS,NUMBER_OF_GEOMETRIC_COMPONENTS,NUMBER_OF_MATERIALS_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_ANALYTIC
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,geometricField
    TYPE(VARYING_STRING) :: localError

    ENTERS("BURGERS_EQUATION_EQUATIONS_SET_SETUP",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE, &
        EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Burgers_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear burgers equation."
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
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_GEOMETRIC_COMPONENTS,err,error,*999)
                NUMBER_OF_DEPENDENT_COMPONENTS=NUMBER_OF_GEOMETRIC_COMPONENTS
              ELSE
                NUMBER_OF_DEPENDENT_COMPONENTS=1
              ENDIF
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_U_VARIABLE_TYPE,NUMBER_OF_DEPENDENT_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                & FIELD_DELUDELN_VARIABLE_TYPE,NUMBER_OF_DEPENDENT_COMPONENTS,err,error,*999)
              !Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ENDDO !component_idx
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_U_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !component_idx
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
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
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
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_GEOMETRIC_COMPONENTS,err,error,*999)
                NUMBER_OF_DEPENDENT_COMPONENTS=NUMBER_OF_GEOMETRIC_COMPONENTS
              ELSE
                NUMBER_OF_DEPENDENT_COMPONENTS=1
              ENDIF
              IF(NUMBER_OF_DEPENDENT_COMPONENTS==1) THEN
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
              ELSE
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
              ENDIF
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DEPENDENT_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_DEPENDENT_COMPONENTS,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                DO component_idx=1,NUMBER_OF_DEPENDENT_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear Burgers equation."
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
              IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                !Not an inviscid Burgers equation
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
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component
                    !i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=1
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                    !3 materials field components
                    !i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=3
                  CASE DEFAULT
                    localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  !Set the number of materials components
                  CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                  !Default the materials components to the 1st geometric component interpolation setup with constant interpolation
                  CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  DO component_idx=1,NUMBER_OF_MATERIALS_COMPONENTS
                    CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                    CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  ENDDO !component_idx
                  !Default the field scaling to that of the geometric field
                  CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
                ELSE
                  !Check the user specified field
                  CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
                  CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                  CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
                  CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                  CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component
                    !i.e., k = viscosity*(-1) in du/dt + k*(d^2u/dx^2)+ u*(du/dx) = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=1
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                      & err,error,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                    !3 materials field components
                    !i.e., a.du/dt + b.(d^2u/dx^2) + c.u*(du/dx)  = 0
                    NUMBER_OF_MATERIALS_COMPONENTS=3
                    CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                      & err,error,*999)
                  CASE DEFAULT
                    localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                  CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                !Not an inviscid Burgers equation
                IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                  !Finish creating the materials field
                  CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                  !Set the default values for the materials field
                  SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    !1 materials field component. Default to
                    !du/dt - d^2u/dx^2 + u*(du/dx) = 0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,-1.0_DP,err,error,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                    !3 materials field components. Default to
                    !du/dt - d^2u/dx^2 + u*(du/dx) = 0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,2,-1.0_DP,err,error,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,3,1.0_DP,err,error,*999)
                  CASE DEFAULT
                    localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                      & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                      & " is invalid for a nonlinear Burgers equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
             ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! S o u r c e   f i e l d
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear Burgers equation."
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
                dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                IF(ASSOCIATED(dependentField)) THEN
                  EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
                  IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
                    IF(EQUATIONS_MATERIALS%MATERIALS_FINISHED) THEN
                      geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                      IF(ASSOCIATED(geometricField)) THEN
                        CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE, &
                          & NUMBER_OF_GEOMETRIC_COMPONENTS,err,error,*999)
                        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
                        CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1)
                            !Check that domain is 1D
                            IF(NUMBER_OF_GEOMETRIC_COMPONENTS/=1) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_COMPONENTS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                & " requires that there be 1 geometric dimension."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Check the materials values are constant
                            CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                              & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_BURGERS_EQUATION_ONE_DIM_1
                            NUMBER_OF_ANALYTIC_COMPONENTS=1
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a Burgers equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                          SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                          CASE(EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_1, &
                            & EQUATIONS_SET_GENERALISED_BURGERS_EQUATION_ONE_DIM_2)
                            !Check that domain is 1D
                            IF(NUMBER_OF_GEOMETRIC_COMPONENTS/=1) THEN
                              localError="The number of geometric dimensions of "// &
                                & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_GEOMETRIC_COMPONENTS,"*",err,error))// &
                                & " is invalid. The analytic function type of "// &
                                & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
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
                            !Set analytic function type
                            EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE
                            NUMBER_OF_ANALYTIC_COMPONENTS=2
                          CASE DEFAULT
                            localError="The specified analytic function type of "// &
                              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                              & " is invalid for a generalised Burgers equation."
                            CALL FlagError(localError,err,error,*999)
                          END SELECT
                        CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE DEFAULT
                          localError="The equation set subtype of "// &
                            & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                            & " is invalid for an analytical nonlinear Burgers equation."
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
                  CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
                    !Default the analytic parameter value to 0.0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,0.0_DP,err,error,*999)
                  CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
                    !Default the analytic parameter values to 1.0
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                    CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_ANALYTIC%ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                      & FIELD_VALUES_SET_TYPE,2,1.0_DP,err,error,*999)
                  CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
                    CALL FlagError("Not implemented.",err,error,*999)
                  CASE DEFAULT
                    localError="The equation set subtype of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                      & " is invalid for an analytical nonlinear Burgers equation."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear Burgers equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        !-----------------------------------------------------------------
        ! E q u a t i o n s    t y p e
        !-----------------------------------------------------------------
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE (EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_BURGERS_SUBTYPE,EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE,EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
                CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                  CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.FALSE.,err,error,*999)
                ELSE
                  CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                ENDIF
                CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                !Set up matrix storage and structure
                IF(equations%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                  !Set up lumping
                  IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                    CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                      & [EQUATIONS_MATRIX_LUMPED],err,error,*999)
                  ELSE
                    CALL EquationsMatrices_DynamicLumpingTypeSet(vectorMatrices, &
                      & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                  ENDIF
                  SELECT CASE(equations%sparsityType)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        [EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                    ELSE
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                      & err,error,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        & [EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                    ELSE
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE],err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                      & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                    CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE, &
                      & err,error,*999)
                  CASE DEFAULT
                    localError="The equations matrices sparsity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(equations%sparsityType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ELSE
                  SELECT CASE(equations%sparsityType)
                  CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                    ELSE
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE, &
                      & err,error,*999)
                  CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                    IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                        & err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                    ELSE
                      CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                        & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                        & err,error,*999)
                      CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                        [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                    ENDIF
                    CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                      & DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                    CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                      EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                  CASE DEFAULT
                    localError="The equations matrices sparsity type of "// &
                      & TRIM(NUMBER_TO_VSTRING(equations%sparsityType,"*",err,error))//" is invalid."
                    CALL FlagError(localError,err,error,*999)
                  END SELECT
                ENDIF
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & err,error,*999)
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
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
            SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
            CASE(EQUATIONS_SET_SETUP_START_ACTION)
              IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
                CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_LinearityTypeSet(equations,EQUATIONS_NONLINEAR,err,error,*999)
                CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
              ELSE
                CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
              ENDIF
            CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Finish the equations
                CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
                CALL Equations_CreateFinish(equations,err,error,*999)
                NULLIFY(vectorEquations)
                CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
                !Create the equations mapping.
                CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
                CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
                CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
                CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
                !Create the equations matrices
                CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE, &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices, &
                    & [MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                  CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                  CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices, &
                    & MATRIX_COMPRESSED_ROW_STORAGE_TYPE,err,error,*999)
                  CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices, &
                    & EQUATIONS_MATRIX_FEM_STRUCTURE,err,error,*999)
                CASE DEFAULT
                  localError="The equations matrices sparsity type of "// &
                    & TRIM(NUMBER_TO_VSTRING(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(localError,err,error,*999)
                END SELECT
                ! Use the analytic Jacobian calculation
                CALL EquationsMatrices_JacobianTypesSet(vectorMatrices,[EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED], &
                  & err,error,*999)
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
                localError="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a nonlinear Burgers equation."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The equation set subtype of "// &
              & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is invalid for an analytical nonlinear Burgers equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a nonlinear Burgers equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not a nonlinear Burgers equation subtype."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("BURGERS_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_EQUATIONS_SET_SETUP",err,error)
    RETURN 1

  END SUBROUTINE BURGERS_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the BURGERS problem pre-solve.
  SUBROUTINE BURGERS_EQUATION_PRE_SOLVE(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(DYNAMIC_SOLVER_TYPE), POINTER :: DYNAMIC_SOLVER
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

    ENTERS("BURGERS_EQUATION_PRE_SOLVE",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification array is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
            END IF
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
              DYNAMIC_SOLVER=>SOLVER%DYNAMIC_SOLVER
              IF(ASSOCIATED(DYNAMIC_SOLVER)) THEN
                IF(DYNAMIC_SOLVER%SOLVER_INITIALISED) &
                  & CALL Burgers_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER,err,error,*999)
              ELSE
                CALL FlagError("Solver dynamic solver is not associated.",err,error,*999)
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Burgers equation type of a fluid mechanics problem class."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF

    EXITS("BURGERS_EQUATION_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_PRE_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE BURGERS_EQUATION_PRE_SOLVE


  !
  !================================================================================================================================
  !
  !updates the boundary conditions and source term to the required analytic values
  SUBROUTINE Burgers_PreSolveUpdateAnalyticValues(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: ANALYTIC_VARIABLE,FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(EquationsType), POINTER :: EQUATIONS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(VARYING_STRING) :: localError
    REAL(DP), POINTER :: ANALYTIC_PARAMETERS(:),GEOMETRIC_PARAMETERS(:),MATERIALS_PARAMETERS(:)
    INTEGER(INTG) :: NUMBER_OF_DIMENSIONS,BOUNDARY_CONDITION_CHECK_VARIABLE

    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    REAL(DP) :: NORMAL(3),TANGENTS(3,3),VALUE,X(3) !<The value to add
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,variable_idx,eqnset_idx
    INTEGER(INTG) :: VARIABLE_TYPE !<The field variable type to add \see FIELD_ROUTINES_VariableTypes,FIELD_ROUTINES
    INTEGER(INTG) :: ANALYTIC_FUNCTION_TYPE
    INTEGER(INTG) :: GLOBAL_DERIV_INDEX

    ENTERS("Burgers_PreSolveUpdateAnalyticValues",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
            SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
            IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
              !Loop over all the equation sets and set the appropriate field variable type BCs and
              !the source field associated with each equation set
              DO eqnset_idx=1,SOLVER_equations%SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                SOLVER_MAPPING=>SOLVER_equations%SOLVER_MAPPING
                EQUATIONS=>SOLVER_MAPPING%EQUATIONS_SET_TO_SOLVER_MAP(eqnset_idx)%EQUATIONS
                IF(ASSOCIATED(EQUATIONS)) THEN
                  EQUATIONS_SET=>equations%equationsSet
                  IF(ASSOCIATED(EQUATIONS_SET)) THEN
                    IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                      dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
                      IF(ASSOCIATED(dependentField)) THEN
                        geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                        IF(ASSOCIATED(geometricField)) THEN
                          ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
                          CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,&
                            & NUMBER_OF_DIMENSIONS,err,error,*999)
                          NULLIFY(GEOMETRIC_VARIABLE)
                          NULLIFY(GEOMETRIC_PARAMETERS)
                          CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
                          CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,&
                            & GEOMETRIC_PARAMETERS,err,error,*999)
                          EQUATIONS_SET%ANALYTIC%ANALYTIC_USER_PARAMS(1)=CURRENT_TIME
                          NULLIFY(ANALYTIC_VARIABLE)
                          NULLIFY(ANALYTIC_PARAMETERS)
                          IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                            CALL Field_VariableGet(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,ANALYTIC_VARIABLE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & ANALYTIC_PARAMETERS,err,error,*999)
                          ENDIF
                          NULLIFY(materialsField)
                          NULLIFY(MATERIALS_VARIABLE)
                          NULLIFY(MATERIALS_PARAMETERS)
                          IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                            materialsField=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                            CALL Field_VariableGet(materialsField,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_DATA_GET(materialsField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                              & MATERIALS_PARAMETERS,err,error,*999)
                          ENDIF
                          DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                            variable_type=dependentField%VARIABLES(variable_idx)%VARIABLE_TYPE
                            FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                            IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                              DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                                IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE== &
                                  & FIELD_NODE_BASED_INTERPOLATION) THEN
                                DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                                IF(ASSOCIATED(DOMAIN)) THEN
                                  IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                                    DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                                    IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                      !Loop over the local nodes excluding the ghosts.
                                      DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
!!TODO \todo We should interpolate the geometric field here and the node position.
                                        DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                          !Default to version 1 of each node derivative
                                          local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                          X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                        ENDDO !dim_idx
                                        !Loop over the derivatives
                                        DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                          ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE
                                          GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)% &
                                            & GLOBAL_DERIVATIVE_INDEX
                                          CALL Burgers_AnalyticFunctionsEvaluate(EQUATIONS_SET, &
                                            & ANALYTIC_FUNCTION_TYPE,X,TANGENTS,NORMAL,CURRENT_TIME,variable_type, &
                                            & GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS, &
                                            & VALUE,err,error,*999)
                                          !Default to version 1 of each node derivative
                                          local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                            & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                          CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField,variable_type, &
                                            & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                          CALL BOUNDARY_CONDITIONS_VARIABLE_GET(SOLVER_equations%BOUNDARY_CONDITIONS, &
                                            & FIELD_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE,err,error,*999)
                                          IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                                            BOUNDARY_CONDITION_CHECK_VARIABLE=BOUNDARY_CONDITIONS_VARIABLE% &
                                              & CONDITION_TYPES(local_ny)
                                            IF(BOUNDARY_CONDITION_CHECK_VARIABLE==BOUNDARY_CONDITION_FIXED) THEN
                                              CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(dependentField, &
                                                & variable_type,FIELD_VALUES_SET_TYPE,local_ny, &
                                                & VALUE,err,error,*999)
                                            ENDIF
                                          ELSE
                                            CALL FlagError("Boundary conditions variable is not associated",err,error,*999)
                                          ENDIF
                                        ENDDO !deriv_idx
                                      ENDDO !node_idx
                                    ELSE
                                      CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                                    ENDIF
                                  ELSE
                                    CALL FlagError("Domain topology is not associated.",err,error,*999)
                                  ENDIF
                                ELSE
                                  CALL FlagError("Domain is not associated.",err,error,*999)
                                ENDIF
                              ELSE
                                CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                              ENDIF
                            ENDDO !component_idx
                            CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type, &
                              & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type, &
                              & FIELD_VALUES_SET_TYPE,err,error,*999)
                            CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type, &
                              & FIELD_VALUES_SET_TYPE,err,error,*999)
                          ELSE
                            CALL FlagError("Field variable is not associated.",err,error,*999)
                          ENDIF

                        ENDDO !variable_idx
                          CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,&
                            & FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS,err,error,*999)
                        ELSE
                          CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
                      ENDIF
                    ENDIF
                  ELSE
                    CALL FlagError("Equations set is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations are not associated.",err,error,*999)
                END IF
                CALL FIELD_PARAMETER_SET_UPDATE_START(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
                CALL FIELD_PARAMETER_SET_UPDATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,err,error,*999)
              ENDDO !eqnset_idx
            ELSE
              CALL FlagError("Solver equations are not associated.",err,error,*999)
            END IF
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a BURGERS equation type of a fluid mechanics problem class."
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

    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN
999 ERRORS("Burgers_PreSolveUpdateAnalyticValues",err,error)
    EXITS("Burgers_PreSolveUpdateAnalyticValues")
    RETURN 1

  END SUBROUTINE Burgers_PreSolveUpdateAnalyticValues


  !
  !================================================================================================================================
  !
  SUBROUTINE Burgers_PreSolveStoreCurrentSolution(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solvers
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_PreSolveStoreCurrentSolution",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN

      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
              ! do nothing ???
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Burgers equation type of a fluid mechanics problem class."
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
  SUBROUTINE BURGERS_EQUATION_POST_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("BURGERS_EQUATION_POST_SOLVE",err,error,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification array is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STATIC_BURGERS_SUBTYPE,PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
            CALL BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*999)
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a BURGERS type of a fluid mechanics problem class."
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

    EXITS("BURGERS_EQUATION_POST_SOLVE")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_POST_SOLVE",err,error)
    RETURN 1

  END SUBROUTINE BURGERS_EQUATION_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Output data post solve
  SUBROUTINE BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA(SOLVER,err,error,*)

    !Argument variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set
    TYPE(FIELDS_TYPE), POINTER :: Fields
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS  !<A pointer to the solver equations
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING !<A pointer to the solver mapping
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError,METHOD,FILENAME
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT,START_TIME,STOP_TIME
    INTEGER(INTG) :: EQUATIONS_SET_IDX,CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER

    CHARACTER(20) :: FILE,OUTPUT_FILE

    ENTERS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%SOLVERS
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%CONTROL_LOOP
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
            IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification array is not allocated.",err,error,*999)
            ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
            END IF
            CALL SYSTEM('mkdir -p ./output')
            SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
              CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
                & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    FILENAME="./output/"//"STATIC_SOLUTION"
                    METHOD="FORTRAN"
                    IF(SOLVER%outputType>=SOLVER_PROGRESS_OUTPUT) THEN
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                      CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                    ENDIF
                    Fields=>EQUATIONS_SET%REGION%FIELDS
                    CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                    CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                    NULLIFY(Fields)
                  ENDDO
                ENDIF
              ENDIF
            CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
              CALL CONTROL_LOOP_TIMES_GET(CONTROL_LOOP,START_TIME,STOP_TIME,CURRENT_TIME,TIME_INCREMENT, &
                & CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER,err,error,*999)
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%SOLVER_MAPPING
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  !Make sure the equations sets are up to date
                  DO equations_set_idx=1,SOLVER_MAPPING%NUMBER_OF_EQUATIONS_SETS
                    EQUATIONS_SET=>SOLVER_MAPPING%EQUATIONS_SETS(equations_set_idx)%ptr
                    IF(OUTPUT_ITERATION_NUMBER/=0) THEN
                      IF(CURRENT_TIME<=STOP_TIME) THEN
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
                        FILENAME="./output/"//"MainTime_"//TRIM(NumberToVString(CURRENT_LOOP_ITERATION,"*",err,error))
                        METHOD="FORTRAN"
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0)  THEN
                          IF(CONTROL_LOOP%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"Now export fields... ",err,error,*999)
                          ENDIF
                          Fields=>EQUATIONS_SET%REGION%FIELDS
                          CALL FIELD_IO_NODES_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                          CALL FIELD_IO_ELEMENTS_EXPORT(Fields,FILENAME,METHOD,err,error,*999)
                          NULLIFY(Fields)
                          IF(CONTROL_LOOP%outputtype >= CONTROL_LOOP_PROGRESS_OUTPUT) THEN
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,FILENAME,err,error,*999)
                            CALL WRITE_STRING(GENERAL_OUTPUT_TYPE,"...",err,error,*999)
                          ENDIF
                        END IF
                        IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
                          CALL AnalyticAnalysis_Output(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FILE,err,error,*999)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDIF
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a BURGERS equation type of a fluid mechanics problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        ELSE
          CALL FlagError("Problem is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Control loop is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver solvers is not associated.",err,error,*999)
    ENDIF
  ELSE
    CALL FlagError("Solver is not associated.",err,error,*999)
  ENDIF

    EXITS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA",err,error)
    RETURN 1

  END SUBROUTINE BURGERS_EQUATION_POST_SOLVE_OUTPUT_DATA

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Burgers problem.
  SUBROUTINE Burgers_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Burgers_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
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
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_BURGERS_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Burgers equation problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

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
  SUBROUTINE BURGERS_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Burgers problem on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)

    ENTERS("BURGERS_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Burgers equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STATIC_BURGERS_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
            CASE(PROBLEM_SETUP_START_ACTION)
              !Set up a simple control loop
              CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
            CASE(PROBLEM_SETUP_FINISH_ACTION)
              !Finish the control loops
              CONTROL_LOOP_ROOT=>PROBLEM%CONTROL_LOOP
              CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
              CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)
            CASE DEFAULT
              localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
                & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a Burgers problem."
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
            !Set the solver to be a static nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_NONLINEAR_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear solver",err,error,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
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
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
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
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Burgers problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(PROBLEM_DYNAMIC_BURGERS_SUBTYPE)
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing????
          CASE DEFAULT
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
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
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
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
            !Set the solver to be a static nonlinear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
            CALL SOLVER_LABEL_SET(SOLVER,"Nonlinear dynamic solver",err,error,*999)
            CALL SOLVER_DYNAMIC_LINEARITY_TYPE_SET(SOLVER,SOLVER_DYNAMIC_NONLINEAR,err,error,*999)
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
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a nonlinear burgers problem."
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
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC, &
              & err,error,*999)
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
            localError="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a Burgers problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a Burgers problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a Burgers problem subtype."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF

    EXITS("BURGERS_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("BURGERS_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1

  END SUBROUTINE BURGERS_EQUATION_PROBLEM_SETUP

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian element stiffness matrices for a BURGERS equation finite element equations set.
  SUBROUTINE Burgers_FiniteElementJacobianEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,nj,ms,nh,nhs,ni,ns,MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: C_PARAM,JGW,SUM1,SUM2,SUM3,DXI_DX,DPHINS_DXI,PHIMS,PHINS,U_VALUE,U_DERIV,VALUE
    LOGICAL :: evaluateJacobian,updateJacobianMatrix
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS1,DEPENDENT_BASIS2,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsJacobianType), POINTER :: jacobianMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME1,QUADRATURE_SCHEME2

    ENTERS("Burgers_FiniteElementJacobianEvaluate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) &
        & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) &
        & CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.",err,error,*999)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        vectorMatrices=>vectorEquations%vectorMatrices
        nonlinearMatrices=>vectorMatrices%nonlinearMatrices
        jacobianMatrix=>nonlinearMatrices%jacobians(1)%ptr
        updateJacobianMatrix=jacobianMatrix%updateJacobian
        evaluateJacobian=updateJacobianMatrix
        IF(evaluateJacobian) THEN
          dependentField=>equations%interpolation%dependentField
          geometricField=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          vectorMapping=>vectorEquations%vectorMapping
          nonlinearMapping=>vectorMapping%nonlinearMapping
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          GEOMETRIC_VARIABLE=>geometricField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) &
            & CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          C_PARAM=1.0_DP
          !Loop over all Gauss points
          DO ng=1,QUADRATURE_SCHEME1%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(SECOND_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              C_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
            ENDIF
            !Loop over rows
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
              DO ms=1,DEPENDENT_BASIS1%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                !Loop over element columns
                DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                  MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(nh)%MESH_COMPONENT_NUMBER
                  DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                    & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                  QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                    &(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                  DO ns=1,DEPENDENT_BASIS2%NUMBER_OF_ELEMENT_PARAMETERS
                    nhs=nhs+1
                    PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                    !Loop over xi directions
                    SUM1=0.0_DP
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      U_DERIV=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr% &
                        & VALUES(mh,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni))
                      DXI_DX=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nh)
                      SUM1=SUM1+U_DERIV*DXI_DX
                    ENDDO !ni
                    SUM1=SUM1*PHINS
                    SUM2=0.0_DP
                    IF(nh==mh) THEN
                      !Loop over spatial directions
                      DO nj=1,GEOMETRIC_VARIABLE%NUMBER_OF_COMPONENTS
                        U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(nj,NO_PART_DERIV)
                        SUM3=0.0_DP
                        DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                          DPHINS_DXI=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          DXI_DX=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                          SUM3=SUM3+DPHINS_DXI*DXI_DX
                        ENDDO !ni
                        SUM2=SUM2+U_VALUE*SUM3
                      ENDDO !nj
                    ENDIF
                    VALUE=C_PARAM*(SUM1+SUM2)*PHIMS
                    jacobianMatrix%elementJacobian%matrix(mhs,nhs)=jacobianMatrix%elementJacobian%matrix(mhs,nhs)+VALUE*JGW
                  ENDDO !ns
                ENDDO !nh
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

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
  SUBROUTINE Burgers_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nj,nh,nhs,ni,mi,ns
    INTEGER(INTG) MESH_COMPONENT1,MESH_COMPONENT2
    REAL(DP) :: A_PARAM,B_PARAM,C_PARAM,JGW,SUM,SUM1,PHIMS,PHINS,DPHIMS_DXI(3),DPHINS_DXI(3),U_VALUE
    LOGICAL :: evaluateAny,evaluateDamping,evaluateLinearDynamic,evaluateResidual,evaluateRHS, &
      & evaluateStiffness,firstDamping,firstRHS,firstStiffness,updateStiffness,updateDamping,updateRHS,updateResidual
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,DEPENDENT_BASIS1,DEPENDENT_BASIS2
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
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: GEOMETRIC_VARIABLE,FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME,QUADRATURE_SCHEME1,QUADRATURE_SCHEME2
    TYPE(VARYING_STRING) :: localError

    ENTERS("Burgers_FiniteElementResidualEvaluate",err,error,*999)

    NULLIFY(dampingMatrix)
    NULLIFY(dynamicMatrices)
    NULLIFY(stiffnessMatrix)

    updateStiffness=.FALSE.
    firstStiffness=.FALSE.
    updateDamping=.FALSE.
    firstDamping=.FALSE.
    updateRHS=.FALSE.
    firstRHS=.FALSE.
    updateResidual=.FALSE.

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) &
        & CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) &
        & CALL FlagError("Equations set specification must have three entries for a Burgers type equations set.",err,error,*999)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        vectorMatrices=>vectorEquations%vectorMatrices
        rhsVector=>vectorMatrices%rhsVector
        nonlinearMatrices=>vectorMatrices%nonlinearMatrices
        vectorMapping=>vectorEquations%vectorMapping
        nonlinearMapping=>vectorMapping%nonlinearMapping
        dependentField=>equations%interpolation%dependentField
        geometricField=>equations%interpolation%geometricField
        GEOMETRIC_VARIABLE=>geometricField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr
        materialsField=>equations%interpolation%materialsField
        GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_BURGERS_SUBTYPE)
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE)
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_STATIC_BURGERS_SUBTYPE)
          linearMatrices=>vectorMatrices%linearMatrices
          stiffnessMatrix=>linearMatrices%matrices(1)%ptr
          linearMapping=>vectorMapping%linearMapping
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE(EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE)
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          dampingMatrix=>dynamicMatrices%matrices(1)%ptr
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>nonlinearMapping%residualVariables(1)%ptr
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a BURGERS equation type of a fluid mechanics equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ASSOCIATED(stiffnessMatrix)) THEN
          updateStiffness=stiffnessMatrix%updateMatrix
          firstStiffness=stiffnessMatrix%firstAssembly
        ENDIF
        IF(ASSOCIATED(dampingMatrix)) THEN
          updateDamping=dampingMatrix%updateMatrix
          firstDamping=dampingMatrix%firstAssembly
        ENDIF
        IF(ASSOCIATED(rhsVector)) THEN
          updateRHS=rhsVector%updateVector
          firstRHS=rhsVector%firstAssembly
        ENDIF
        IF(ASSOCIATED(nonlinearMatrices)) updateResidual=nonlinearMatrices%updateResidual
        evaluateResidual=updateResidual
        evaluateStiffness=firstStiffness.OR.updateStiffness
        evaluateDamping=firstDamping.OR.updateDamping
        evaluateRHS=firstRHS.OR.updateRHS
        evaluateLinearDynamic=evaluateStiffness.OR.evaluateDamping.OR.evaluateRHS
        evaluateAny=evaluateLinearDynamic.OR.updateResidual
        IF(evaluateAny) THEN
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) &
            & CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
          A_PARAM=1.0_DP
          B_PARAM=1.0_DP
          C_PARAM=1.0_DP
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & dependentInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_VAR_TYPE)%ptr,err,error,*999)
            IF(EQUATIONS_SET%SPECIFICATION(3)/=EQUATIONS_SET_INVISCID_BURGERS_SUBTYPE) THEN
              CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_VAR_TYPE)%ptr,err,error,*999)
              IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_GENERALISED_BURGERS_SUBTYPE) THEN
                A_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
                B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV)
                C_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV)
              ELSE
                B_PARAM=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,NO_PART_DERIV)
              ENDIF
            ENDIF
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              MESH_COMPONENT1=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
              DEPENDENT_BASIS1=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT1)%ptr% &
                & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
              QUADRATURE_SCHEME1=>DEPENDENT_BASIS1%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
              JGW=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
                & QUADRATURE_SCHEME1%GAUSS_WEIGHTS(ng)
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                IF(evaluateLinearDynamic) THEN
                  IF(evaluateStiffness.OR.evaluateDamping) THEN
                    nhs=0
                    !Loop over element columns
                    DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      MESH_COMPONENT2=FIELD_VARIABLE%COMPONENTS(mh)%MESH_COMPONENT_NUMBER
                      DEPENDENT_BASIS2=>dependentField%DECOMPOSITION%DOMAIN(MESH_COMPONENT2)%ptr% &
                        & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
                      QUADRATURE_SCHEME2=>DEPENDENT_BASIS2%QUADRATURE%QUADRATURE_SCHEME_MAP&
                        &(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
                      DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                        nhs=nhs+1
                        !Diffusion matrix
                        IF(evaluateStiffness) THEN
                          DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DPHIMS_DXI(ni)=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                            DPHINS_DXI(ni)=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                          END DO !ni
                          SUM=0.0_DP
                          !Calculate SUM
                          DO mi=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                            DO ni=1,DEPENDENT_BASIS2%NUMBER_OF_XI
                              SUM=SUM+DPHINS_DXI(ni)*DPHIMS_DXI(mi)*equations%interpolation% &
                                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%GU(mi,ni)
                            ENDDO !ni
                          ENDDO !mi
                          stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)- &
                            & B_PARAM*SUM*JGW
                        ENDIF
                        !Mass matrix
                        IF(evaluateDamping) THEN
                          PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                          PHINS=QUADRATURE_SCHEME2%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                          dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                            & A_PARAM*PHIMS*PHINS*JGW
                        ENDIF
                      ENDDO !ns
                    ENDDO !nh
                  ENDIF !Stiffness or Damping
                  !Calculate RHS
                  IF(evaluateRHS) THEN
                    rhsVector%elementVector%vector(mhs)=0.0_DP
                  ENDIF
                ENDIF !Evaluate linear dynamic
                !Calculate nonlinear vector
                IF(evaluateResidual) THEN
                  PHIMS=QUADRATURE_SCHEME1%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                  SUM=0.0_DP
                  DO nj=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    U_VALUE=equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(nj,NO_PART_DERIV)
                    SUM1=0.0_DP
                    !Calculate SUM
                    DO ni=1,DEPENDENT_BASIS1%NUMBER_OF_XI
                      SUM1=SUM1+equations%interpolation%dependentInterpPoint(FIELD_VAR_TYPE)%ptr%VALUES(nj, &
                        & PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni))*equations%interpolation% &
                        & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%DXI_DX(ni,nj)
                    ENDDO !ni
                    SUM=SUM+U_VALUE*SUM1
                  ENDDO !nj
                  nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices% &
                    & elementResidual%vector(mhs)+C_PARAM*SUM*PHIMS*JGW
                ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng

          !Scale factor adjustment
          IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,err,error,*999)
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(evaluateStiffness.OR.evaluateDamping) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      IF(evaluateStiffness) &
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                      IF(evaluateDamping) &
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(evaluateRHS) &
                  & rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
                IF(evaluateResidual) &
                  & nonlinearMatrices%elementResidual%vector(mhs)=nonlinearMatrices%elementResidual%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF

        ENDIF !Evaluate any

      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN
999 ERRORS("Burgers_FiniteElementResidualEvaluate",err,error)
    EXITS("Burgers_FiniteElementResidualEvaluate")
    RETURN 1

  END SUBROUTINE Burgers_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

END MODULE BURGERS_EQUATION_ROUTINES
