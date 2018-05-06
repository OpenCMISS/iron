!> \file
!> \author Chris Bradley
!> \brief This module handles all Poiseuille equations routines.
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

!>This module handles all Poiseuille equations routines.
MODULE POISEUILLE_EQUATIONS_ROUTINES

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
! temporary input for vector-source
  USE FLUID_MECHANICS_IO_ROUTINES

#include "macros.h"  

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC POISEUILLE_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC Poiseuille_EquationsSetSolutionMethodSet

  PUBLIC Poiseuille_BoundaryConditionsAnalyticCalculate

  PUBLIC Poiseuille_EquationsSetSpecificationSet

  PUBLIC Poiseuille_FiniteElementCalculate

  PUBLIC POISEUILLE_EQUATION_PROBLEM_SETUP

  PUBLIC POISEUILLE_PRE_SOLVE,POISEUILLE_POST_SOLVE

  PUBLIC Poiseuille_ProblemSpecificationSet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Poiseuille_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,deriv_idx,dim_idx,local_ny,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type
    INTEGER(INTG) :: COUNT_DOF
    REAL(DP) :: VALUE,X(3)
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(VARYING_STRING) :: localError    
    
    ENTERS("Poiseuille_BoundaryConditionsAnalyticCalculate",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          geometricField=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(geometricField)) THEN            
            CALL FIELD_NUMBER_OF_COMPONENTS_GET(geometricField,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
            NULLIFY(GEOMETRIC_VARIABLE)
            CALL Field_VariableGet(geometricField,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
            CALL FIELD_PARAMETER_SET_DATA_GET(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
              & err,error,*999)

            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN

              COUNT_DOF=0

              DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES !U and deludeln
                variable_type=dependentField%VARIABLES(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  CALL FIELD_PARAMETER_SET_CREATE(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS !u,v,w
                    IF(FIELD_VARIABLE%COMPONENTS(component_idx)%interpolation_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
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
                                SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                CASE(EQUATIONS_SET_POISEUILLE_EQUATION_TWO_DIM_1)
                                  !u=ln(4/(x+y+1^2))
                                  SELECT CASE(variable_type)
                                  CASE(FIELD_U_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      VALUE=LOG(4.0_DP/((X(1)+X(2)+1.0_DP)**2))
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE DEFAULT
                                      localError="The global derivative index of "//TRIM(NumberToVString( &
                                      DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & err,error))//" is invalid."
                                      CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                    SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                    CASE(NO_GLOBAL_DERIV)
                                      !This is believed to be incorrect, should be: VALUE=-2.0_DP/(X(1)+X(2)+1.0_DP)
                                      VALUE=-2.0_DP*(X(1)+X(2))/(X(1)+X(2)+1.0_DP)
                                    CASE(GLOBAL_DERIV_S1)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S2)
                                      CALL FlagError("Not implemented.",err,error,*999)
                                    CASE(GLOBAL_DERIV_S1_S2)
                                     CALL FlagError("Not implemented.",err,error,*999)
                                    CASE DEFAULT
                                      localError="The global derivative index of "//TRIM(NumberToVString( &
                                      DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                        & err,error))//" is invalid."
                                      CALL FlagError(localError,err,error,*999)
                                    END SELECT
                                  END SELECT
                                CASE DEFAULT
                                  IF(variable_type==FIELD_U_VARIABLE_TYPE.AND.DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN !For Dirichlet
                                    CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,dependentField,variable_type, &
                                      & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                  ENDIF
                                END SELECT

                                IF(variable_type==FIELD_DELUDELN_VARIABLE_TYPE.and.node_idx/=1) THEN
                                  IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                    !If we are a boundary node then set the analytic value on the boundary
                                    !CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,variable_type,local_ny, &
                                    !  & BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
                                    !Do nothing at present
                                  ENDIF
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
                  CALL FIELD_PARAMETER_SET_UPDATE_START(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                  CALL FIELD_PARAMETER_SET_UPDATE_FINISH(dependentField,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                    & err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO
              CALL FIELD_PARAMETER_SET_DATA_RESTORE(geometricField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                & GEOMETRIC_PARAMETERS,err,error,*999)
            ELSE
              CALL FlagError("Boundary conditions is not associated.",err,error,*999)
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
    
    EXITS("Poiseuille_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("Poiseuille_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("Poiseuille_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Poiseuille_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Sets up the Poiseuille equation type of a fluid mechanics equations set class.
  SUBROUTINE POISEUILLE_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Poiseuille equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("POISEUILLE_EQUATION_EQUATIONS_SET_SETUP",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Poiseuille type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE)
        CALL Poiseuille_EquationsSetStaticSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_DYNAMIC_POISEUILLE_SUBTYPE)
        CALL Poiseuille_EquationsSetStaticSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE DEFAULT
        localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Poiseuille equation type of a fluid mechanics equation set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("POISEUILLE_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("POISEUILLE_EQUATION_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE POISEUILLE_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Poiseuille equation type of an fluid mechanics equations set class.
  SUBROUTINE Poiseuille_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Poiseuille_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Poiseuille type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE,EQUATIONS_SET_DYNAMIC_POISEUILLE_SUBTYPE)        
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
          & " is not valid for a Poiseuille equation type of an fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Poiseuille_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("Poiseuille_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE Poiseuille_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Poiseuille fluid mechanics equations set.
  SUBROUTINE Poiseuille_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Poiseuille_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Poiseuille equation type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE, &
          & EQUATIONS_SET_DYNAMIC_POISEUILLE_SUBTYPE)
        !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a Poiseuille equation type of a fluid mechanics equations set class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_POISEUILLE_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Poiseuille_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Poiseuille_EquationsSetSpecificationSet",err,error)
    EXITS("Poiseuille_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Poiseuille_EquationsSetSpecificationSet


  !
  !================================================================================================================================
  !

  !>Sets up the standard Poiseuille equation for linear sources.
  SUBROUTINE Poiseuille_EquationsSetStaticSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("Poiseuille_EquationsSetStaticSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Poiseuille type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Poiseuille_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,ERR, &
              & ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!TODO: Check valid setup
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing???
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
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
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2, &
                & err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2, &
                & GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                & GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2, &
                  & FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,2,FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
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
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2, &
                  & FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                 CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                  & FIELD_ELEMENT_BASED_INTERPOLATION,err,error,*999)
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
              & " is invalid for a static Poiseuille equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Create the auto created materials field
                CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & MATERIALS_FIELD,err,error,*999)
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
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                !Set the number of materials components
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 3,err,error,*999)
                !Default the 3 materials components to the first geometric interpolation setup with constant interpolation
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                DO component_idx=1,3
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
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
                CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3,err,error,*999)
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
                !Default component 1 (viscosity) 
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,0.005_DP,err,error,*999)
                !Default component 2 (radius) 
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.5_DP,err,error,*999)
                !Default component 3 (length) 
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,3,1.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE)
              SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
               !Set start action
                CASE(EQUATIONS_SET_SETUP_START_ACTION)
                  !Add in gravity source field
                !Specify finish action
                CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
                  !Add in gravity source field
                CASE DEFAULT
                  localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
                  & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
                  & " is invalid for a linear Poiseuille subtype"
                  CALL FlagError(localError,err,error,*999)
              END SELECT
            CASE DEFAULT
              localError="The equation set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
                & " is invalid for a linear Poiseuille equation."
               CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)

          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            
          CASE(EQUATIONS_SET_SETUP_GENERATE_ACTION)
             
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              !Create the equations
              CALL Equations_CreateStart(EQUATIONS_SET,equations,err,error,*999)
              CALL Equations_LinearityTypeSet(equations,EQUATIONS_LINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(equations,EQUATIONS_STATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the creation of the equations
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,equations,err,error,*999)
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
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a static Poiseuille equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not a static Poiseuille equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Poiseuille_EquationsSetStaticSetup")
    RETURN
999 ERRORSEXITS("Poiseuille_EquationsSetStaticSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Poiseuille_EquationsSetStaticSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the Poiseuille problem.
  SUBROUTINE POISEUILLE_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Poiseuille equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("POISEUILLE_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Poiseuille problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STATIC_POISEUILLE_SUBTYPE)
        CALL POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE)
        CALL POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE DEFAULT
        localError="Problem subtype "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Poiseuille equation type of a fluid mechanics problem class."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("POISEUILLE_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("POISEUILLE_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE POISEUILLE_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Poiseuille equation finite element equations set.
  SUBROUTINE Poiseuille_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField,materialsField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: localError

    ENTERS("Poiseuille_FiniteElementCalculate",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) &
        & CALL FlagError("Equations set specification must have three entries for a Poiseuille type equations set.",err,error,*999)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STATIC_POISEUILLE_SUBTYPE,EQUATIONS_SET_DYNAMIC_POISEUILLE_SUBTYPE)
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
          FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%variable
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(dependentField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)

!TODO

            
          ENDDO !ng
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a Poiseuille equation type of a fluid mechanics equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Poiseuille_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Poiseuille_FiniteElementCalculate",err,error)
    RETURN 1
    
  END SUBROUTINE Poiseuille_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Poiseuille problem.
  SUBROUTINE Poiseuille_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Poiseuille_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STATIC_POISEUILLE_SUBTYPE, &
            & PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE)
          !ALl ok
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Poiseuille fluid mechanics problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_FLUID_MECHANICS_CLASS,PROBLEM_POISEUILLE_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Poiseuille problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("Poiseuille_ProblemSpecificationSet")
    RETURN
999 ERRORS("Poiseuille_ProblemSpecificationSet",err,error)
    EXITS("Poiseuille_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Poiseuille_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the static Poiseuille equations problem.
  SUBROUTINE POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: localError
   
    ENTERS("POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(problem%specification)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(problem%specification,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Poiseuille problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_STATIC_POISEUILLE_SUBTYPE) THEN
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
              & " is invalid for a static Poiseuille equation."
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
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
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
            !Set the solver to be a linear solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
             !Start the linear solver creation
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
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
            localError="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a static Poiseuille equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a static Poiseuille equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        localError="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a static Poiseuille equation subtype."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP")
    RETURN
999 ERRORSEXITS("POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP",err,error)
    RETURN 1
  END SUBROUTINE POISEUILLE_EQUATION_PROBLEM_STATIC_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets up the Poiseuille problem post solve.
  SUBROUTINE POISEUILLE_POST_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solver
    TYPE(VARYING_STRING) :: localError

    ENTERS("POISEUILLE_POST_SOLVE",err,error,*999)
    NULLIFY(SOLVER2)
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Poiseuille problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_POISEUILLE_SUBTYPE)
!               do nothing
            CASE(PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE)
!               do nothing
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Poiseuille type of a fluid mechanics problem class."
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

    EXITS("POISEUILLE_POST_SOLVE")
    RETURN
999 ERRORSEXITS("POISEUILLE_POST_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE POISEUILLE_POST_SOLVE

  !
  !================================================================================================================================
  !


  !>Sets up the Poiseuille problem pre solve.
  SUBROUTINE POISEUILLE_PRE_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(SOLVER_TYPE), POINTER :: SOLVER2 !<A pointer to the solvers
    TYPE(VARYING_STRING) :: localError

    ENTERS("POISEUILLE_PRE_SOLVE",err,error,*999)
    NULLIFY(SOLVER2)
 
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
          IF(.NOT.ALLOCATED(CONTROL_LOOP%problem%specification)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%problem%specification,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Poiseuille problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_STATIC_POISEUILLE_SUBTYPE)
!               do nothing
            CASE(PROBLEM_DYNAMIC_POISEUILLE_SUBTYPE)
!               do nothing
            CASE DEFAULT
              localError="Problem subtype "//TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is not valid for a Poiseuille type of a fluid mechanics problem class."
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

    EXITS("POISEUILLE_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("POISEUILLE_PRE_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE POISEUILLE_PRE_SOLVE

  !
  !================================================================================================================================
  !
 

END MODULE POISEUILLE_EQUATIONS_ROUTINES
