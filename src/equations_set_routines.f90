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
MODULE EQUATIONS_SET_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BIOELECTRIC_ROUTINES
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE CLASSICAL_FIELD_ROUTINES
  USE CmissMPI
  USE ComputationEnvironment
  USE Constants
  USE COORDINATE_ROUTINES
  USE DistributedMatrixVector
  USE DOMAIN_MAPPINGS
  USE ELASTICITY_ROUTINES
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingAccessRoutines
  USE EquationsMatricesRoutines
  USE EquationsMatricesAccessRoutines
  USE EquationsSetAccessRoutines
  USE EquationsSetConstants
  USE FIELD_ROUTINES
  USE FieldAccessRoutines
  USE FittingRoutines
  USE FLUID_MECHANICS_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Lists
  USE MatrixVector
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

  PUBLIC EQUATIONS_SET_ANALYTIC_CREATE_START,EQUATIONS_SET_ANALYTIC_CREATE_FINISH

  PUBLIC EQUATIONS_SET_ANALYTIC_DESTROY

  PUBLIC EQUATIONS_SET_ANALYTIC_EVALUATE

  PUBLIC EQUATIONS_SET_ASSEMBLE
  
  PUBLIC EQUATIONS_SET_BACKSUBSTITUTE,EQUATIONS_SET_NONLINEAR_RHS_UPDATE
  
  PUBLIC EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC

  PUBLIC EQUATIONS_SET_CREATE_START,EQUATIONS_SET_CREATE_FINISH

  PUBLIC EQUATIONS_SET_DESTROY

  PUBLIC EQUATIONS_SETS_FINALISE,EQUATIONS_SETS_INITIALISE

  PUBLIC EQUATIONS_SET_EQUATIONS_CREATE_FINISH,EQUATIONS_SET_EQUATIONS_CREATE_START

  PUBLIC EQUATIONS_SET_EQUATIONS_DESTROY
  
  PUBLIC EQUATIONS_SET_MATERIALS_CREATE_START,EQUATIONS_SET_MATERIALS_CREATE_FINISH

  PUBLIC EQUATIONS_SET_MATERIALS_DESTROY
  
  PUBLIC EQUATIONS_SET_DEPENDENT_CREATE_START,EQUATIONS_SET_DEPENDENT_CREATE_FINISH

  PUBLIC EQUATIONS_SET_DEPENDENT_DESTROY

  PUBLIC EquationsSet_DerivedCreateStart,EquationsSet_DerivedCreateFinish

  PUBLIC EquationsSet_DerivedDestroy
  
  PUBLIC EQUATIONS_SET_INDEPENDENT_CREATE_START,EQUATIONS_SET_INDEPENDENT_CREATE_FINISH

  PUBLIC EQUATIONS_SET_INDEPENDENT_DESTROY
  
  PUBLIC EquationsSet_JacobianEvaluate,EquationsSet_ResidualEvaluate

  PUBLIC EquationsSet_OutputTypeGet,EquationsSet_OutputTypeSet
  
  PUBLIC EQUATIONS_SET_SOLUTION_METHOD_GET,EQUATIONS_SET_SOLUTION_METHOD_SET
  
  PUBLIC EQUATIONS_SET_SOURCE_CREATE_START,EQUATIONS_SET_SOURCE_CREATE_FINISH

  PUBLIC EQUATIONS_SET_SOURCE_DESTROY

  PUBLIC EquationsSet_SpecificationGet,EquationsSet_SpecificationSizeGet

  PUBLIC EquationsSet_TensorInterpolateGaussPoint

  PUBLIC EquationsSet_TensorInterpolateXi

  PUBLIC EquationsSet_DerivedVariableCalculate,EquationsSet_DerivedVariableSet

  PUBLIC EQUATIONS_SET_LOAD_INCREMENT_APPLY
  
  PUBLIC EQUATIONS_SET_ANALYTIC_USER_PARAM_SET,EQUATIONS_SET_ANALYTIC_USER_PARAM_GET

  PUBLIC EquationsSet_TimesGet,EquationsSet_TimesSet

CONTAINS

  !
  !================================================================================================================================
  !
      
  !>Finish the creation of a analytic solution for equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticCreateFinish
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_FINISH(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create the analytic for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD

    ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
          CALL FlagError("Equations set analytic has already been finished.",err,error,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
          IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=ANALYTIC_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>ANALYTIC_FIELD
          ENDIF
          !Finish the equations set specific analytic setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finish the analytic creation
          EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FlagError("The equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a analytic solution for a equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticCreateStart
  SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_START(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,ANALYTIC_FIELD_USER_NUMBER,ANALYTIC_FIELD, &
    & err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of an analytic for.
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The analytic function type to setup \see EquationsSetConstants_AnalyticFunctionTypes,EquationsSetConstants
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FIELD_USER_NUMBER !<The user specified analytic field number
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD !<If associated on entry, a pointer to the user created analytic field which has the same user number as the specified analytic field user number. If not associated on entry, on exit, a pointer to the created analytic field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,ANALYTIC_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EQUATIONS_SET_ANALYTIC_CREATE_START",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        CALL FlagError("The equations set analytic is already associated.",err,error,*998)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
            !Check the analytic field has been finished
            IF(ANALYTIC_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(ANALYTIC_FIELD_USER_NUMBER/=ANALYTIC_FIELD%USER_NUMBER) THEN
                localError="The specified analytic field user number of "// &
                  & TRIM(NumberToVString(ANALYTIC_FIELD_USER_NUMBER,"*",err,error))// &
                  & " does not match the user number of the specified analytic field of "// &
                  & TRIM(NumberToVString(ANALYTIC_FIELD%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              ANALYTIC_FIELD_REGION=>ANALYTIC_FIELD%REGION
              IF(ASSOCIATED(ANALYTIC_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(ANALYTIC_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified analytic field has been created on region number "// &
                    & TRIM(NumberToVString(ANALYTIC_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Check the specified analytic field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,ANALYTIC_FIELD%DECOMPOSITION)) THEN
                    CALL FlagError("The specified analytic field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified analytic field region is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The specified analytic field has not been finished.",err,error,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(ANALYTIC_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
            IF(ASSOCIATED(FIELD)) THEN
              localError="The specified analytic field user number of "// &
                & TRIM(NumberToVString(ANALYTIC_FIELD_USER_NUMBER,"*",err,error))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
          !Initialise the equations set analytic
          CALL EQUATIONS_SET_ANALYTIC_INITIALISE(EQUATIONS_SET,err,error,*999)
          IF(.NOT.ASSOCIATED(ANALYTIC_FIELD)) EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_ANALYTIC_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=ANALYTIC_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>ANALYTIC_FIELD
          EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=ANALYTIC_FUNCTION_TYPE
          !Start the equations set specific analytic setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Set pointers
          IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
            ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
          ELSE
            EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD=>ANALYTIC_FIELD
          ENDIF
        ELSE
          CALL FlagError("Equations set region is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroy the analytic solution for an equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticDestroy
  SUBROUTINE EQUATIONS_SET_ANALYTIC_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the analytic solutins for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_ANALYTIC_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN        
        CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,err,error,*999)
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the current analytic solution for an equations set. \see OpenCMISS::cmfe_EquationsSet_AnalyticEvaluate
  SUBROUTINE EQUATIONS_SET_ANALYTIC_EVALUATE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the current analytic solutins for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,derivative_idx,element_idx,Gauss_idx,GLOBAL_DERIV_INDEX,local_ny,node_idx, &
      & NUMBER_OF_ANALYTIC_COMPONENTS,NUMBER_OF_DIMENSIONS,variable_idx, &
      & variable_type,version_idx
    REAL(DP) :: NORMAL(3),POSITION(3),TANGENTS(3,3),VALUE
    REAL(DP) :: ANALYTIC_DUMMY_VALUES(1)=0.0_DP
    REAL(DP) :: MATERIALS_DUMMY_VALUES(1)=0.0_DP
    LOGICAL :: reverseNormal=.FALSE.
    TYPE(BASIS_TYPE), POINTER :: BASIS
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: DOMAIN_ELEMENTS
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,dependentField,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_INTERPOLATION_PARAMETERS_PTR_TYPE), POINTER :: ANALYTIC_INTERP_PARAMETERS(:),GEOMETRIC_INTERP_PARAMETERS(:), &
      & MATERIALS_INTERP_PARAMETERS(:)
    TYPE(FIELD_INTERPOLATED_POINT_PTR_TYPE), POINTER :: ANALYTIC_INTERP_POINT(:),GEOMETRIC_INTERP_POINT(:), &
      & MATERIALS_INTERP_POINT(:)
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_PTR_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS(:)
    TYPE(FIELD_PHYSICAL_POINT_PTR_TYPE), POINTER :: ANALYTIC_PHYSICAL_POINT(:),MATERIALS_PHYSICAL_POINT(:)
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_ANALYTIC_EVALUATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED) THEN
          dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          IF(ASSOCIATED(dependentField)) THEN
            GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
            IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
              CALL Field_NumberOfComponentsGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
              CALL Field_InterpolationParametersInitialise(GEOMETRIC_FIELD,GEOMETRIC_INTERP_PARAMETERS,err,error,*999)
              CALL Field_InterpolatedPointsInitialise(GEOMETRIC_INTERP_PARAMETERS,GEOMETRIC_INTERP_POINT,err,error,*999)
              CALL Field_InterpolatedPointsMetricsInitialise(GEOMETRIC_INTERP_POINT,GEOMETRIC_INTERPOLATED_POINT_METRICS, &
                & err,error,*999)
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                CALL Field_NumberOfComponentsGet(ANALYTIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
                  & err,error,*999)
                CALL Field_InterpolationParametersInitialise(ANALYTIC_FIELD,ANALYTIC_INTERP_PARAMETERS,err,error,*999)
                CALL Field_InterpolatedPointsInitialise(ANALYTIC_INTERP_PARAMETERS,ANALYTIC_INTERP_POINT,err,error,*999)
                CALL Field_PhysicalPointsInitialise(ANALYTIC_INTERP_POINT,GEOMETRIC_INTERP_POINT,ANALYTIC_PHYSICAL_POINT, &
                  & err,error,*999)
              ENDIF
              NULLIFY(MATERIALS_FIELD)
              IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
                MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
                CALL Field_NumberOfComponentsGet(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_ANALYTIC_COMPONENTS, &
                  & err,error,*999)
                CALL Field_InterpolationParametersInitialise(MATERIALS_FIELD,MATERIALS_INTERP_PARAMETERS,err,error,*999)
                CALL Field_InterpolatedPointsInitialise(MATERIALS_INTERP_PARAMETERS,MATERIALS_INTERP_POINT,err,error,*999)
                CALL Field_PhysicalPointsInitialise(MATERIALS_INTERP_POINT,GEOMETRIC_INTERP_POINT,MATERIALS_PHYSICAL_POINT, &
                  & err,error,*999)
              ENDIF
              DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
                variable_type=dependentField%variables(variable_idx)%VARIABLE_TYPE
                FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(variable_type)%ptr
                IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                  DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                    IF(ASSOCIATED(DOMAIN)) THEN
                      IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                        SELECT CASE(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE)
                        CASE(FIELD_CONSTANT_INTERPOLATION)
                          CALL FlagError("Cannot evaluate an analytic solution for a constant interpolation components.", &
                            & err,error,*999)
                        CASE(FIELD_ELEMENT_BASED_INTERPOLATION)
                          DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
                          IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                            !Loop over the local elements excluding the ghosts
                            DO element_idx=1,DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                              BASIS=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)%BASIS
                              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                  & ANALYTIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                  & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              CALL FIELD_INTERPOLATE_XI(FIRST_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
                                & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE, &
                                & GEOMETRIC_INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              CALL Field_PositionNormalTangentsCalculateIntPtMetric( &
                                & GEOMETRIC_INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%ptr,reverseNormal, &
                                & POSITION,NORMAL,TANGENTS,err,error,*999)
                              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
                                  & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                CALL FIELD_INTERPOLATE_XI(NO_PART_DERIV,[0.5_DP,0.5_DP,0.5_DP], &
                                  & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
!! \todo Maybe do this with optional arguments?
                              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                  CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                    & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                    & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                    & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                    & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                    & VALUE,err,error,*999)
                                ELSE
                                  CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                    & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                    & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                    & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                    & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                ENDIF
                              ELSE
                                IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                  CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                    & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                    & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                    & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                    & VALUE,err,error,*999)
                                ELSE
                                  CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                    & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                    & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                    & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                ENDIF
                              ENDIF
                              local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                & ELEMENT_PARAM2DOF_MAP%ELEMENTS(element_idx)
                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                            ENDDO !element_idx
                          ELSE
                            CALL FlagError("Domain topology elements is not associated.",err,error,*999)
                          ENDIF
                        CASE(FIELD_NODE_BASED_INTERPOLATION)
                          DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                          IF(ASSOCIATED(DOMAIN_NODES)) THEN
                            !Loop over the local nodes excluding the ghosts.
                            DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                              CALL Field_PositionNormalTangentsCalculateNode(dependentField,variable_type,component_idx, &
                                & node_idx,POSITION,NORMAL,TANGENTS,err,error,*999)
                              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                CALL Field_InterpolateFieldNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE,ANALYTIC_FIELD, &
                                  & FIELD_U_VARIABLE_TYPE,component_idx,node_idx,ANALYTIC_PHYSICAL_POINT( &
                                  & FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                CALL Field_InterpolateFieldNode(NO_PHYSICAL_DERIV,FIELD_VALUES_SET_TYPE,MATERIALS_FIELD, &
                                  & FIELD_U_VARIABLE_TYPE,component_idx,node_idx,MATERIALS_PHYSICAL_POINT( &
                                  & FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              !Loop over the derivatives
                              DO derivative_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES                                
                                GLOBAL_DERIV_INDEX=DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)% &
                                  & GLOBAL_DERIVATIVE_INDEX
!! \todo Maybe do this with optional arguments?
                                IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                  IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                      & ANALYTIC_PHYSICAL_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES, &
                                      & MATERIALS_PHYSICAL_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES,VALUE,err,error,*999)
                                  ELSE
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                      & ANALYTIC_PHYSICAL_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES, &
                                      & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                  ENDIF
                                ELSE
                                  IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                      & MATERIALS_PHYSICAL_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES,VALUE,err,error,*999)
                                  ELSE
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                      & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                  ENDIF
                                ENDIF
                                !Loop over the versions
                                DO version_idx=1,DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(derivative_idx)%numberOfVersions
                                  local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                    & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(derivative_idx)%VERSIONS(version_idx)
                                  CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                    & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                ENDDO !version_idx
                              ENDDO !deriv_idx
                            ENDDO !node_idx
                          ELSE
                            CALL FlagError("Domain topology nodes is not associated.",err,error,*999)
                          ENDIF
                        CASE(FIELD_GRID_POINT_BASED_INTERPOLATION)
                          CALL FlagError("Not implemented.",err,error,*999)
                        CASE(FIELD_GAUSS_POINT_BASED_INTERPOLATION)
                          DOMAIN_ELEMENTS=>DOMAIN%TOPOLOGY%ELEMENTS
                          IF(ASSOCIATED(DOMAIN_ELEMENTS)) THEN
                            !Loop over the local elements excluding the ghosts
                            DO element_idx=1,DOMAIN_ELEMENTS%NUMBER_OF_ELEMENTS
                              BASIS=>DOMAIN_ELEMENTS%ELEMENTS(element_idx)%BASIS
                              CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                & GEOMETRIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                  & ANALYTIC_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,element_idx, &
                                  & MATERIALS_INTERP_PARAMETERS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                              ENDIF
                              !Loop over the Gauss points in the element
                              DO gauss_idx=1,BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr% &
                                & NUMBER_OF_GAUSS
                                CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                                  & GEOMETRIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                CALL Field_InterpolatedPointMetricsCalculate(COORDINATE_JACOBIAN_NO_TYPE, &
                                  & GEOMETRIC_INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                CALL Field_PositionNormalTangentsCalculateIntPtMetric( &
                                  & GEOMETRIC_INTERPOLATED_POINT_METRICS(FIELD_U_VARIABLE_TYPE)%ptr,reverseNormal, &
                                  & POSITION,NORMAL,TANGENTS,err,error,*999)
                                IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                                    & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ENDIF
                                IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                  CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,gauss_idx, &
                                    & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                                ENDIF
!! \todo Maybe do this with optional arguments?
                                IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                                  IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                      & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                      & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                      & VALUE,err,error,*999)
                                  ELSE
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx, &
                                      & ANALYTIC_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                      & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                  ENDIF
                                ELSE
                                  IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                      & MATERIALS_INTERP_POINT(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(:,NO_PART_DERIV), &
                                      & VALUE,err,error,*999)
                                  ELSE
                                    CALL EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%ANALYTIC% &
                                      & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME, &
                                      & variable_type,GLOBAL_DERIV_INDEX,component_idx,ANALYTIC_DUMMY_VALUES, &
                                      & MATERIALS_DUMMY_VALUES,VALUE,err,error,*999)
                                  ENDIF
                                ENDIF
                                local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                  & GAUSS_POINT_PARAM2DOF_MAP%GAUSS_POINTS(Gauss_idx,element_idx)
                                CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                  & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                              ENDDO !Gauss_idx
                            ENDDO !element_idx
                          ELSE
                            CALL FlagError("Domain topology elements is not associated.",err,error,*999)
                          ENDIF
                        CASE DEFAULT
                          localError="The interpolation type of "//TRIM(NumberToVString(FIELD_VARIABLE% &
                            & COMPONENTS(component_idx)%INTERPOLATION_TYPE,"*",err,error))// &
                            & " for component "//TRIM(NumberToVString(component_idx,"*",err,error))//" of variable type "// &
                            & TRIM(NumberToVString(variable_type,"*",err,error))//" is invalid."
                          CALL FlagError(localError,err,error,*999)
                        END SELECT
                      ELSE
                        CALL FlagError("Domain topology is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Domain is not associated.",err,error,*999)
                    ENDIF
                  ENDDO !component_idx
                  CALL Field_ParameterSetUpdateStart(dependentField,variable_type, &
                    & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                  CALL Field_ParameterSetUpdateFinish(dependentField,variable_type, &
                    & FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                ELSE
                  CALL FlagError("Field variable is not associated.",err,error,*999)
                ENDIF
              ENDDO !variable_idx
              IF(ASSOCIATED(MATERIALS_FIELD)) THEN
                CALL Field_PhysicalPointsFinalise(MATERIALS_PHYSICAL_POINT,err,error,*999)
                CALL Field_InterpolatedPointsFinalise(MATERIALS_INTERP_POINT,err,error,*999)
                CALL Field_InterpolationParametersFinalise(MATERIALS_INTERP_PARAMETERS,err,error,*999)
              ENDIF
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                CALL Field_PhysicalPointsFinalise(ANALYTIC_PHYSICAL_POINT,err,error,*999)
                CALL Field_InterpolatedPointsFinalise(ANALYTIC_INTERP_POINT,err,error,*999)
                CALL Field_InterpolationParametersFinalise(ANALYTIC_INTERP_PARAMETERS,err,error,*999)
              ENDIF
              CALL Field_InterpolatedPointsMetricsFinalise(GEOMETRIC_INTERPOLATED_POINT_METRICS,err,error,*999)
              CALL Field_InterpolatedPointsFinalise(GEOMETRIC_INTERP_POINT,err,error,*999)
              CALL Field_InterpolationParametersFinalise(GEOMETRIC_INTERP_PARAMETERS,err,error,*999)
              
            ELSE
              CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set analytic has not been finished.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_EVALUATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_EVALUATE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_EVALUATE

  !
  !================================================================================================================================
  !

  !>Finalise the analytic solution for an equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET_ANALYTIC,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: EQUATIONS_SET_ANALYTIC !<A pointer to the equations set analytic to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_ANALYTIC_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET_ANALYTIC)) THEN        
      DEALLOCATE(EQUATIONS_SET_ANALYTIC)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_FINALISE

  !
  !================================================================================================================================
  !

  !>Evaluate the analytic solution for an equations set.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,TIME, &
    & VARIABLE_TYPE,GLOBAL_DERIVATIVE,COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to evaluate the analytic for
    INTEGER(INTG), INTENT(IN) :: ANALYTIC_FUNCTION_TYPE !<The type of analytic function to evaluate
    REAL(DP), INTENT(IN) :: POSITION(:) !<POSITION(dimention_idx). The geometric position to evaluate at
    REAL(DP), INTENT(IN) :: TANGENTS(:,:) !<TANGENTS(dimention_idx,xi_idx). The geometric tangents at the point to evaluate at.
    REAL(DP), INTENT(IN) :: NORMAL(:) !<NORMAL(dimension_idx). The normal vector at the point to evaluate at.
    REAL(DP), INTENT(IN) :: TIME !<The time to evaluate at
    INTEGER(INTG), INTENT(IN) :: VARIABLE_TYPE !<The field variable type to evaluate at
    INTEGER(INTG), INTENT(IN) :: GLOBAL_DERIVATIVE !<The global derivative direction to evaluate at
    INTEGER(INTG), INTENT(IN) :: COMPONENT_NUMBER !<The dependent field component number to evaluate
    REAL(DP), INTENT(IN) :: ANALYTIC_PARAMETERS(:) !<A pointer to any analytic field parameters
    REAL(DP), INTENT(IN) :: MATERIALS_PARAMETERS(:) !<A pointer to any materials field parameters
    REAL(DP), INTENT(OUT) :: VALUE !<On return, the analtyic function value.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<1) THEN
        CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(1))
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
          CALL FlagError("Equations set specification must have at least two entries for a "// &
            & "classical field equations set.",err,error,*999)
        END IF
        CALL CLASSICAL_FIELD_ANALYTIC_FUNCTIONS_EVALUATE(EQUATIONS_SET,EQUATIONS_SET%SPECIFICATION(2), &
          & ANALYTIC_FUNCTION_TYPE,POSITION,TANGENTS,NORMAL,TIME,VARIABLE_TYPE,GLOBAL_DERIVATIVE, &
          & COMPONENT_NUMBER,ANALYTIC_PARAMETERS,MATERIALS_PARAMETERS,VALUE,err,error,*999)
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
          & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(1),"*",err,error))//" is not valid."
        CALL FLAG_ERROR(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_FUNCTIONS_EVALUATE

  !
  !================================================================================================================================
  !

  !>Initialises the analytic solution for an equations set.
  SUBROUTINE EQUATIONS_SET_ANALYTIC_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("EQUATIONS_SET_ANALYTIC_INITIALISE",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        CALL FlagError("Analytic is already associated for this equations set.",err,error,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%ANALYTIC,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate equations set analytic.",err,error,*999)
        EQUATIONS_SET%ANALYTIC%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%ANALYTIC%ANALYTIC_FINISHED=.FALSE.
        EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD)
        EQUATIONS_SET%ANALYTIC%ANALYTIC_TIME=0.0_DP
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_ANALYTIC_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_INITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_INITIALISE

  !
  !================================================================================================================================
  !

  !>Sets the analytic problem user parameter
  SUBROUTINE EQUATIONS_SET_ANALYTIC_USER_PARAM_SET(EQUATIONS_SET,PARAM_IDX,PARAM,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(IN) :: PARAM_IDX !<Index of the user parameter
    REAL(DP), INTENT(IN) :: PARAM !<Value of the parameter
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: ANALYTIC

    ENTERS("EQUATIONS_SET_ANALYTIC_USER_PARAM_SET",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      ANALYTIC=>EQUATIONS_SET%ANALYTIC
      IF(ASSOCIATED(ANALYTIC)) THEN
        IF(PARAM_IDX>0.AND.PARAM_IDX<=SIZE(ANALYTIC%ANALYTIC_USER_PARAMS)) THEN
          !Set the value
          ANALYTIC%ANALYTIC_USER_PARAMS(PARAM_IDX)=PARAM
        ELSE
          CALL FlagError("Invalid parameter index.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF    
    ELSE 
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SET_ANALYTIC_USER_PARAM_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_USER_PARAM_SET",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_USER_PARAM_SET

  !
  !================================================================================================================================
  !

  !>Sets the analytic problem user parameter
  SUBROUTINE EQUATIONS_SET_ANALYTIC_USER_PARAM_GET(EQUATIONS_SET,PARAM_IDX,PARAM,err,error,*)
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the analytic solution for.
    INTEGER(INTG), INTENT(IN) :: PARAM_IDX !<Index of the user parameter
    REAL(DP), INTENT(OUT) :: PARAM !<Value of the parameter
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local variables
    TYPE(EQUATIONS_SET_ANALYTIC_TYPE), POINTER :: ANALYTIC

    ENTERS("EQUATIONS_SET_ANALYTIC_USER_PARAM_GET",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      ANALYTIC=>EQUATIONS_SET%ANALYTIC
      IF(ASSOCIATED(ANALYTIC)) THEN
        IF(PARAM_IDX>0.AND.PARAM_IDX<=SIZE(ANALYTIC%ANALYTIC_USER_PARAMS)) THEN
          !Set the value
          PARAM=ANALYTIC%ANALYTIC_USER_PARAMS(PARAM_IDX)
        ELSE
          CALL FlagError("Invalid parameter index.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set analytic is not associated.",err,error,*999)
      ENDIF    
    ELSE 
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SET_ANALYTIC_USER_PARAM_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ANALYTIC_USER_PARAM_GET",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_ANALYTIC_USER_PARAM_GET

  !
  !================================================================================================================================
  !

  !>Assembles the equations for an equations set.
  SUBROUTINE EQUATIONS_SET_ASSEMBLE(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EquationsType), POINTER :: equations
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EQUATIONS_SET_ASSEMBLE",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    NULLIFY(equations)
    CALL EquationsSet_EquationsGet(equationsSet,equations,err,error,*999)
   
    IF(.NOT.equations%equationsFinished) CALL FlagError("Equations have not been finished.",err,error,*999)
    
    IF(equationsSet%outputType>=EQUATIONS_SET_PROGRESS_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Equations set assemble: ",equationsSet%label,err,error,*999)
    ENDIF
    
    SELECT CASE(equations%timeDependence)
    CASE(EQUATIONS_STATIC)
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_NONLINEAR)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
      ! chrm, 17/09/09
      SELECT CASE(equations%linearity)
      CASE(EQUATIONS_LINEAR)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
    
    EXITS("EQUATIONS_SET_ASSEMBLE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_ASSEMBLE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_ASSEMBLE

  !
  !================================================================================================================================
  !
  
  !>Assembles the equations stiffness matrix and rhs for a dynamic linear equations set using the finite element method.
  SUBROUTINE EquationsSet_AssembleDynamicLinearFEM(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
      & mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    
!#ifdef TAUPROF
!    CHARACTER(28) :: CVAR
!    INTEGER :: PHASE(2)= [ 0, 0 ]
!    SAVE PHASE
!#endif

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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_VectorValuesInitialise()")
#endif
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_VectorValuesInitialise()")
#endif
    !Assemble the elements
    !Allocate the element matrices 
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_ElementInitialise()")
#endif
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
      & mappings%elements
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_ElementInitialise()")
#endif
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
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("Internal Elements Loop")
#endif
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      !#ifdef TAUPROF
      !              WRITE (CVAR,'(a23,i3)') 'Internal Elements Loop ',elementIdx
      !              CALL TAU_PHASE_CREATE_DYNAMIC(PHASE,CVAR)
      !              CALL TAU_PHASE_START(PHASE)
      !#endif
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
      !#ifdef TAUPROF
      !              CALL TAU_PHASE_STOP(PHASE)
      !#endif
    ENDDO !elementIdx
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("Internal Elements Loop")
#endif
    
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
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("Boundary and Ghost Elements Loop")
#endif
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
    ENDDO !elementIdx
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("Boundary and Ghost Elements Loop")
#endif
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
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsMatrices_ElementFinalise()")
#endif
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsMatrices_ElementFinalise()")
#endif
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
      & mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_AssembleQuasistaticNonlinearFEM",err,error,*999)

    ! currently no difference
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_LINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
      & mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementCalculate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
  SUBROUTINE EQUATIONS_SET_BACKSUBSTITUTE(equationsSet,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to backsubstitute
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<The boundary conditions to use for the backsubstitution
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_column_idx,equations_column_number,equations_matrix_idx,equations_row_number, &
      & EQUATIONS_STORAGE_TYPE,rhs_boundary_condition,rhs_global_dof,rhs_variable_dof,rhsVariableType,variable_dof,VARIABLE_TYPE
    INTEGER(INTG), POINTER :: COLUMN_INDICES(:),ROW_INDICES(:)
    REAL(DP) :: DEPENDENT_VALUE,MATRIX_VALUE,RHS_VALUE,SOURCE_VALUE
    REAL(DP), POINTER :: DEPENDENT_PARAMETERS(:),equationsMatrixData(:),sourceVectorData(:)
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: RHS_BOUNDARY_CONDITIONS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: COLUMN_DOMAIN_MAPPING,RHS_DOMAIN_MAPPING
    TYPE(DistributedMatrixType), POINTER :: EQUATIONS_DISTRIBUTED_MATRIX
    TYPE(DistributedVectorType), POINTER :: SOURCE_DISTRIBUTED_VECTOR
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMappingSourceType), POINTER :: sourceMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE,rhsVariable
    TYPE(VARYING_STRING) :: localError

    NULLIFY(DEPENDENT_PARAMETERS)
    NULLIFY(equationsMatrixData)
    NULLIFY(sourceVectorData)

    ENTERS("EQUATIONS_SET_BACKSUBSTITUTE",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%EQUATIONS_SET_FINISHED) THEN
        dependentField=>equationsSet%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          equations=>equationsSet%EQUATIONS
          IF(ASSOCIATED(EQUATIONS)) THEN
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            vectorMatrices=>vectorEquations%vectorMatrices
            IF(ASSOCIATED(vectorMatrices)) THEN
              dynamicMatrices=>vectorMatrices%dynamicMatrices
              IF(ASSOCIATED(dynamicMatrices)) THEN
                !CALL FlagError("Not implemented.",err,error,*999)
              ELSE
                linearMatrices=>vectorMatrices%linearMatrices
                IF(ASSOCIATED(linearMatrices)) THEN
                  vectorMapping=>vectorEquations%vectorMapping
                  IF(ASSOCIATED(vectorMapping)) THEN
                    linearMapping=>vectorMapping%linearMapping
                    IF(ASSOCIATED(linearMapping)) THEN
                      rhsMapping=>vectorMapping%rhsMapping
                      sourceMapping=>vectorMapping%sourceMapping
                      IF(ASSOCIATED(rhsMapping)) THEN
                        IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                          IF(ASSOCIATED(sourceMapping)) THEN
                            sourceVector=>vectorMatrices%sourceVector
                            IF(ASSOCIATED(sourceVector)) THEN
                              SOURCE_DISTRIBUTED_VECTOR=>sourceVector%vector
                              IF(ASSOCIATED(SOURCE_DISTRIBUTED_VECTOR)) THEN
                                CALL DistributedVector_DataGet(SOURCE_DISTRIBUTED_VECTOR,sourceVectorData,err,error,*999)
                              ELSE
                                CALL FlagError("Source distributed vector is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("Source vector is not associated.",err,error,*999)
                            ENDIF
                          ENDIF
                          rhsVariable=>rhsMapping%rhsVariable
                          IF(ASSOCIATED(rhsVariable)) THEN
                            rhsVariableType=rhsVariable%VARIABLE_TYPE
                            RHS_DOMAIN_MAPPING=>rhsVariable%DOMAIN_MAPPING
                            IF(ASSOCIATED(RHS_DOMAIN_MAPPING)) THEN
                              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,rhsVariable,RHS_BOUNDARY_CONDITIONS, &
                                & err,error,*999)
                              IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                                !Loop over the equations matrices
                                DO equations_matrix_idx=1,linearMatrices%numberOfLinearMatrices
                                  DEPENDENT_VARIABLE=>linearMapping%equationsMatrixToVarMaps(equations_matrix_idx)%VARIABLE
                                  IF(ASSOCIATED(DEPENDENT_VARIABLE)) THEN
                                    VARIABLE_TYPE=DEPENDENT_VARIABLE%VARIABLE_TYPE
                                    !Get the dependent field variable parameters
                                    CALL Field_ParameterSetDataGet(dependentField,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                      & DEPENDENT_PARAMETERS,err,error,*999)
                                    equationsMatrix=>linearMatrices%matrices(equations_matrix_idx)%ptr
                                    IF(ASSOCIATED(equationsMatrix)) THEN
                                      COLUMN_DOMAIN_MAPPING=>linearMapping%equationsMatrixToVarMaps(equations_matrix_idx)% &
                                        & columnDOFSMapping
                                      IF(ASSOCIATED(COLUMN_DOMAIN_MAPPING)) THEN
                                        EQUATIONS_DISTRIBUTED_MATRIX=>equationsMatrix%MATRIX
                                        IF(ASSOCIATED(EQUATIONS_DISTRIBUTED_MATRIX)) THEN
                                          CALL DistributedMatrix_StorageTypeGet(EQUATIONS_DISTRIBUTED_MATRIX, &
                                            & EQUATIONS_STORAGE_TYPE,err,error,*999)
                                          CALL DistributedMatrix_DataGet(EQUATIONS_DISTRIBUTED_MATRIX,equationsMatrixData, &
                                            & err,error,*999)
                                          SELECT CASE(EQUATIONS_STORAGE_TYPE)
                                          CASE(DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE)
                                            !Loop over the non ghosted rows in the equations set
                                            DO equations_row_number=1,vectorMapping%numberOfRows
                                              RHS_VALUE=0.0_DP
                                              rhs_variable_dof=rhsMapping%equationsRowToRHSDOFMap(equations_row_number)
                                              rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                              rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%DOF_TYPES(rhs_global_dof)
                                              !For free RHS DOFs, set the right hand side field values by multiplying the
                                              !row by the dependent variable value
                                              SELECT CASE(rhs_boundary_condition)
                                              CASE(BOUNDARY_CONDITION_DOF_FREE)
                                                !Back substitute
                                                !Loop over the local columns of the equations matrix
                                                DO equations_column_idx=1,COLUMN_DOMAIN_MAPPING%TOTAL_NUMBER_OF_LOCAL
                                                  equations_column_number=COLUMN_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP( &
                                                    & equations_column_idx)
                                                  variable_dof=equations_column_idx
                                                  MATRIX_VALUE=equationsMatrixData(equations_row_number+ &
                                                    & (equations_column_number-1)*vectorMatrices%totalNumberOfRows)
                                                  DEPENDENT_VALUE=DEPENDENT_PARAMETERS(variable_dof)
                                                  RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                                ENDDO !equations_column_idx
                                              CASE(BOUNDARY_CONDITION_DOF_FIXED)
                                                !Do nothing
                                              CASE(BOUNDARY_CONDITION_DOF_MIXED)
                                                !Robin or is it Cauchy??? boundary conditions
                                                CALL FlagError("Not implemented.",err,error,*999)
                                              CASE DEFAULT
                                                localError="The RHS variable boundary condition of "// &
                                                  & TRIM(NumberToVString(rhs_boundary_condition,"*",err,error))// &
                                                  & " for RHS variable dof number "// &
                                                  & TRIM(NumberToVString(rhs_variable_dof,"*",err,error))//" is invalid."
                                                CALL FlagError(localError,err,error,*999)
                                              END SELECT
                                              IF(ASSOCIATED(sourceMapping)) THEN
                                                SOURCE_VALUE=sourceVectorData(equations_row_number)
                                                RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                              ENDIF
                                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,rhsVariableType, &
                                                & FIELD_VALUES_SET_TYPE,rhs_variable_dof,RHS_VALUE,err,error,*999)
                                            ENDDO !equations_row_number
                                          CASE(DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(DISTRIBUTED_MATRIX_COLUMN_MAJOR_STORAGE_TYPE)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(DISTRIBUTED_MATRIX_ROW_MAJOR_STORAGE_TYPE)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE)
                                            CALL DistributedMatrix_StorageLocationsGet(EQUATIONS_DISTRIBUTED_MATRIX, &
                                              & ROW_INDICES,COLUMN_INDICES,err,error,*999)
                                            !Loop over the non-ghosted rows in the equations set
                                            DO equations_row_number=1,vectorMapping%numberOfRows
                                              RHS_VALUE=0.0_DP
                                              rhs_variable_dof=rhsMapping%equationsRowToRHSDOFMap(equations_row_number)
                                              rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(rhs_variable_dof)
                                              rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%DOF_TYPES(rhs_global_dof)
                                              SELECT CASE(rhs_boundary_condition)
                                              CASE(BOUNDARY_CONDITION_DOF_FREE)
                                                !Back substitute
                                                !Loop over the local columns of the equations matrix
                                                DO equations_column_idx=ROW_INDICES(equations_row_number), &
                                                  ROW_INDICES(equations_row_number+1)-1
                                                  equations_column_number=COLUMN_INDICES(equations_column_idx)
                                                  variable_dof=equations_column_idx-ROW_INDICES(equations_row_number)+1
                                                  MATRIX_VALUE=equationsMatrixData(equations_column_idx)
                                                  DEPENDENT_VALUE=DEPENDENT_PARAMETERS(variable_dof)
                                                  RHS_VALUE=RHS_VALUE+MATRIX_VALUE*DEPENDENT_VALUE
                                                ENDDO !equations_column_idx
                                              CASE(BOUNDARY_CONDITION_DOF_FIXED)
                                                !Do nothing
                                              CASE(BOUNDARY_CONDITION_DOF_MIXED)
                                                !Robin or is it Cauchy??? boundary conditions
                                                CALL FlagError("Not implemented.",err,error,*999)
                                              CASE DEFAULT
                                                localError="The global boundary condition of "// &
                                                  & TRIM(NumberToVString(rhs_boundary_condition,"*",err,error))// &
                                                  & " for RHS variable dof number "// &
                                                  & TRIM(NumberToVString(rhs_variable_dof,"*",err,error))//" is invalid."
                                                CALL FlagError(localError,err,error,*999)
                                              END SELECT
                                              IF(ASSOCIATED(sourceMapping)) THEN
                                                SOURCE_VALUE=sourceVectorData(equations_row_number)
                                                RHS_VALUE=RHS_VALUE-SOURCE_VALUE
                                              ENDIF
                                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,rhsVariableType, &
                                                & FIELD_VALUES_SET_TYPE,rhs_variable_dof,RHS_VALUE,err,error,*999)
                                            ENDDO !equations_row_number
                                          CASE(DISTRIBUTED_MATRIX_COMPRESSED_COLUMN_STORAGE_TYPE)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE(DISTRIBUTED_MATRIX_ROW_COLUMN_STORAGE_TYPE)
                                            CALL FlagError("Not implemented.",err,error,*999)
                                          CASE DEFAULT
                                            localError="The matrix storage type of "// &
                                              & TRIM(NumberToVString(EQUATIONS_STORAGE_TYPE,"*",err,error))//" is invalid."
                                            CALL FlagError(localError,err,error,*999)
                                          END SELECT
                                          CALL DistributedMatrix_DataRestore(EQUATIONS_DISTRIBUTED_MATRIX,equationsMatrixData, &
                                            & err,error,*999)
                                        ELSE
                                          CALL FlagError("Equations matrix distributed matrix is not associated.",err,error,*999)
                                        ENDIF
                                      ELSE
                                        CALL FlagError("Equations column domain mapping is not associated.",err,error,*999)
                                      ENDIF
                                    ELSE
                                      CALL FlagError("Equations equations matrix is not associated.",err,error,*999)
                                    ENDIF
                                    !Restore the dependent field variable parameters
                                    CALL Field_ParameterSetDataRestore(dependentField,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                      & DEPENDENT_PARAMETERS,err,error,*999)
                                  ELSE
                                    CALL FlagError("Dependent variable is not associated.",err,error,*999)
                                  ENDIF
                                ENDDO !equations_matrix_idx
                                !Start the update of the field parameters
                                CALL Field_ParameterSetUpdateStart(dependentField,rhsVariableType,FIELD_VALUES_SET_TYPE, &
                                  & err,error,*999)
                                !Finish the update of the field parameters
                                CALL Field_ParameterSetUpdateFinish(dependentField,rhsVariableType,FIELD_VALUES_SET_TYPE, &
                                  & err,error,*999)
                              ELSE
                                CALL FlagError("RHS boundary conditions variable is not associated.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("RHS variable domain mapping is not associated.",err,error,*999)
                            ENDIF
                          ELSE
                            CALL FlagError("RHS variable is not associated.",err,error,*999)
                          ENDIF
                          IF(ASSOCIATED(sourceMapping)) THEN
                            CALL DistributedVector_DataRestore(SOURCE_DISTRIBUTED_VECTOR,sourceVectorData,err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("Boundary conditions are not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Equations mapping RHS mappings is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Equations mapping linear mapping is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations mapping is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("Equations matrices linear matrices is not associated.",err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations matrices is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Dependent field is not associated.",err,error,*999)
        ENDIF
      ELSE            
        CALL FlagError("Equations set has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
          
    EXITS("EQUATIONS_SET_BACKSUBSTITUTE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_BACKSUBSTITUTE",err,error)
    RETURN 1
   
  END SUBROUTINE EQUATIONS_SET_BACKSUBSTITUTE
  
  !
  !================================================================================================================================
  !

  !>Updates the right hand side variable from the equations residual vector
  SUBROUTINE EQUATIONS_SET_NONLINEAR_RHS_UPDATE(equationsSet,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<Boundary conditions to use for the RHS update
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: variable_dof,row_idx,VARIABLE_TYPE,rhs_global_dof,rhs_boundary_condition,equations_matrix_idx
    REAL(DP) :: VALUE
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMappingRHSType), POINTER :: rhsMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(DistributedVectorType), POINTER :: residualVector
    TYPE(FIELD_TYPE), POINTER :: RHS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rhsVariable,residualVariable
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: RHS_BOUNDARY_CONDITIONS
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: RHS_DOMAIN_MAPPING
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_NONLINEAR_RHS_UPDATE",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      equations=>equationsSet%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        vectorMapping=>vectorEquations%vectorMapping
        IF(ASSOCIATED(vectorMapping)) THEN
          rhsMapping=>vectorMapping%rhsMapping
          IF(ASSOCIATED(rhsMapping)) THEN
            rhsVariable=>rhsMapping%rhsVariable
            IF(ASSOCIATED(rhsVariable)) THEN
              !Get the right hand side variable
              RHS_FIELD=>rhsVariable%FIELD
              VARIABLE_TYPE=rhsVariable%VARIABLE_TYPE
            ELSE
              CALL FlagError("RHS mapping RHS variable is not associated.",err,error,*999)
            ENDIF
          ELSE
            CALL FlagError("Equations mapping RHS mapping is not associated.",err,error,*999)
          ENDIF
          IF(ASSOCIATED(RHS_FIELD)) THEN
            IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
              RHS_DOMAIN_MAPPING=>rhsVariable%DOMAIN_MAPPING
              IF(ASSOCIATED(RHS_DOMAIN_MAPPING)) THEN
                CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,rhsVariable,RHS_BOUNDARY_CONDITIONS, &
                  & err,error,*999)
                IF(ASSOCIATED(RHS_BOUNDARY_CONDITIONS)) THEN
                  !Get the equations residual vector
                  vectorMatrices=>vectorEquations%vectorMatrices
                  IF(ASSOCIATED(vectorMatrices)) THEN
                    nonlinearMatrices=>vectorMatrices%nonlinearMatrices
                    IF(ASSOCIATED(nonlinearMatrices)) THEN
                      residualVector=>nonlinearMatrices%residual
                      IF(ASSOCIATED(residualVector)) THEN
                        !Get mapping from equations rows to field dofs
                        nonlinearMapping=>vectorMapping%nonlinearMapping
                        IF(ASSOCIATED(nonlinearMapping)) THEN
                          DO equations_matrix_idx=1,nonlinearMapping%numberOfResidualVariables
                            residualVariable=>nonlinearMapping%jacobianToVarMap(equations_matrix_idx)%VARIABLE
                            IF(ASSOCIATED(residualVariable)) THEN
                              DO row_idx=1,vectorMapping%numberOfRows
                                variable_dof=rhsMapping%equationsRowToRHSDOFMap(row_idx)
                                rhs_global_dof=RHS_DOMAIN_MAPPING%LOCAL_TO_GLOBAL_MAP(variable_dof)
                                rhs_boundary_condition=RHS_BOUNDARY_CONDITIONS%DOF_TYPES(rhs_global_dof)
                                SELECT CASE(rhs_boundary_condition)
                                CASE(BOUNDARY_CONDITION_DOF_FREE)
                                  !Add residual to field value
                                  CALL DistributedVector_ValuesGet(residualVector,row_idx,VALUE,err,error,*999)
                                  CALL Field_ParameterSetUpdateLocalDOF(RHS_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                                    & variable_dof,VALUE,err,error,*999)
                                CASE(BOUNDARY_CONDITION_DOF_FIXED)
                                  !Do nothing
                                CASE(BOUNDARY_CONDITION_DOF_MIXED)
                                  CALL FlagError("Not implemented.",err,error,*999)
                                CASE DEFAULT
                                  localError="The RHS variable boundary condition of "// &
                                    & TRIM(NumberToVString(rhs_boundary_condition,"*",err,error))// &
                                    & " for RHS variable dof number "// &
                                    & TRIM(NumberToVString(variable_dof,"*",err,error))//" is invalid."
                                  CALL FlagError(localError,err,error,*999)
                                END SELECT
                              ENDDO
                            ELSE
                              CALL FlagError("Residual variable is not associated.",err,error,*999)
                            ENDIF
                          ENDDO !equations_matrix_idx
                        ELSE
                          CALL FlagError("Nonlinear mapping is not associated.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Residual vector is not associated.",err,error,*999)
                      ENDIF
                    ELSE
                      CALL FlagError("Nonlinear matrices is not associated.",err,error,*999)
                    ENDIF
                  ELSE
                    CALL FlagError("Equations matrices is not associated.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("RHS boundary conditions variable is not associated.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("RHS variable domain mapping is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Boundary conditions are not associated.",err,error,*999)
            ENDIF
            CALL Field_ParameterSetUpdateStart(RHS_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
            CALL Field_ParameterSetUpdateFinish(RHS_FIELD,VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,err,error,*999)
          ELSE
            CALL FlagError("RHS variable field is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations mapping is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SET_NONLINEAR_RHS_UPDATE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_NONLINEAR_RHS_UPDATE",err,error)
    RETURN 1

  END SUBROUTINE EQUATIONS_SET_NONLINEAR_RHS_UPDATE

  !
  !================================================================================================================================
  !

  !>Set boundary conditions for an equation set according to the analytic equations. \see OpenCMISS::cmfe_EquationsSet_BoundaryConditionsAnalytic
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC(equationsSet,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the analyticboundary conditions for.
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<A pointer to the boundary conditions to set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(.NOT.ALLOCATED(equationsSet%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(equationsSet%SPECIFICATION,1)<1) THEN
        CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
      END IF
      IF(equationsSet%DEPENDENT%DEPENDENT_FINISHED) THEN
        IF(ASSOCIATED(equationsSet%ANALYTIC)) THEN
          IF(equationsSet%ANALYTIC%ANALYTIC_FINISHED) THEN
            SELECT CASE(equationsSet%SPECIFICATION(1))
            CASE(EQUATIONS_SET_ELASTICITY_CLASS)
              CALL Elasticity_BoundaryConditionsAnalyticCalculate(equationsSet,BOUNDARY_CONDITIONS,err,error,*999)
            CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
              CALL FluidMechanics_BoundaryConditionsAnalyticCalculate(equationsSet,BOUNDARY_CONDITIONS,err,error,*999)
            CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
              CALL ClassicalField_BoundaryConditionsAnalyticCalculate(equationsSet,BOUNDARY_CONDITIONS,err,error,*999)
            CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_MODAL_CLASS)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
              CALL FlagError("Not implemented.",err,error,*999)
            CASE DEFAULT
              localError="The first equations set specification of "//TRIM(NumberToVString(equationsSet%SPECIFICATION(1),"*", &
                & err,error))//" is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Equations set analytic has not been finished.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Equations set analytic is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Equations set dependent has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_ANALYTIC

  !
  !================================================================================================================================
  !

  !>Finishes the process of creating an equation set on a region. \see OpenCMISS::cmfe_EquationsSet_CreateStart
  SUBROUTINE EQUATIONS_SET_CREATE_FINISH(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to finish creating
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    
    ENTERS("EQUATIONS_SET_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%EQUATIONS_SET_FINISHED) THEN
        CALL FlagError("Equations set has already been finished.",err,error,*999)
      ELSE            
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INITIAL_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        !Finish the equations set specific setup
        CALL EQUATIONS_SET_SETUP(equationsSet,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        !Finish the equations set specific geometry setup
        CALL EQUATIONS_SET_SETUP(equationsSet,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
        !Finish the equations set creation
        equationsSet%EQUATIONS_SET_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_CREATE_FINISH",err,error)
    RETURN 1
   
  END SUBROUTINE EQUATIONS_SET_CREATE_FINISH
        
  !
  !================================================================================================================================
  !

  !>Starts the process of creating an equations set defined by USER_NUMBER in the region identified by REGION. \see OpenCMISS::cmfe_EquationsSet_CreateStart
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
  SUBROUTINE EQUATIONS_SET_CREATE_START(USER_NUMBER,REGION,GEOM_FIBRE_FIELD,EQUATIONS_SET_SPECIFICATION,&
      & EQUATIONS_SET_FIELD_USER_NUMBER,EQUATIONS_SET_FIELD_FIELD,EQUATIONS_SET,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the equations set
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to create the equations set on
    TYPE(FIELD_TYPE), POINTER :: GEOM_FIBRE_FIELD !<A pointer to the either the geometry or, if appropriate, the fibre field for the equation set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_SPECIFICATION(:) !<The equations set specification array to set
    INTEGER(INTG), INTENT(IN) :: EQUATIONS_SET_FIELD_USER_NUMBER !<The user number of the equations set field
    TYPE(FIELD_TYPE), POINTER :: EQUATIONS_SET_FIELD_FIELD !<On return, a pointer to the equations set field
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<On return, a pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr,equations_set_idx
    TYPE(EQUATIONS_SET_TYPE), POINTER :: NEW_EQUATIONS_SET
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(REGION_TYPE), POINTER :: GEOM_FIBRE_FIELD_REGION,EQUATIONS_SET_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE), POINTER :: EQUATIONS_EQUATIONS_SET_FIELD
    TYPE(FIELD_TYPE), POINTER :: FIELD

    NULLIFY(NEW_EQUATIONS_SET)
    NULLIFY(NEW_EQUATIONS_SETS)
    NULLIFY(EQUATIONS_EQUATIONS_SET_FIELD)

    ENTERS("EQUATIONS_SET_CREATE_START",err,error,*997)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        CALL EQUATIONS_SET_USER_NUMBER_FIND(USER_NUMBER,REGION,NEW_EQUATIONS_SET,err,error,*997)
        IF(ASSOCIATED(NEW_EQUATIONS_SET)) THEN
          localError="Equations set user number "//TRIM(NumberToVString(USER_NUMBER,"*",err,error))// &
            & " has already been created on region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*997)
        ELSE
          NULLIFY(NEW_EQUATIONS_SET)
          IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
            IF(GEOM_FIBRE_FIELD%FIELD_FINISHED) THEN
              IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE.OR.GEOM_FIBRE_FIELD%TYPE==FIELD_FIBRE_TYPE) THEN
                GEOM_FIBRE_FIELD_REGION=>GEOM_FIBRE_FIELD%REGION
                IF(ASSOCIATED(GEOM_FIBRE_FIELD_REGION)) THEN
                  IF(GEOM_FIBRE_FIELD_REGION%USER_NUMBER==REGION%USER_NUMBER) THEN
                      IF(ASSOCIATED(EQUATIONS_SET_FIELD_FIELD)) THEN
                        !Check the equations set field has been finished
                        IF(EQUATIONS_SET_FIELD_FIELD%FIELD_FINISHED.eqv..TRUE.) THEN
                          !Check the user numbers match
                          IF(EQUATIONS_SET_FIELD_USER_NUMBER/=EQUATIONS_SET_FIELD_FIELD%USER_NUMBER) THEN
                            localError="The specified equations set field user number of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_FIELD_USER_NUMBER,"*",err,error))// &
                              & " does not match the user number of the specified equations set field of "// &
                              & TRIM(NumberToVString(EQUATIONS_SET_FIELD_FIELD%USER_NUMBER,"*",err,error))//"."
                            CALL FlagError(localError,err,error,*999)
                          ENDIF
                          EQUATIONS_SET_FIELD_REGION=>EQUATIONS_SET_FIELD_FIELD%REGION
                          IF(ASSOCIATED(EQUATIONS_SET_FIELD_REGION)) THEN                
                            !Check the field is defined on the same region as the equations set
                            IF(EQUATIONS_SET_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                              localError="Invalid region setup. The specified equations set field was created on region no. "// &
                                & TRIM(NumberToVString(EQUATIONS_SET_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                                & " and the specified equations set has been created on region number "// &
                                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                              CALL FlagError(localError,err,error,*999)
                            ENDIF
                            !Check the specified equations set field has the same decomposition as the geometric field
                            IF(ASSOCIATED(GEOM_FIBRE_FIELD)) THEN
                              IF(.NOT.ASSOCIATED(GEOM_FIBRE_FIELD%DECOMPOSITION,EQUATIONS_SET_FIELD_FIELD%DECOMPOSITION)) THEN
                                CALL FlagError("The specified equations set field does not have the same decomposition "// &
                                  & "as the geometric field for the specified equations set.",err,error,*999)
                              ENDIF
                            ELSE
                              CALL FlagError("The geom. field is not associated for the specified equations set.",err,error,*999)
                            ENDIF
                              
                          ELSE
                            CALL FlagError("The specified equations set field region is not associated.",err,error,*999)
                          ENDIF
                        ELSE
                          CALL FlagError("The specified equations set field has not been finished.",err,error,*999)
                        ENDIF
                      ELSE
                        !Check the user number has not already been used for a field in this region.
                        NULLIFY(FIELD)
                        CALL FIELD_USER_NUMBER_FIND(EQUATIONS_SET_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
                        IF(ASSOCIATED(FIELD)) THEN
                          localError="The specified equations set field user number of "// &
                            & TRIM(NumberToVString(EQUATIONS_SET_FIELD_USER_NUMBER,"*",err,error))// &
                            & "has already been used to create a field on region number "// &
                            & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                          CALL FlagError(localError,err,error,*999)
                        ENDIF
                      ENDIF
                      !Initalise equations set
                      CALL EQUATIONS_SET_INITIALISE(NEW_EQUATIONS_SET,err,error,*999)
                      !Set default equations set values
                      NEW_EQUATIONS_SET%USER_NUMBER=USER_NUMBER
                      NEW_EQUATIONS_SET%GLOBAL_NUMBER=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1
                      NEW_EQUATIONS_SET%EQUATIONS_SETS=>REGION%EQUATIONS_SETS
                      NEW_EQUATIONS_SET%label="Equations Set "//TRIM(NumberToVString(USER_NUMBER,"*",err,error))
                      NEW_EQUATIONS_SET%REGION=>REGION
                      !Set the equations set class, type and subtype
                      CALL EquationsSet_SpecificationSet(NEW_EQUATIONS_SET,EQUATIONS_SET_SPECIFICATION,err,error,*999)
                      NEW_EQUATIONS_SET%EQUATIONS_SET_FINISHED=.FALSE.
                      !Initialise the setup
                      CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
                      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INITIAL_TYPE
                      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
                      !Here, we get a pointer to the equations_set_field; default is null
                      EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=EQUATIONS_SET_FIELD_USER_NUMBER
                      EQUATIONS_SET_SETUP_INFO%FIELD=>EQUATIONS_SET_FIELD_FIELD
                      !Start equations set specific setup
                      CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
                      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
                      !Set up the equations set geometric fields
                      CALL EQUATIONS_SET_GEOMETRY_INITIALISE(NEW_EQUATIONS_SET,err,error,*999)
                      IF(GEOM_FIBRE_FIELD%TYPE==FIELD_GEOMETRIC_TYPE) THEN
                        NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD
                        NULLIFY(NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)
                      ELSE
                        NEW_EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD=>GEOM_FIBRE_FIELD%GEOMETRIC_FIELD
                        NEW_EQUATIONS_SET%GEOMETRY%FIBRE_FIELD=>GEOM_FIBRE_FIELD
                      ENDIF
                      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_GEOMETRY_TYPE
                      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
                      EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=GEOM_FIBRE_FIELD%USER_NUMBER
                      EQUATIONS_SET_SETUP_INFO%FIELD=>GEOM_FIBRE_FIELD
                      !Set up equations set specific geometry
                      CALL EQUATIONS_SET_SETUP(NEW_EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
                      !Finalise the setup
                      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
                      !Add new equations set into list of equations set in the region
                      ALLOCATE(NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1),STAT=ERR)
                      IF(ERR/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
                      DO equations_set_idx=1,REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
                        NEW_EQUATIONS_SETS(equations_set_idx)%ptr=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr
                      ENDDO !equations_set_idx
                      NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1)%ptr=>NEW_EQUATIONS_SET
                      IF(ASSOCIATED(REGION%EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
                      REGION%EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
                      REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS+1
                      EQUATIONS_SET=>NEW_EQUATIONS_SET
                      EQUATIONS_EQUATIONS_SET_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD
                      !\todo check pointer setup
                      IF(EQUATIONS_EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED) THEN
                        EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD
                      ELSE
                        EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD=>EQUATIONS_SET_FIELD_FIELD
                      ENDIF
                  ELSE
                    localError="The geometric field region and the specified region do not match. "// &
                      & "The geometric field was created on region number "// &
                      & TRIM(NumberToVString(GEOM_FIBRE_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                      & " and the specified region number is "// &
                      & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                    CALL FlagError(localError,err,error,*997)
                  ENDIF
                ELSE
                  CALL FlagError("The specified geometric fields region is not associated.",err,error,*997)
                ENDIF
              ELSE
                CALL FlagError("The specified geometric field is not a geometric or fibre field.",err,error,*997)
              ENDIF
            ELSE
              CALL FlagError("The specified geometric field is not finished.",err,error,*997)
            ENDIF
          ELSE
            CALL FlagError("The specified geometric field is not associated.",err,error,*997)
          ENDIF
        ENDIF
      ELSE
        localError="The equations sets on region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))// &
          & " are not associated."
        CALL FlagError(localError,err,error,*997)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*997)
    ENDIF
    
    EXITS("EQUATIONS_SET_CREATE_START")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SET))CALL EQUATIONS_SET_FINALISE(NEW_EQUATIONS_SET,dummyErr,dummyError,*998)
998 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
997 ERRORSEXITS("EQUATIONS_SET_CREATE_START",err,error)
    RETURN 1   
  END SUBROUTINE EQUATIONS_SET_CREATE_START
  
  !
  !================================================================================================================================
  !

  !>Destroys an equations set identified by a user number on the give region and deallocates all memory. \see OpenCMISS::cmfe_EquationsSet_Destroy
  SUBROUTINE EQUATIONS_SET_DESTROY_NUMBER(USER_NUMBER,REGION,err,error,*)

    !Argument variables
    INTEGER(INTG), INTENT(IN) :: USER_NUMBER !<The user number of the equations set to destroy
    TYPE(REGION_TYPE), POINTER :: REGION !<The region of the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    LOGICAL :: FOUND
    TYPE(VARYING_STRING) :: localError
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)

    NULLIFY(NEW_EQUATIONS_SETS)

    ENTERS("EQUATIONS_SET_DESTROY_NUMBER",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        
        !Find the equations set identified by the user number
        FOUND=.FALSE.
        equations_set_position=0
        DO WHILE(equations_set_position<REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS.AND..NOT.FOUND)
          equations_set_position=equations_set_position+1
          IF(REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_position)%ptr%USER_NUMBER==USER_NUMBER)FOUND=.TRUE.
        ENDDO
        
        IF(FOUND) THEN
          
          EQUATIONS_SET=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_position)%ptr
          
          !Destroy all the equations set components
          CALL EQUATIONS_SET_FINALISE(EQUATIONS_SET,err,error,*999)
          
          !Remove the equations set from the list of equations set
          IF(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>1) THEN
            ALLOCATE(NEW_EQUATIONS_SETS(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1),STAT=ERR)
            IF(ERR/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
            DO equations_set_idx=1,REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
              IF(equations_set_idx<equations_set_position) THEN
                NEW_EQUATIONS_SETS(equations_set_idx)%ptr=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr
              ELSE IF(equations_set_idx>equations_set_position) THEN
                REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr%GLOBAL_NUMBER=REGION%EQUATIONS_SETS% &
                  & EQUATIONS_SETS(equations_set_idx)%ptr%GLOBAL_NUMBER-1
                NEW_EQUATIONS_SETS(equations_set_idx-1)%ptr=>REGION%EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr
              ENDIF
            ENDDO !equations_set_idx
            IF(ASSOCIATED(REGION%EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
            REGION%EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
            REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1
          ELSE
            DEALLOCATE(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
            REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
          ENDIF
          
        ELSE
          localError="Equations set number "//TRIM(NumberToVString(USER_NUMBER,"*",err,error))// &
            & " has not been created on region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        ENDIF
      ELSE
        localError="The equations sets on region number "//TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))// &
          & " are not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*998)
    ENDIF    

    EXITS("EQUATIONS_SET_DESTROY_NUMBER")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
998 ERRORSEXITS("EQUATIONS_SET_DESTROY_NUMBER",err,error)
    RETURN 1   
  END SUBROUTINE EQUATIONS_SET_DESTROY_NUMBER
  
  !
  !================================================================================================================================
  !

  !>Destroys an equations set identified by a pointer and deallocates all memory. \see OpenCMISS::cmfe_EquationsSet_Destroy
  SUBROUTINE EQUATIONS_SET_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx,equations_set_position
    TYPE(EQUATIONS_SETS_TYPE), POINTER :: EQUATIONS_SETS
    TYPE(EQUATIONS_SET_PTR_TYPE), POINTER :: NEW_EQUATIONS_SETS(:)

    NULLIFY(NEW_EQUATIONS_SETS)

    ENTERS("EQUATIONS_SET_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SETS=>EQUATIONS_SET%EQUATIONS_SETS
      IF(ASSOCIATED(EQUATIONS_SETS)) THEN
        equations_set_position=EQUATIONS_SET%GLOBAL_NUMBER

        !Destroy all the equations set components
        CALL EQUATIONS_SET_FINALISE(EQUATIONS_SET,err,error,*999)
        
        !Remove the equations set from the list of equations set
        IF(EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>1) THEN
          ALLOCATE(NEW_EQUATIONS_SETS(EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1),STAT=ERR)
          IF(ERR/=0) CALL FlagError("Could not allocate new equations sets.",err,error,*999)
          DO equations_set_idx=1,EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS
            IF(equations_set_idx<equations_set_position) THEN
              NEW_EQUATIONS_SETS(equations_set_idx)%ptr=>EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr
            ELSE IF(equations_set_idx>equations_set_position) THEN
              EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr%GLOBAL_NUMBER=EQUATIONS_SETS% &
                & EQUATIONS_SETS(equations_set_idx)%ptr%GLOBAL_NUMBER-1
              NEW_EQUATIONS_SETS(equations_set_idx-1)%ptr=>EQUATIONS_SETS%EQUATIONS_SETS(equations_set_idx)%ptr
            ENDIF
          ENDDO !equations_set_idx
          IF(ASSOCIATED(EQUATIONS_SETS%EQUATIONS_SETS)) DEALLOCATE(EQUATIONS_SETS%EQUATIONS_SETS)
          EQUATIONS_SETS%EQUATIONS_SETS=>NEW_EQUATIONS_SETS
          EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS-1
        ELSE
          DEALLOCATE(EQUATIONS_SETS%EQUATIONS_SETS)
          EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
        ENDIF
        
      ELSE
        CALL FlagError("Equations set equations set is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF    

    EXITS("EQUATIONS_SET_DESTROY")
    RETURN
999 IF(ASSOCIATED(NEW_EQUATIONS_SETS)) DEALLOCATE(NEW_EQUATIONS_SETS)
998 ERRORSEXITS("EQUATIONS_SET_DESTROY",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalise the equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_FINALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finalise.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      CALL EQUATIONS_SET_GEOMETRY_FINALISE(EQUATIONS_SET%GEOMETRY,err,error,*999)
      CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,err,error,*999)
      CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,err,error,*999)
      CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,err,error,*999)
      CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,err,error,*999)
      CALL EQUATIONS_SET_ANALYTIC_FINALISE(EQUATIONS_SET%ANALYTIC,err,error,*999)
      CALL EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE(EQUATIONS_SET%EQUATIONS_SET_FIELD,err,error,*999)
      CALL EquationsSet_DerivedFinalise(EQUATIONS_SET%derived,err,error,*999)
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) CALL Equations_Destroy(EQUATIONS_SET%EQUATIONS,err,error,*999)
      IF(ALLOCATED(EQUATIONS_SET%SPECIFICATION)) DEALLOCATE(EQUATIONS_SET%SPECIFICATION)
      EQUATIONS_SET%label=""
      DEALLOCATE(EQUATIONS_SET)
    ENDIF
       
    EXITS("EQUATIONS_SET_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_FINALISE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_FINALISE

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matries and rhs vector for the given element number for a finite element equations set.
  SUBROUTINE EquationsSet_FiniteElementCalculate(equationsSet,elementNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
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
    TYPE(EquationsMatricesSourceType), POINTER :: sourceVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(VARYING_STRING) :: localError
    
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_START("EquationsSet_FiniteElementCalculate()")
#endif

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
      CALL ELASTICITY_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
      CALL FLUID_MECHANICS_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
    CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
      CALL FlagError("Not implemented.",err,error,*999)
    CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
      CALL CLASSICAL_FIELD_FINITE_ELEMENT_CALCULATE(equationsSet,elementNumber,err,error,*999)
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
 
#ifdef TAUPROF
    CALL TAU_STATIC_PHASE_STOP("EquationsSet_FiniteElementCalculate()")
#endif
       
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
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

    DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
      SELECT CASE(nonlinearMatrices%jacobians(matrixIdx)%ptr%jacobianCalculationType)
      CASE(EQUATIONS_JACOBIAN_ANALYTIC_CALCULATED)
        ! None of these routines currently support calculating off diagonal terms for coupled problems,
        ! but when one does we will have to pass through the matrixIdx parameter
        IF(matrixIdx>1) CALL FlagError("Analytic off-diagonal Jacobian calculation not implemented.",err,error,*999)
        SELECT CASE(equationsSet%specification(1))
        CASE(EQUATIONS_SET_ELASTICITY_CLASS)
          CALL ELASTICITY_FINITE_ELEMENT_JACOBIAN_EVALUATE(equationsSet,elementNumber,err,error,*999)
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
        CALL EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,elementNumber,matrixIdx,err,error,*999)
      CASE DEFAULT
        localError="The Jacobian calculation type of "//TRIM(NumberToVString(nonlinearMatrices%jacobians(matrixIdx)%ptr% &
          & jacobianCalculationType,"*",err,error))//" is not valid for matrix index "// &
          & TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    END DO !matrixIdx
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
  SUBROUTINE EquationsSet_FiniteElementJacobianEvaluateFD(equationsSet,elementNumber,jacobianNumber,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet  !<A pointer to the equations set to evaluate the element Jacobian for
    INTEGER(INTG), INTENT(IN) :: elementNumber  !<The element number to calculate the Jacobian for
    INTEGER(INTG), INTENT(IN) :: jacobianNumber  !<The Jacobian number to calculate when there are coupled problems
    INTEGER(INTG), INTENT(OUT) :: err  !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error  !<The error string
    !Local Variables
    TYPE(BASIS_TYPE), POINTER :: basis
    TYPE(DOMAIN_ELEMENTS_TYPE), POINTER :: elementsTopology
    TYPE(DistributedVectorType), POINTER :: parameters
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: rowVariable,columnVariable
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
    NULLIFY(nonlinearMapping)
    CALL EquationsMappingVector_NonlinearMappingGet(vectorMapping,nonlinearMapping,err,error,*999)
    
    ! The first residual variable is always the row variable, which is the variable the
    ! residual is calculated for
    rowVariable=>nonlinearMapping%residualVariables(1)%ptr
    ! For coupled problems this routine will be called multiple times if multiple Jacobians use finite
    ! differencing, so make sure we only calculate the residual vector once, to save time and because
    ! it would otherwise add together
    IF(nonlinearMatrices%elementResidualCalculated/=elementNumber) &
      & CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
    ! Make a temporary copy of the unperturbed residuals
    elementVector=nonlinearMatrices%elementResidual
    IF(jacobianNumber<=nonlinearMatrices%numberOfJacobians) THEN
      ! For coupled nonlinear problems there will be multiple Jacobians
      ! For this equations set, we calculate the residual for the row variable
      ! while pertubing parameters from the column variable.
      ! For non coupled problems these two variables will be the same
      columnVariable=>nonlinearMapping%residualVariables(jacobianNumber)%ptr
      parameters=>columnVariable%PARAMETER_SETS%PARAMETER_SETS(FIELD_VALUES_SET_TYPE)%ptr%PARAMETERS  
      numberOfRows=nonlinearMatrices%jacobians(jacobianNumber)%ptr%elementJacobian%numberOfRows
      IF(numberOfRows/=nonlinearMatrices%elementResidual%numberOfRows) &
        & CALL FlagError("Element matrix number of rows does not match element residual vector size.",err,error,*999)
      ! determine step size
      CALL DistributedVector_L2Norm(parameters,delta,err,error,*999)
      delta=(1.0_DP+delta)*1.0E-6_DP
      ! the actual finite differencing algorithm is about 4 lines but since the parameters are all
      ! distributed out, have to use proper field accessing routines..
      ! so let's just loop over component, node/el, derivative
      column=0  ! element jacobian matrix column number
      DO componentIdx=1,columnVariable%NUMBER_OF_COMPONENTS
        elementsTopology=>columnVariable%COMPONENTS(componentIdx)%DOMAIN%TOPOLOGY%ELEMENTS
        componentInterpolationType=columnVariable%COMPONENTS(componentIdx)%INTERPOLATION_TYPE
        SELECT CASE(componentInterpolationType)
        CASE (FIELD_NODE_BASED_INTERPOLATION)
          basis=>elementsTopology%ELEMENTS(elementNumber)%BASIS
          DO nodeIdx=1,basis%NUMBER_OF_NODES
            node=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_NODES(nodeIdx)
            DO derivativeIdx=1,basis%NUMBER_OF_DERIVATIVES(nodeIdx)
              derivative=elementsTopology%ELEMENTS(elementNumber)%ELEMENT_DERIVATIVES(derivativeIdx,nodeIdx)
              version=elementsTopology%ELEMENTS(elementNumber)%elementVersions(derivativeIdx,nodeIdx)
              localDOF=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP%NODES(node)% &
                & DERIVATIVES(derivative)%VERSIONS(version)
              ! one-sided finite difference
              CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
              CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
              nonlinearMatrices%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
              CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
              CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
              column=column+1
              nonlinearMatrices%jacobians(jacobianNumber)%ptr%elementJacobian%matrix(1:numberOfRows,column)= &
                & (nonlinearMatrices%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
            ENDDO !derivativeIdx
          ENDDO !nodeIdx
        CASE (FIELD_ELEMENT_BASED_INTERPOLATION)
          localDOF=columnVariable%COMPONENTS(componentIdx)%PARAM_TO_DOF_MAP%ELEMENT_PARAM2DOF_MAP%ELEMENTS(elementNumber)
          ! one-sided finite difference
          CALL DistributedVector_ValuesGet(parameters,localDOF,origDepVar,err,error,*999)
          CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar+delta,err,error,*999)
          nonlinearMatrices%elementResidual%vector=0.0_DP ! must remember to flush existing results, otherwise they're added
          CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,elementNumber,err,error,*999)
          CALL DistributedVector_ValuesSet(parameters,localDOF,origDepVar,err,error,*999)
          column=column+1
          nonlinearMatrices%jacobians(jacobianNumber)%ptr%elementJacobian%matrix(1:numberOfRows,column)= &
            & (nonlinearMatrices%elementResidual%vector(1:numberOfRows)-elementVector%vector(1:numberOfRows))/delta
        CASE DEFAULT
          CALL FlagError("Unsupported type of interpolation.",err,error,*999)
        END SELECT
      END DO !componentIdx
      ! put the original residual back in
      nonlinearMatrices%elementResidual=elementVector
    ELSE
      CALL FlagError("Invalid Jacobian number of "//TRIM(NumberToVString(jacobianNumber,"*",err,error))// &
        & ". The number should be <= "//TRIM(NumberToVString(nonlinearMatrices%numberOfJacobians,"*",err,error))// &
        & ".",err,error,*999)
    END IF

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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
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
      CALL ELASTICITY_FINITE_ELEMENT_RESIDUAL_EVALUATE(equationsSet,elementNumber,err,error,*999)
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

  !>Finish the creation of independent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_IndependentCreateFinish
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_FINISH(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the independent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: independentField

    ENTERS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED) THEN
          CALL FlagError("Equations set independent field has already been finished.",err,error,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
          IF(ASSOCIATED(independentField)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=independentField%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>independentField
            !Finish equations set specific startup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          ELSE
            CALL FlagError("Equations set independent independent field is not associated.",err,error,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finish independent creation
          EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FlagError("The equations set independent is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_INDEPENDENT_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of independent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_IndependentCreateStart
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_START(EQUATIONS_SET,INDEPENDENT_FIELD_USER_NUMBER,independentField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: INDEPENDENT_FIELD_USER_NUMBER !<The user specified independent field number
    TYPE(FIELD_TYPE), POINTER :: independentField !<If associated on entry, a pointer to the user created independent field which has the same user number as the specified independent field user number. If not associated on entry, on exit, a pointer to the created independent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,INDEPENDENT_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EQUATIONS_SET_INDEPENDENT_CREATE_START",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL FlagError("The equations set independent is already associated",err,error,*998)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(independentField)) THEN
            !Check the independent field has been finished
            IF(independentField%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(INDEPENDENT_FIELD_USER_NUMBER/=independentField%USER_NUMBER) THEN
                localError="The specified independent field user number of "// &
                  & TRIM(NumberToVString(INDEPENDENT_FIELD_USER_NUMBER,"*",err,error))// &
                  & " does not match the user number of the specified independent field of "// &
                  & TRIM(NumberToVString(independentField%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              INDEPENDENT_FIELD_REGION=>independentField%REGION
              IF(ASSOCIATED(INDEPENDENT_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(INDEPENDENT_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified independent field has been created on region number "// &
                    & TRIM(NumberToVString(INDEPENDENT_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Check the specified independent field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,independentField%DECOMPOSITION)) THEN
                    CALL FlagError("The specified independent field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified independent field region is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The specified independent field has not been finished.",err,error,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(INDEPENDENT_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
            IF(ASSOCIATED(FIELD)) THEN
              localError="The specified independent field user number of "// &
                & TRIM(NumberToVString(INDEPENDENT_FIELD_USER_NUMBER,"*",err,error))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
          !Initialise the equations set independent
          CALL EQUATIONS_SET_INDEPENDENT_INITIALISE(EQUATIONS_SET,err,error,*999)
          IF(.NOT.ASSOCIATED(independentField)) EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_INDEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=INDEPENDENT_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>independentField
          !Start equations set specific startup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Set pointers
          IF(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED) THEN            
            independentField=>EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD
          ELSE
            EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD=>independentField
          ENDIF
        ELSE
          CALL FlagError("Equation set region is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_INDEPENDENT_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_INDEPENDENT_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the independent field for an equations set. \see OpenCMISS::cmfe_EquationsSet_IndependentDestroy
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the independent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_INDEPENDENT_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,err,error,*999)
      ELSE
        CALL FlagError("Equations set indpendent is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_INDEPENDENT_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_INDEPENDENT_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the independent field for an equations set.
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET_INDEPENDENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_INDEPENDENT_TYPE), POINTER :: EQUATIONS_SET_INDEPENDENT !<A pointer to the equations set independent to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_INDEPENDENT_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET_INDEPENDENT)) THEN
      DEALLOCATE(EQUATIONS_SET_INDEPENDENT)
    ENDIF
       
    EXITS("EQUATIONS_SET_INDEPENDENT_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_INDEPENDENT_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the independent field for an equations set.
  SUBROUTINE EQUATIONS_SET_INDEPENDENT_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the independent for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EQUATIONS_SET_INDEPENDENT_INITIALISE",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%INDEPENDENT)) THEN
        CALL FlagError("Independent field is already associated for these equations sets.",err,error,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%INDEPENDENT,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate equations set independent field.",err,error,*999)
        EQUATIONS_SET%INDEPENDENT%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FINISHED=.FALSE.
        EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%INDEPENDENT%INDEPENDENT_FIELD)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_INDEPENDENT_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_INDEPENDENT_FINALISE(EQUATIONS_SET%INDEPENDENT,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_INDEPENDENT_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INDEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Initialises an equations set.
  SUBROUTINE EQUATIONS_SET_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<The pointer to the equations set to initialise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
 
    ENTERS("EQUATIONS_SET_INITIALISE",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      CALL FlagError("Equations set is already associated.",err,error,*998)
    ELSE
      ALLOCATE(EQUATIONS_SET,STAT=ERR)
      IF(ERR/=0) CALL FlagError("Could not allocate equations set.",err,error,*999)
      EQUATIONS_SET%USER_NUMBER=0
      EQUATIONS_SET%GLOBAL_NUMBER=0
      EQUATIONS_SET%EQUATIONS_SET_FINISHED=.FALSE.
      NULLIFY(EQUATIONS_SET%EQUATIONS_SETS)
      EQUATIONS_SET%label=""
      NULLIFY(EQUATIONS_SET%REGION)
      EQUATIONS_SET%currentTime=0.0_DP
      EQUATIONS_SET%deltaTime=0.0_DP
      EQUATIONS_SET%outputType=EQUATIONS_SET_NO_OUTPUT
      EQUATIONS_SET%SOLUTION_METHOD=0
      CALL EQUATIONS_SET_GEOMETRY_INITIALISE(EQUATIONS_SET,err,error,*999)
      CALL EQUATIONS_SET_DEPENDENT_INITIALISE(EQUATIONS_SET,err,error,*999)
      CALL EquationsSet_EquationsSetFieldInitialise(EQUATIONS_SET,err,error,*999)
      NULLIFY(EQUATIONS_SET%INDEPENDENT)
      NULLIFY(EQUATIONS_SET%MATERIALS)
      NULLIFY(EQUATIONS_SET%SOURCE)
      NULLIFY(EQUATIONS_SET%ANALYTIC)
      NULLIFY(EQUATIONS_SET%derived)
      NULLIFY(EQUATIONS_SET%EQUATIONS)
      NULLIFY(EQUATIONS_SET%BOUNDARY_CONDITIONS)
    ENDIF
       
    EXITS("EQUATIONS_SET_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_FINALISE(EQUATIONS_SET,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finalise the geometry for an equations set
  SUBROUTINE EQUATIONS_SET_GEOMETRY_FINALISE(EQUATIONS_SET_GEOMETRY,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_GEOMETRY_TYPE) :: EQUATIONS_SET_GEOMETRY !<A pointer to the equations set geometry to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_GEOMETRY_FINALISE",err,error,*999)
    
    NULLIFY(EQUATIONS_SET_GEOMETRY%GEOMETRIC_FIELD)
    NULLIFY(EQUATIONS_SET_GEOMETRY%FIBRE_FIELD)
       
    EXITS("EQUATIONS_SET_GEOMETRY_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_GEOMETRY_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_GEOMETRY_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the geometry for an equation set
  SUBROUTINE EQUATIONS_SET_GEOMETRY_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the geometry for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_SET_GEOMETRY_INITIALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%GEOMETRY%EQUATIONS_SET=>EQUATIONS_SET
      NULLIFY(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD)
      NULLIFY(EQUATIONS_SET%GEOMETRY%FIBRE_FIELD)
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_GEOMETRY_INITIALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_GEOMETRY_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_GEOMETRY_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of materials for an equations set. \see OpenCMISS::cmfe_EquationsSet_MaterialsCreateFinish
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_FINISH(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the materials field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD

    ENTERS("EQUATIONS_SET_MATERIALS_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        IF(EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED) THEN
          CALL FlagError("Equations set materials has already been finished.",err,error,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_MATERIALS_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          IF(ASSOCIATED(MATERIALS_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=MATERIALS_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>MATERIALS_FIELD
            !Finish equations set specific startup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          ELSE
            CALL FlagError("Equations set materials materials field is not associated.",err,error,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finish materials creation
          EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FlagError("The equations set materials is not associated",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_MATERIALS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_MATERIALS_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of materials for a problem. \see OpenCMISS::cmfe_EquationsSet_MaterialsCreateStart
  SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_START(EQUATIONS_SET,MATERIALS_FIELD_USER_NUMBER,MATERIALS_FIELD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of the materials field for
    INTEGER(INTG), INTENT(IN) :: MATERIALS_FIELD_USER_NUMBER !<The user specified materials field number
    TYPE(FIELD_TYPE), POINTER :: MATERIALS_FIELD !<If associated on entry, a pointer to the user created materials field which has the same user number as the specified materials field user number. If not associated on entry, on exit, a pointer to the created materials field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,MATERIALS_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EQUATIONS_SET_MATERIALS_CREATE_START",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL FlagError("The equations set materials is already associated",err,error,*998)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(MATERIALS_FIELD)) THEN
            !Check the materials field has been finished
            IF(MATERIALS_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(MATERIALS_FIELD_USER_NUMBER/=MATERIALS_FIELD%USER_NUMBER) THEN
                localError="The specified materials field user number of "// &
                  & TRIM(NumberToVString(MATERIALS_FIELD_USER_NUMBER,"*",err,error))// &
                  & " does not match the user number of the specified materials field of "// &
                  & TRIM(NumberToVString(MATERIALS_FIELD%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              MATERIALS_FIELD_REGION=>MATERIALS_FIELD%REGION
              IF(ASSOCIATED(MATERIALS_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(MATERIALS_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified materials field has been created on region number "// &
                    & TRIM(NumberToVString(MATERIALS_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Check the specified materials field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,MATERIALS_FIELD%DECOMPOSITION)) THEN
                    CALL FlagError("The specified materials field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified materials field region is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The specified materials field has not been finished.",err,error,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(MATERIALS_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
            IF(ASSOCIATED(FIELD)) THEN
              localError="The specified materials field user number of "// &
                & TRIM(NumberToVString(MATERIALS_FIELD_USER_NUMBER,"*",err,error))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
          !Initialise the equations set materials
          CALL EQUATIONS_SET_MATERIALS_INITIALISE(EQUATIONS_SET,err,error,*999)
          IF(.NOT.ASSOCIATED(MATERIALS_FIELD)) EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_MATERIALS_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=MATERIALS_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>MATERIALS_FIELD
          !Start equations set specific startup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Set pointers
          IF(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN            
            MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
          ELSE
            EQUATIONS_SET%MATERIALS%MATERIALS_FIELD=>MATERIALS_FIELD
          ENDIF
        ELSE
          CALL FlagError("Equation set region is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_MATERIALS_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_MATERIALS_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the materials for an equations set. \see OpenCMISS::cmfe_EquationsSet_MaterialsDestroy
  SUBROUTINE EQUATIONS_SET_MATERIALS_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the materials for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_MATERIALS_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,err,error,*999)
      ELSE
        CALL FlagError("Equations set materials is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_MATERIALS_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_MATERIALS_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the materials for an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET_MATERIALS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_SET_MATERIALS !<A pointer to the equations set materials to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_MATERIALS_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET_MATERIALS)) THEN
      DEALLOCATE(EQUATIONS_SET_MATERIALS)
    ENDIF
       
    EXITS("EQUATIONS_SET_MATERIALS_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_MATERIALS_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the materials for an equations set.
  SUBROUTINE EQUATIONS_SET_MATERIALS_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the materials for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EQUATIONS_SET_MATERIALS_INITIALISE",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%MATERIALS)) THEN
        CALL FlagError("Materials is already associated for these equations sets.",err,error,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%MATERIALS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate equations set materials.",err,error,*999)
        EQUATIONS_SET%MATERIALS%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%MATERIALS%MATERIALS_FINISHED=.FALSE.
        EQUATIONS_SET%MATERIALS%MATERIALS_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%MATERIALS%MATERIALS_FIELD)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_MATERIALS_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_MATERIALS_FINALISE(EQUATIONS_SET%MATERIALS,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_MATERIALS_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_MATERIALS_INITIALISE

  !
  !
  !================================================================================================================================
  !

  !>Finish the creation of a dependent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_DependentCreateFinish
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_FINISH(EQUATIONS_SET,err,error,*)
    
    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: dependentField

    ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FlagError("Equations set dependent has already been finished",err,error,*999)
      ELSE
        !Initialise the setup
        CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=dependentField%USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>dependentField
          !Finish equations set specific setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        ELSE
          CALL FlagError("Equations set dependent dependent field is not associated.",err,error,*999)
        ENDIF
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
        !Finish the equations set creation
        EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED=.TRUE.
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_DEPENDENT_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of dependent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_DependentCreateStart
  SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_START(EQUATIONS_SET,DEPENDENT_FIELD_USER_NUMBER,dependentField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a dependent field on
    INTEGER(INTG), INTENT(IN) :: DEPENDENT_FIELD_USER_NUMBER !<The user specified dependent field number
    TYPE(FIELD_TYPE), POINTER :: dependentField !<If associated on entry, a pointer to the user created dependent field which has the same user number as the specified dependent field user number. If not associated on entry, on exit, a pointer to the created dependent field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,DEPENDENT_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError
    
    ENTERS("EQUATIONS_SET_DEPENDENT_CREATE_START",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
        CALL FlagError("The equations set dependent has been finished.",err,error,*999)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(dependentField)) THEN
            !Check the dependent field has been finished
            IF(dependentField%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(DEPENDENT_FIELD_USER_NUMBER/=dependentField%USER_NUMBER) THEN
                localError="The specified dependent field user number of "// &
                  & TRIM(NumberToVString(DEPENDENT_FIELD_USER_NUMBER,"*",err,error))// &
                  & " does not match the user number of the specified dependent field of "// &
                  & TRIM(NumberToVString(dependentField%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              DEPENDENT_FIELD_REGION=>dependentField%REGION
              IF(ASSOCIATED(DEPENDENT_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(DEPENDENT_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified dependent field has been created on region number "// &
                    & TRIM(NumberToVString(DEPENDENT_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Check the specified dependent field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,dependentField%DECOMPOSITION)) THEN
                    CALL FlagError("The specified dependent field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified dependent field region is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The specified dependent field has not been finished.",err,error,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(DEPENDENT_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
            IF(ASSOCIATED(FIELD)) THEN
              localError="The specified dependent field user number of "// &
                & TRIM(NumberToVString(DEPENDENT_FIELD_USER_NUMBER,"*",err,error))// &
                & " has already been used to create a field on region number "// &
                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
            EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.TRUE.
          ENDIF
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_DEPENDENT_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=DEPENDENT_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>dependentField
          !Start the equations set specfic solution setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Set pointers
          IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
            dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
          ELSE
            EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD=>dependentField
          ENDIF
        ELSE
          CALL FlagError("Equation set region is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations_set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_DEPENDENT_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_DEPENDENT_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_CREATE_START

  !
  !================================================================================================================================
  !
  
  !>Destroy the dependent variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_DependentDestroy
  SUBROUTINE EQUATIONS_SET_DEPENDENT_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<The pointer to the equations set to destroy
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_DEPENDENT_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      CALL EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET%DEPENDENT,err,error,*999)
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
    
    EXITS("EQUATIONS_SET_DEPENDENT_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_DEPENDENT_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_DESTROY
  
  !
  !================================================================================================================================
  !

  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_FINALISE(EQUATIONS_SET_DEPENDENT,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_DEPENDENT_TYPE) :: EQUATIONS_SET_DEPENDENT !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_DEPENDENT_FINALISE",err,error,*999)

    NULLIFY(EQUATIONS_SET_DEPENDENT%EQUATIONS_SET)
    EQUATIONS_SET_DEPENDENT%DEPENDENT_FINISHED=.FALSE.
    EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
    NULLIFY(EQUATIONS_SET_DEPENDENT%DEPENDENT_FIELD)
    
    EXITS("EQUATIONS_SET_DEPENDENT_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_DEPENDENT_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialises the dependent variables for a equations set.
  SUBROUTINE EQUATIONS_SET_DEPENDENT_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EQUATIONS_SET_DEPENDENT_INITIALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%DEPENDENT%EQUATIONS_SET=>EQUATIONS_SET
      EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED=.FALSE.
      EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED=.FALSE.
      NULLIFY(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD)
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_DEPENDENT_INITIALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_DEPENDENT_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_DEPENDENT_INITIALISE

  !
  !================================================================================================================================
  !

  !>Finish the creation of a derived variables field for an equations set. \see OpenCMISS::cmfe_EquationsSet_DerivedCreateFinish
  SUBROUTINE EquationsSet_DerivedCreateFinish(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to finish the derived variable creation for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: equationsSetSetupInfo
    TYPE(FIELD_TYPE), POINTER :: derivedField

    ENTERS("EquationsSet_DerivedCreateFinish",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(ASSOCIATED(equationsSet%derived)) THEN
        IF(equationsSet%derived%derivedFinished) THEN
          CALL FlagError("Equations set derived field information has already been finished",err,error,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(equationsSetSetupInfo,err,error,*999)
          equationsSetSetupInfo%SETUP_TYPE=EQUATIONS_SET_SETUP_DERIVED_TYPE
          equationsSetSetupInfo%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          derivedField=>equationsSet%derived%derivedField
          IF(ASSOCIATED(derivedField)) THEN
            equationsSetSetupInfo%FIELD_USER_NUMBER=derivedField%USER_NUMBER
            equationsSetSetupInfo%field=>derivedField
            !Finish equations set specific setup
            CALL EQUATIONS_SET_SETUP(equationsSet,equationsSetSetupInfo,err,error,*999)
          ELSE
            CALL FlagError("Equations set derived field is not associated.",err,error,*999)
          END IF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(equationsSetSetupInfo,err,error,*999)
          !Finish the equations set derived creation
          equationsSet%derived%derivedFinished=.TRUE.
        END IF
      ELSE
        CALL FlagError("Equations set derived is not associated",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    EXITS("EquationsSet_DerivedCreateFinish")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedCreateFinish",err,error)
    RETURN 1
  END SUBROUTINE EquationsSet_DerivedCreateFinish

  !
  !================================================================================================================================
  !

  !>Start the creation of derived variables field for an equations set. \see OpenCMISS::cmfe_EquationsSet_DerivedCreateStart
  SUBROUTINE EquationsSet_DerivedCreateStart(equationsSet,derivedFieldUserNumber,derivedField,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to start the creation of a derived field on
    INTEGER(INTG), INTENT(IN) :: derivedFieldUserNumber !<The user specified derived field number
    TYPE(FIELD_TYPE), POINTER :: derivedField !<If associated on entry, a pointer to the user created derived field which has the same user number as the specified derived field user number. If not associated on entry, on exit, a pointer to the created derived field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: equationsSetSetupInfo
    TYPE(FIELD_TYPE), POINTER :: field,geometricField
    TYPE(REGION_TYPE), POINTER :: region,derivedFieldRegion
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EquationsSet_DerivedCreateStart",err,error,*998)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(ASSOCIATED(equationsSet%derived)) THEN
        CALL FlagError("Equations set derived is already associated.",err,error,*998)
      ELSE
        region=>equationsSet%REGION
        IF(ASSOCIATED(region)) THEN
          IF(ASSOCIATED(derivedField)) THEN
            !Check the derived field has been finished
            IF(derivedField%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(derivedFieldUserNumber/=derivedField%USER_NUMBER) THEN
                localError="The specified derived field user number of "// &
                  & TRIM(NumberToVString(derivedFieldUserNumber,"*",err,error))// &
                  & " does not match the user number of the specified derived field of "// &
                  & TRIM(NumberToVString(derivedField%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              END IF
              derivedFieldRegion=>derivedField%REGION
              IF(ASSOCIATED(derivedFieldRegion)) THEN
                !Check the field is defined on the same region as the equations set
                IF(derivedFieldRegion%USER_NUMBER/=region%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified derived field has been created on region number "// &
                    & TRIM(NumberToVString(derivedFieldRegion%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                END IF
                !Check the specified derived field has the same decomposition as the geometric field
                geometricField=>equationsSet%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(geometricField)) THEN
                  IF(.NOT.ASSOCIATED(geometricField%DECOMPOSITION,derivedField%DECOMPOSITION)) THEN
                    CALL FlagError("The specified derived field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  END IF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                END IF
              ELSE
                CALL FlagError("The specified derived field region is not associated.",err,error,*999)
              END IF
            ELSE
              CALL FlagError("The specified derived field has not been finished.",err,error,*999)
            END IF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(field)
            CALL FIELD_USER_NUMBER_FIND(derivedFieldUserNumber,region,field,err,error,*999)
            IF(ASSOCIATED(field)) THEN
              localError="The specified derived field user number of "// &
                & TRIM(NumberToVString(derivedFieldUserNumber,"*",err,error))// &
                & " has already been used to create a field on region number "// &
                & TRIM(NumberToVString(region%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            END IF
            equationsSet%derived%derivedFieldAutoCreated=.TRUE.
          END IF
          CALL EquationsSet_DerivedInitialise(equationsSet,err,error,*999)
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(equationsSetSetupInfo,err,error,*999)
          equationsSetSetupInfo%SETUP_TYPE=EQUATIONS_SET_SETUP_DERIVED_TYPE
          equationsSetSetupInfo%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          equationsSetSetupInfo%FIELD_USER_NUMBER=derivedFieldUserNumber
          equationsSetSetupInfo%FIELD=>derivedField
          !Start the equations set specfic solution setup
          CALL EQUATIONS_SET_SETUP(equationsSet,equationsSetSetupInfo,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(equationsSetSetupInfo,err,error,*999)
          !Set pointers
          IF(.NOT.equationsSet%derived%derivedFieldAutoCreated) THEN
            equationsSet%derived%derivedField=>derivedField
          END IF
        ELSE
          CALL FlagError("Equation set region is not associated.",err,error,*999)
        END IF
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    END IF

    EXITS("EquationsSet_DerivedCreateStart")
    RETURN
999 CALL EquationsSet_DerivedFinalise(equationsSet%derived,dummyErr,dummyError,*998)
998 ERRORSEXITS("EquationsSet_DerivedCreateStart",err,error)
    RETURN 1
  END SUBROUTINE EquationsSet_DerivedCreateStart

  !
  !================================================================================================================================
  !

  !>Destroy the derived variables for an equations set. \see OpenCMISS::cmfe_EquationsSet_DerivedDestroy
  SUBROUTINE EquationsSet_DerivedDestroy(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<The pointer to the equations set to destroy the derived fields for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_DerivedDestroy",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      CALL EquationsSet_DerivedFinalise(equationsSet%derived,err,error,*999)
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

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
    END IF

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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to initialise the derived field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_DerivedInitialise",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(ASSOCIATED(equationsSet%derived)) THEN
        CALL FlagError("Derived information is already associated for this equations set.",err,error,*998)
      ELSE
        ALLOCATE(equationsSet%derived,stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set derived information.",err,error,*998)
        ALLOCATE(equationsSet%derived%variableTypes(EQUATIONS_SET_NUMBER_OF_TENSOR_TYPES),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set derived variable types.",err,error,*999)
        equationsSet%derived%variableTypes=0
        equationsSet%derived%numberOfVariables=0
        equationsSet%derived%equationsSet=>equationsSet
        equationsSet%derived%derivedFinished=.FALSE.
        equationsSet%derived%derivedFieldAutoCreated=.FALSE.
        NULLIFY(equationsSet%derived%derivedField)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("EquationsSet_DerivedInitialise")
    RETURN
999 CALL EquationsSet_DerivedFinalise(equationsSet%derived,err,error,*999)
998 ERRORSEXITS("EquationsSet_DerivedInitialise",err,error)
    RETURN 1
  END SUBROUTINE EquationsSet_DerivedInitialise

  !
  !================================================================================================================================
  !

  !>Finalises the dependent variables for an equation set and deallocates all memory.
  SUBROUTINE EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE(EQUATIONS_SET_FIELD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_EQUATIONS_SET_FIELD_TYPE) :: EQUATIONS_SET_FIELD !<The pointer to the equations set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE",err,error,*999)

    NULLIFY(EQUATIONS_SET_FIELD%EQUATIONS_SET)
    EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FINISHED=.FALSE.
    EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED=.FALSE.
    NULLIFY(EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD)
    
    EXITS("EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_SET_FIELD_FINALISE
  
  !
  !================================================================================================================================
  !
  !>Initialises the equations set field for a equations set.
  SUBROUTINE EquationsSet_EquationsSetFieldInitialise(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the dependent field for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_EquationsSetFieldInitialise",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET=>EQUATIONS_SET
      EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FINISHED=.FALSE.
      EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_AUTO_CREATED=.TRUE.
      NULLIFY(EQUATIONS_SET%EQUATIONS_SET_FIELD%EQUATIONS_SET_FIELD_FIELD)
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_EquationsSetFieldInitialise")
    RETURN
999 ERRORSEXITS("EquationsSet_EquationsSetFieldInitialise",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_EquationsSetFieldInitialise

  !
  !================================================================================================================================
  !



  !>Sets up the specifices for an equation set.
  SUBROUTINE EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the setup on
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_SETUP",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<1) THEN
        CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(1))
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
        CALL FLUID_MECHANICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
        CALL CLASSICAL_FIELD_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
        IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
          CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation class.", &
            & err,error,*999)
        END IF
        IF(EQUATIONS_SET%SPECIFICATION(2) == EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
          CALL MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        ELSE
          CALL BIOELECTRIC_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        END IF
      CASE(EQUATIONS_SET_FITTING_CLASS)
        CALL Fitting_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      CASE(EQUATIONS_SET_MODAL_CLASS)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
        CALL MULTI_PHYSICS_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      CASE DEFAULT
        localError="The first equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(1),"*",err,error))//" is not valid."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

 !>Finish the creation of equations for the equations set. \see OpenCMISS::cmfe_EquationsSet_EquationsCreateFinish
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_FINISH(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to finish the creation of the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    
    ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Initialise the setup
      CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
      EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
      EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
      !Finish the equations specific solution setup.
      CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
      !Finalise the setup
      CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_EQUATIONS_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of equations for the equation set. \see CMISSEquationsSetEquationsCreateStart
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
  SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_START(EQUATIONS_SET,equations,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to create equations for
    TYPE(EquationsType), POINTER :: equations !<On exit, a pointer to the created equations. Must not be associated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO

    ENTERS("EQUATIONS_SET_EQUATIONS_CREATE_START",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS)) THEN
        CALL FlagError("Equations is already associated.",err,error,*999)
      ELSE
        !Initialise the setup
        CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
        EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_EQUATIONS_TYPE
        EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
        !Start the equations set specific solution setup
        CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
        !Finalise the setup
        CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
        !Return the pointer
        equations=>EQUATIONS_SET%EQUATIONS
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_EQUATIONS_CREATE_START")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_EQUATIONS_CREATE_START",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the equations for an equations set. \see OpenCMISS::cmfe_EquationsSet_EquationsDestroy
  SUBROUTINE EQUATIONS_SET_EQUATIONS_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the equations for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_EQUATIONS_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%EQUATIONS)) THEN
        CALL EQUATIONS_FINALISE(EQUATIONS_SET%EQUATIONS,err,error,*999)
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_EQUATIONS_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_EQUATIONS_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_EQUATIONS_DESTROY

  !
  !================================================================================================================================
  !

  !>Evaluates the Jacobian for a nonlinear equations set.
  SUBROUTINE EquationsSet_JacobianEvaluate(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
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
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_QUASISTATIC)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
            & TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_QUASISTATIC)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        ! sebk 15/09/09
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
  
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_JacobianElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_JacobianElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
  
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_JacobianElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementJacobianEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_JacobianElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: residualVariableIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingNonlinearType), POINTER :: nonlinearMapping
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_PARAMETER_SET_TYPE), POINTER :: residualParameterSet
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: residualVariable
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
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
            & " is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_FIRST_ORDER_DYNAMIC,EQUATIONS_SECOND_ORDER_DYNAMIC)
        SELECT CASE(equationsSet%SOLUTION_METHOD)
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
          localError="The equations set solution method  of "//TRIM(NumberToVString(equationsSet%SOLUTION_METHOD,"*",err,error))// &
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
    DO residualVariableIdx=1,nonlinearMapping%numberOfResidualVariables
      residualVariable=>nonlinearMapping%residualVariables(residualVariableIdx)%ptr
      IF(.NOT.ASSOCIATED(residualVariable)) THEN
        localError="Nonlinear mapping residual variable for residual variable index "// &
          & TRIM(NumberToVString(residualVariableIdx,"*",err,error))//" is not associated."
        CALL FlagError(localError,err,error,*999)
      ENDIF
      residualParameterSet=>residualVariable%PARAMETER_SETS%SET_TYPE(FIELD_RESIDUAL_SET_TYPE)%ptr
      IF(ASSOCIATED(residualParameterSet)) THEN
        !Residual parameter set exists. Copy the residual vector to the residuals parameter set.
        CALL DistributedVector_Copy(nonlinearMatrices%residual,residualParameterSet%parameters,1.0_DP,err,error,*999)
      ENDIF
    ENDDO !residualVariableIdx
       
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
 
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: elementIdx,element,numberOfTimes
    REAL(SP) :: elementUserElapsed,elementSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: elementsMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
 
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
!!Do we need to transfer parameter sets???
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Assemble the elements
    !Allocate the element matrices 
    CALL EquationsMatrices_ElementInitialise(vectorMatrices,err,error,*999)
    elementsMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr% &
      & mappings%elements
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
    DO elementIdx=elementsMapping%INTERNAL_START,elementsMapping%INTERNAL_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    DO elementIdx=elementsMapping%BOUNDARY_START,elementsMapping%GHOST_FINISH
      element=elementsMapping%DOMAIN_LIST(elementIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_ElementCalculate(vectorMatrices,element,err,error,*999)
      CALL EquationsSet_FiniteElementResidualEvaluate(equationsSet,element,err,error,*999)
      CALL EquationsMatrices_ElementAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_ElementFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
  SUBROUTINE EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(OUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup to be finalised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_SETUP_FINALISE",err,error,*999)

    EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=0
    EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=0
    EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=0
    NULLIFY(EQUATIONS_SET_SETUP_INFO%FIELD)
    EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=0
    
    EXITS("EQUATIONS_SET_SETUP_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SETUP_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP_FINALISE
  
  !
  !================================================================================================================================
  !

  !>Initialise the equations set setup.
  SUBROUTINE EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(OUT) :: EQUATIONS_SET_SETUP_INFO !<The equations set setup to be initialised
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_SETUP_INITIALISE",err,error,*999)

    EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=0
    EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=0
    EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=0
    NULLIFY(EQUATIONS_SET_SETUP_INFO%FIELD)
    EQUATIONS_SET_SETUP_INFO%ANALYTIC_FUNCTION_TYPE=0
    
    EXITS("EQUATIONS_SET_SETUP_INITIALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SETUP_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SETUP_INITIALISE
  
  !
  !================================================================================================================================
  !

  !>Gets the output type for an equations set.
  SUBROUTINE EquationsSet_OutputTypeGet(equationsSet,outputType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the output type for
    INTEGER(INTG), INTENT(OUT) :: outputType !<On exit, the output type of the equations set \see EquationsSetConstants_OutputTypes,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
 
    ENTERS("EquationsSet_OutputTypeGet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
    
    outputType=equationsSet%outputType
       
    EXITS("EquationsSet_OutputTypeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_OutputTypeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_OutputTypeGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the output type for an equations set.
  SUBROUTINE EquationsSet_OutputTypeSet(equationsSet,outputType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the output type for
    INTEGER(INTG), INTENT(IN) :: outputType !<The output type to set \see EquationsSetConstants_OutputTypes,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
 
    ENTERS("EquationsSet_OutputTypeSet",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has already been finished.",err,error,*999)

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

  !>Sets/changes the solution method for an equations set. \see OpenCMISS::cmfe_EquationsSet_SolutionMethodSet
  SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_SET(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The equations set solution method to set \see EquationsSetConstants_SolutionMethods,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_SOLUTION_METHOD_SET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        CALL FlagError("Equations set has already been finished.",err,error,*999)
      ELSE
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<1) THEN
          CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
        END IF
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(1))
        CASE(EQUATIONS_SET_ELASTICITY_CLASS)
          CALL Elasticity_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
        CASE(EQUATIONS_SET_FLUID_MECHANICS_CLASS)
          CALL FluidMechanics_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
        CASE(EQUATIONS_SET_ELECTROMAGNETICS_CLASS)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_CLASSICAL_FIELD_CLASS)
          CALL ClassicalField_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
        CASE(EQUATIONS_SET_BIOELECTRICS_CLASS)
          IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
            CALL FlagError("Equations set specification must have at least two entries for a bioelectrics equation set.", &
              & err,error,*999)
          END IF
          IF(EQUATIONS_SET%SPECIFICATION(2) == EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) THEN
            CALL Monodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
          ELSE
            CALL Bioelectric_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
          END IF
        CASE(EQUATIONS_SET_MODAL_CLASS)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_MULTI_PHYSICS_CLASS)
          CALL MultiPhysics_EquationsSetSolnMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*999)
        CASE DEFAULT
          localError="The first equations set specification of "// &
            & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(1),"*",err,error))//" is invalid."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EQUATIONS_SET_SOLUTION_METHOD_SET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SOLUTION_METHOD_SET",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_SET
  
  !
  !================================================================================================================================
  !

  !>Returns the solution method for an equations set. \see OpenCMISS::cmfe_EquationsSet_SolutionMethodGet
  SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_GET(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to get the solution method for
    INTEGER(INTG), INTENT(OUT) :: SOLUTION_METHOD !<On return, the equations set solution method \see EquationsSetConstants_SolutionMethods,EquationsSetConstants
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_SOLUTION_METHOD_GET",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(EQUATIONS_SET%EQUATIONS_SET_FINISHED) THEN
        SOLUTION_METHOD=EQUATIONS_SET%SOLUTION_METHOD
      ELSE
        CALL FlagError("Equations set has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("EQUATIONS_SET_SOLUTION_METHOD_GET")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SOLUTION_METHOD_GET",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOLUTION_METHOD_GET
  
  !
  !================================================================================================================================
  !

  !>Finish the creation of a source for an equation set. \see OpenCMISS::cmfe_EquationsSet_SourceCreateFinish
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_FINISH(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a souce for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD

    ENTERS("EQUATIONS_SET_SOURCE_CREATE_FINISH",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        IF(EQUATIONS_SET%SOURCE%SOURCE_FINISHED) THEN
          CALL FlagError("Equations set source has already been finished.",err,error,*999)
        ELSE
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_SOURCE_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_FINISH_ACTION
          SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=SOURCE_FIELD%USER_NUMBER
            EQUATIONS_SET_SETUP_INFO%FIELD=>SOURCE_FIELD
            !Finish the equation set specific source setup
            CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          ELSE
            CALL FlagError("Equations set source source field is not associated.",err,error,*999)
          ENDIF
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finish the source creation
          EQUATIONS_SET%SOURCE%SOURCE_FINISHED=.TRUE.
        ENDIF
      ELSE
        CALL FlagError("The equations set source is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_SOURCE_CREATE_FINISH")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SOURCE_CREATE_FINISH",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_FINISH

  !
  !================================================================================================================================
  !

  !>Start the creation of a source for an equations set. \see OpenCMISS::cmfe_EquationsSet_SourceCreateStart
  SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_START(EQUATIONS_SET,SOURCE_FIELD_USER_NUMBER,SOURCE_FIELD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to start the creation of a source for
    INTEGER(INTG), INTENT(IN) :: SOURCE_FIELD_USER_NUMBER !<The user specified source field number
    TYPE(FIELD_TYPE), POINTER :: SOURCE_FIELD !<If associated on entry, a pointer to the user created source field which has the same user number as the specified source field user number. If not associated on entry, on exit, a pointer to the created source field for the equations set.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(EQUATIONS_SET_SETUP_TYPE) :: EQUATIONS_SET_SETUP_INFO
    TYPE(FIELD_TYPE), POINTER :: FIELD,GEOMETRIC_FIELD
    TYPE(REGION_TYPE), POINTER :: REGION,SOURCE_FIELD_REGION
    TYPE(VARYING_STRING) :: dummyError,localError

    ENTERS("EQUATIONS_SET_SOURCE_CREATE_START",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL FlagError("The equations set source is already associated.",err,error,*998)
      ELSE
        REGION=>EQUATIONS_SET%REGION
        IF(ASSOCIATED(REGION)) THEN
          IF(ASSOCIATED(SOURCE_FIELD)) THEN
            !Check the source field has been finished
            IF(SOURCE_FIELD%FIELD_FINISHED) THEN
              !Check the user numbers match
              IF(SOURCE_FIELD_USER_NUMBER/=SOURCE_FIELD%USER_NUMBER) THEN
                localError="The specified source field user number of "// &
                  & TRIM(NumberToVString(SOURCE_FIELD_USER_NUMBER,"*",err,error))// &
                  & " does not match the user number of the specified source field of "// &
                  & TRIM(NumberToVString(SOURCE_FIELD%USER_NUMBER,"*",err,error))//"."
                CALL FlagError(localError,err,error,*999)
              ENDIF
              SOURCE_FIELD_REGION=>SOURCE_FIELD%REGION
              IF(ASSOCIATED(SOURCE_FIELD_REGION)) THEN                
                !Check the field is defined on the same region as the equations set
                IF(SOURCE_FIELD_REGION%USER_NUMBER/=REGION%USER_NUMBER) THEN
                  localError="Invalid region setup. The specified source field has been created on region number "// &
                    & TRIM(NumberToVString(SOURCE_FIELD_REGION%USER_NUMBER,"*",err,error))// &
                    & " and the specified equations set has been created on region number "// &
                    & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
                  CALL FlagError(localError,err,error,*999)
                ENDIF
                !Check the specified source field has the same decomposition as the geometric field
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  IF(.NOT.ASSOCIATED(GEOMETRIC_FIELD%DECOMPOSITION,SOURCE_FIELD%DECOMPOSITION)) THEN
                    CALL FlagError("The specified source field does not have the same decomposition as the geometric "// &
                      & "field for the specified equations set.",err,error,*999)
                  ENDIF
                ELSE
                  CALL FlagError("The geometric field is not associated for the specified equations set.",err,error,*999)
                ENDIF
              ELSE
                CALL FlagError("The specified source field region is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("The specified source field has not been finished.",err,error,*999)
            ENDIF
          ELSE
            !Check the user number has not already been used for a field in this region.
            NULLIFY(FIELD)
            CALL FIELD_USER_NUMBER_FIND(SOURCE_FIELD_USER_NUMBER,REGION,FIELD,err,error,*999)
            IF(ASSOCIATED(FIELD)) THEN
              localError="The specified source field user number of "// &
                & TRIM(NumberToVString(SOURCE_FIELD_USER_NUMBER,"*",err,error))// &
                & "has already been used to create a field on region number "// &
                & TRIM(NumberToVString(REGION%USER_NUMBER,"*",err,error))//"."
              CALL FlagError(localError,err,error,*999)
            ENDIF
          ENDIF
          !Initialise the equations set source
          CALL EQUATIONS_SET_SOURCE_INITIALISE(EQUATIONS_SET,err,error,*999)
          IF(.NOT.ASSOCIATED(SOURCE_FIELD)) EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED=.TRUE.
          !Initialise the setup
          CALL EQUATIONS_SET_SETUP_INITIALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          EQUATIONS_SET_SETUP_INFO%SETUP_TYPE=EQUATIONS_SET_SETUP_SOURCE_TYPE
          EQUATIONS_SET_SETUP_INFO%ACTION_TYPE=EQUATIONS_SET_SETUP_START_ACTION
          EQUATIONS_SET_SETUP_INFO%FIELD_USER_NUMBER=SOURCE_FIELD_USER_NUMBER
          EQUATIONS_SET_SETUP_INFO%FIELD=>SOURCE_FIELD
          !Start the equation set specific source setup
          CALL EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Finalise the setup
          CALL EQUATIONS_SET_SETUP_FINALISE(EQUATIONS_SET_SETUP_INFO,err,error,*999)
          !Set pointers
          IF(EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED) THEN            
            SOURCE_FIELD=>EQUATIONS_SET%SOURCE%SOURCE_FIELD
          ELSE
            EQUATIONS_SET%SOURCE%SOURCE_FIELD=>SOURCE_FIELD
          ENDIF
        ELSE
          CALL FlagError("Equation set region is not associated.",err,error,*999)
        ENDIF
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_SOURCE_CREATE_START")
    RETURN
999 CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_SOURCE_CREATE_START",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_CREATE_START

  !
  !================================================================================================================================
  !

  !>Destroy the source for an equations set. \see OpenCMISS::cmfe_EquationsSet_SourceDestroy
  SUBROUTINE EQUATIONS_SET_SOURCE_DESTROY(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to destroy the source for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_SOURCE_DESTROY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,err,error,*999)
      ELSE
        CALL FlagError("Equations set source is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    ENDIF
       
    EXITS("EQUATIONS_SET_SOURCE_DESTROY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SOURCE_DESTROY",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_DESTROY

  !
  !================================================================================================================================
  !

  !>Finalise the source for a equations set and deallocate all memory.
  SUBROUTINE EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET_SOURCE,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_SOURCE_TYPE), POINTER :: EQUATIONS_SET_SOURCE !<A pointer to the equations set source to finalise
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SET_SOURCE_FINALISE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET_SOURCE)) THEN
       DEALLOCATE(EQUATIONS_SET_SOURCE)
    ENDIF
       
    EXITS("EQUATIONS_SET_SOURCE_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_SOURCE_FINALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SET_SOURCE_FINALISE

  !
  !================================================================================================================================
  !

  !>Initialises the source for an equations set.
  SUBROUTINE EQUATIONS_SET_SOURCE_INITIALISE(EQUATIONS_SET,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to initialise the source field for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: dummyErr
    TYPE(VARYING_STRING) :: dummyError
    
    ENTERS("EQUATIONS_SET_SOURCE_INITIALISE",err,error,*998)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%SOURCE)) THEN
        CALL FlagError("Source is already associated for this equations set.",err,error,*998)
      ELSE
        ALLOCATE(EQUATIONS_SET%SOURCE,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate equations set source.",err,error,*999)
        EQUATIONS_SET%SOURCE%EQUATIONS_SET=>EQUATIONS_SET
        EQUATIONS_SET%SOURCE%SOURCE_FINISHED=.FALSE.
        EQUATIONS_SET%SOURCE%SOURCE_FIELD_AUTO_CREATED=.FALSE.
        NULLIFY(EQUATIONS_SET%SOURCE%SOURCE_FIELD)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*998)
    ENDIF
       
    EXITS("EQUATIONS_SET_SOURCE_INITIALISE")
    RETURN
999 CALL EQUATIONS_SET_SOURCE_FINALISE(EQUATIONS_SET%SOURCE,dummyErr,dummyError,*998)
998 ERRORSEXITS("EQUATIONS_SET_SOURCE_INITIALISE",err,error)
    RETURN 1
    
  END SUBROUTINE EQUATIONS_SET_SOURCE_INITIALISE

  !
  !================================================================================================================================
  !

  !>Returns the equations set specification i.e., equations set class, type and subtype for an equations set. \see OpenCMISS::cmfe_EquationsSet_SpecificationGet
  SUBROUTINE EquationsSet_SpecificationGet(equationsSet,equationsSetSpecification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the specification for
    INTEGER(INTG), INTENT(INOUT) :: equationsSetSpecification(:) !<On return, The equations set specifcation array. Must be allocated on entry.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: specificationLength,specificationIdx
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_SpecificationGet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%equations_set_finished) THEN
        specificationLength=0
        DO specificationIdx=1,SIZE(equationsSet%specification,1)
          IF(equationsSet%specification(specificationIdx)>0) THEN
            specificationLength=specificationIdx
          END IF
        END DO
        IF(SIZE(equationsSetSpecification,1)>=specificationLength) THEN
          equationsSetSpecification(1:specificationLength)=equationsSet%specification(1:specificationLength)
        ELSE
          localError="The equations set specification array size is "//TRIM(NumberToVstring(specificationLength,"*",err,error))// &
            & " and it needs to be >= "//TRIM(NumberToVstring(SIZE(equationsSetSpecification,1),"*",err,error))//"."
          CALL FlagError(localError,err,error,*999)
        END IF
      ELSE
        CALL FlagError("Equations set has not been finished.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("EquationsSet_SpecificationGet")
    RETURN
999 ERRORS("EquationsSet_SpecificationGet",err,error)
    EXITS("EquationsSet_SpecificationGet")
    RETURN 1
    
  END SUBROUTINE EquationsSet_SpecificationGet

  !
  !================================================================================================================================
  !

  !>Gets the size of the equations set specification array for a problem identified by a pointer. \see OpenCMISS::cmfe_EquationsSet_SpecificationSizeGet
  SUBROUTINE EquationsSet_SpecificationSizeGet(equationsSet,specificationSize,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations Set to get the specification for.
    INTEGER(INTG), INTENT(OUT) :: specificationSize !<On return, the size of the problem specifcation array.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EquationsSet_SpecificationSizeGet",err,error,*999)

    specificationSize=0
    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%equations_set_finished) THEN
        IF(.NOT.ALLOCATED(equationsSet%specification)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        END IF
        specificationSize=SIZE(equationsSet%specification,1)
      ELSE
        CALL FlagError("Equations set has not been finished.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EquationsSet_SpecificationSizeGet")
    RETURN
999 ERRORSEXITS("EquationsSet_SpecificationSizeGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_SpecificationSizeGet

  !
  !================================================================================================================================
  !

  !>Gets the current times for an equations set. \see OpenCMISS::cmfe_EquationsSet_TimesGet
  SUBROUTINE EquationsSet_TimesGet(equationsSet,currentTime,deltaTime,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to get the times for
    REAL(DP), INTENT(OUT) :: currentTime !<The current time for the equations set to get.
    REAL(DP), INTENT(OUT) :: deltaTime !<The current time incremenet for the equations set to get
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_TimesGet",err,error,*999) 

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)

    currentTime=equationsSet%currentTime
    deltaTime=equationsSet%deltaTime
      
    EXITS("EquationsSet_TimesGet")
    RETURN
999 ERRORSEXITS("EquationsSet_TimesGet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_TimesGet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the current times for an equations set. \see OpenCMISS::cmfe_EquationsSet_TimesSet
  SUBROUTINE EquationsSet_TimesSet(equationsSet,currentTime,deltaTime,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the times for
    REAL(DP), INTENT(IN) :: currentTime !<The current time for the equations set to set.
    REAL(DP), INTENT(IN) :: deltaTime !<The current time incremenet for the equations set to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_TimesSet",err,error,*999) 

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
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

  !>Calculates a derived variable value for the equations set. \see OpenCMISS::cmfe_EquationsSet_DerivedVariableCalculate
  SUBROUTINE EquationsSet_DerivedVariableCalculate(equationsSet,derivedType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate output for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to calculate. \see EquationsSetConstants_DerivedTypes.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_DerivedVariableCalculate",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
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

  !>Sets the field variable type of the derived field to be used to store a derived variable. \see OpenCMISS::cmfe_EquationsSet_DerivedVariableSet
  SUBROUTINE EquationsSet_DerivedVariableSet(equationsSet,derivedType,fieldVariableType,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to calculate a derived field for
    INTEGER(INTG), INTENT(IN) :: derivedType !<The derived value type to calculate. \see EquationsSetConstants_DerivedTypes.
    INTEGER(INTG), INTENT(IN) :: fieldVariableType !<The field variable type used to store the calculated derived value
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(FIELD_TYPE), POINTER :: derivedField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: fieldVariable
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_DerivedVariableSet",err,error,*999)

    !Check pointers and finished state
    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%EQUATIONS_SET_FINISHED) CALL FlagError("Equations set has not been finished.",err,error,*999)
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
    CALL Field_VariableGet(derivedField,fieldVariableType,fieldVariable,err,error,*999)
    
    IF(equationsSet%derived%variableTypes(derivedType)==0) &
      & equationsSet%derived%numberOfVariables=equationsSet%derived%numberOfVariables+1
    equationsSet%derived%variableTypes(derivedType)=fieldVariableType

    EXITS("EquationsSet_DerivedVariableSet")
    RETURN
999 ERRORSEXITS("EquationsSet_DerivedVariableSet",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_DerivedVariableSet
  
  !
  !================================================================================================================================
  !

  !>Sets/changes the equations set specification i.e., equations set class, type and subtype for an equations set. \see OpenCMISS::cmfe_EquationsSet_SpecificationSet
  SUBROUTINE EquationsSet_SpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification array to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_SpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(equationsSet%equations_set_finished) THEN
        CALL FlagError("Equations set has been finished.",err,error,*999)
      ELSE
        IF(SIZE(specification,1)<1) THEN
          CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
        END IF
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
      END IF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

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
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: gaussPointNumber !<The Gauss point number of the field to interpolate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("EquationsSet_TensorInterpolateGaussPoint",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%equations_set_finished) CALL FlagError("Equations set has not been finished.",err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER, INTENT(IN) :: equationsSet !<A pointer to the equations set to interpolate the tensor for.
    INTEGER(INTG), INTENT(IN) :: tensorEvaluateType !<The type of tensor to evaluate.
    INTEGER(INTG), INTENT(IN) :: userElementNumber !<The user element number of the field to interpolate.
    REAL(DP), INTENT(IN) :: xi(:) !<The element xi to interpolate the field at.
    REAL(DP), INTENT(OUT) :: values(:,:) !<On exit, the interpolated tensor values.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code.
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EquationsSet_TensorInterpolateXi",err,error,*999)

    IF(.NOT.ASSOCIATED(equationsSet)) CALL FlagError("Equations set is not associated.",err,error,*999)
    IF(.NOT.equationsSet%equations_set_finished) CALL FlagError("Equations set has not been finished.",err,error,*999)
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
  SUBROUTINE EQUATIONS_SETS_FINALISE(REGION,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to finalise the problems for.
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SETS_FINALISE",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        DO WHILE(REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS>0)
          CALL EQUATIONS_SET_DESTROY(REGION%EQUATIONS_SETS%EQUATIONS_SETS(1)%ptr,err,error,*999)
        ENDDO !problem_idx
        DEALLOCATE(REGION%EQUATIONS_SETS)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SETS_FINALISE")
    RETURN
999 ERRORSEXITS("EQUATIONS_SETS_FINALISE",err,error)
    RETURN 1   
  END SUBROUTINE EQUATIONS_SETS_FINALISE

  !
  !================================================================================================================================
  !

  !>Intialises all equations sets on a region.
  SUBROUTINE EQUATIONS_SETS_INITIALISE(REGION,err,error,*)

    !Argument variables
    TYPE(REGION_TYPE), POINTER :: REGION !<A pointer to the region to initialise the equations sets for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables

    ENTERS("EQUATIONS_SETS_INITIALISE",err,error,*999)

    IF(ASSOCIATED(REGION)) THEN
      IF(ASSOCIATED(REGION%EQUATIONS_SETS)) THEN
        CALL FlagError("Region already has associated equations sets",err,error,*998)
      ELSE
!!TODO: Inherit any equations sets from the parent region???
        ALLOCATE(REGION%EQUATIONS_SETS,STAT=ERR)
        IF(ERR/=0) CALL FlagError("Could not allocate region equations sets",err,error,*999)
        REGION%EQUATIONS_SETS%REGION=>REGION
        REGION%EQUATIONS_SETS%NUMBER_OF_EQUATIONS_SETS=0
        NULLIFY(REGION%EQUATIONS_SETS%EQUATIONS_SETS)
      ENDIF
    ELSE
      CALL FlagError("Region is not associated.",err,error,*998)
    ENDIF

    EXITS("EQUATIONS_SETS_INITIALISE")
    RETURN
999 IF(ASSOCIATED(REGION%EQUATIONS_SETS)) DEALLOCATE(REGION%EQUATIONS_SETS)
998 ERRORSEXITS("EQUATIONS_SETS_INITIALISE",err,error)
    RETURN 1
  END SUBROUTINE EQUATIONS_SETS_INITIALISE
  
  !
  !================================================================================================================================
  !

  !> Apply the boundary condition load increment to dependent field
  SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT(EQUATIONS_SET,BOUNDARY_CONDITIONS,ITERATION_NUMBER, &
    & MAXIMUM_NUMBER_OF_ITERATIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<The boundary conditions to apply the increment to
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    !Local variables
    TYPE(FIELD_TYPE), POINTER :: dependentField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: DEPENDENT_VARIABLE
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DOMAIN_MAPPING
    TYPE(BOUNDARY_CONDITIONS_VARIABLE_TYPE), POINTER :: BOUNDARY_CONDITIONS_VARIABLE
    TYPE(BOUNDARY_CONDITIONS_DIRICHLET_TYPE), POINTER :: DIRICHLET_BOUNDARY_CONDITIONS
    TYPE(BOUNDARY_CONDITIONS_PRESSURE_INCREMENTED_TYPE), POINTER :: PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
    INTEGER(INTG) :: variable_idx,variable_type,dirichlet_idx,dirichlet_dof_idx,neumann_point_dof
    INTEGER(INTG) :: condition_idx, condition_global_dof, condition_local_dof, myComputationalNodeNumber
    REAL(DP), POINTER :: FULL_LOADS(:),CURRENT_LOADS(:), PREV_LOADS(:)
    REAL(DP) :: FULL_LOAD, CURRENT_LOAD, NEW_LOAD, PREV_LOAD
    TYPE(VARYING_STRING) :: localError

    ENTERS("EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT",err,error,*999)

    NULLIFY(dependentField)
    NULLIFY(DEPENDENT_VARIABLE)
    NULLIFY(BOUNDARY_CONDITIONS_VARIABLE)
    NULLIFY(DIRICHLET_BOUNDARY_CONDITIONS)
    NULLIFY(FULL_LOADS)
    NULLIFY(PREV_LOADS)
    NULLIFY(CURRENT_LOADS)

    myComputationalNodeNumber=ComputationalEnvironment_NodeNumberGet(err,error)
    
    !Take the stored load, scale it down appropriately then apply to the unknown variables
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(DIAGNOSTICS1) THEN
        CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  equations set",EQUATIONS_SET%USER_NUMBER,err,error,*999)
      ENDIF
      IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
        dependentField=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(dependentField)) THEN
          IF(ALLOCATED(dependentField%VARIABLES)) THEN
            !Loop over the variables associated with this equations set
            !\todo: Looping over all field variables is not safe when volume-coupled problem is solved. Look at matrix and rhs mapping instead?
            DO variable_idx=1,dependentField%NUMBER_OF_VARIABLES
              DEPENDENT_VARIABLE=>dependentField%variables(variable_idx)
              variable_type=DEPENDENT_VARIABLE%VARIABLE_TYPE
              CALL BOUNDARY_CONDITIONS_VARIABLE_GET(BOUNDARY_CONDITIONS,DEPENDENT_VARIABLE,BOUNDARY_CONDITIONS_VARIABLE, &
                & err,error,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE)) THEN
                DOMAIN_MAPPING=>DEPENDENT_VARIABLE%DOMAIN_MAPPING
                IF(ASSOCIATED(DOMAIN_MAPPING)) THEN

                  ! Check if there are any incremented conditions applied for this boundary conditions variable
                  IF(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_FIXED_INCREMENTED)>0.OR. &
                      & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)>0) THEN
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS)) THEN
                      DIRICHLET_BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS_VARIABLE%DIRICHLET_BOUNDARY_CONDITIONS
                      !Get the pointer to vector holding the full and current loads
                      !   full load: FIELD_BOUNDARY_CONDITIONS_SET_TYPE - holds the target load values
                      !   current load: FIELD_VALUES_SET_TYPE - holds the current increment values
                      CALL Field_ParameterSetDataGet(dependentField,variable_type,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                        & FULL_LOADS,err,error,*999)
                      !chrm 22/06/2010: 'FIELD_BOUNDARY_CONDITIONS_SET_TYPE' does not get updated with time (update_BCs)
                      !\ToDo: How can this be achieved ???
  !                     write(*,*)'FULL_LOADS = ',FULL_LOADS
                      CALL Field_ParameterSetDataGet(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                        & CURRENT_LOADS,err,error,*999)
  !                     write(*,*)'CURRENT_LOADS = ',CURRENT_LOADS
                      !Get full increment, calculate new load, then apply to dependent field
                      DO dirichlet_idx=1,BOUNDARY_CONDITIONS_VARIABLE%NUMBER_OF_DIRICHLET_CONDITIONS
                        dirichlet_dof_idx=DIRICHLET_BOUNDARY_CONDITIONS%DIRICHLET_DOF_INDICES(dirichlet_idx)
                        !Check whether we have an incremented boundary condition type
                        SELECT CASE(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(dirichlet_dof_idx))
                        CASE(BOUNDARY_CONDITION_FIXED_INCREMENTED, &
                            & BOUNDARY_CONDITION_MOVED_WALL_INCREMENTED)
                          !Convert dof index to local index
                          IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(dirichlet_dof_idx)%DOMAIN_NUMBER(1)== &
                            & myComputationalNodeNumber) THEN
                            dirichlet_dof_idx=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(dirichlet_dof_idx)%LOCAL_NUMBER(1)
                            IF(0<dirichlet_dof_idx.AND.dirichlet_dof_idx<DOMAIN_MAPPING%GHOST_START) THEN
                              FULL_LOAD=FULL_LOADS(dirichlet_dof_idx)
                              ! Apply full load if last step, or fixed BC
                              IF(ITERATION_NUMBER==MAXIMUM_NUMBER_OF_ITERATIONS) THEN
                                CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                                  & dirichlet_dof_idx,FULL_LOAD,err,error,*999)
                              ELSE
                                !Calculate new load and apply to dependent field
                                CURRENT_LOAD=CURRENT_LOADS(dirichlet_dof_idx)
                                NEW_LOAD=CURRENT_LOAD+(FULL_LOAD-CURRENT_LOAD)/(MAXIMUM_NUMBER_OF_ITERATIONS-ITERATION_NUMBER+1)
                                CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                                  & dirichlet_dof_idx,NEW_LOAD,err,error,*999)
                                IF(DIAGNOSTICS1) THEN
                                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx",dirichlet_dof_idx,err,error,*999)
                                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load",CURRENT_LOAD,err,error,*999)
                                  CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",NEW_LOAD,err,error,*999)
                                ENDIF
                              ENDIF !Full or intermediate load
                            ENDIF !non-ghost dof
                          ENDIF !current domain
                        CASE DEFAULT
                          !Do nothing for non-incremented boundary conditions
                        END SELECT
                      ENDDO !dirichlet_idx
  !---tob
                      !\ToDo: What happens if the call below is issued
                      !without actually that the dependent field has been modified in above conditional ?
                      CALL Field_ParameterSetUpdateStart(dependentField, &
                        & variable_type, FIELD_VALUES_SET_TYPE,err,error,*999)
                      CALL Field_ParameterSetUpdateFinish(dependentField, &
                        & variable_type, FIELD_VALUES_SET_TYPE,err,error,*999)
  !---toe
                      !Restore the vector handles
                      CALL Field_ParameterSetDataRestore(dependentField,variable_type,FIELD_BOUNDARY_CONDITIONS_SET_TYPE, &
                        & FULL_LOADS,err,error,*999)
                      CALL Field_ParameterSetDataRestore(dependentField,variable_type,FIELD_VALUES_SET_TYPE, &
                        & CURRENT_LOADS,err,error,*999)
                    ELSE
                      localError="Dirichlet boundary condition for variable type "// &
                        & TRIM(NumberToVString(variable_type,"*",err,error))//" is not associated."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF

                  ! Also increment any incremented Neumann point conditions
                  IF(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)>0) THEN
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions)) THEN
                      ! The boundary conditions parameter set contains the full values and the
                      ! current incremented values are transferred to the point values vector
                      DO condition_idx=1,BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED)+ &
                          & BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_NEUMANN_POINT)
                        condition_global_dof=BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions%setDofs(condition_idx)
                        ! condition_global_dof could be for non-incremented point Neumann condition
                        IF(BOUNDARY_CONDITIONS_VARIABLE%CONDITION_TYPES(condition_global_dof)/= &
                          & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED) CYCLE
                        IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)%DOMAIN_NUMBER(1)== &
                          & myComputationalNodeNumber) THEN
                          condition_local_dof=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)% &
                            & LOCAL_NUMBER(1)
                          neumann_point_dof=BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions%pointDofMapping% &
                            & GLOBAL_TO_LOCAL_MAP(condition_idx)%LOCAL_NUMBER(1)
                          CALL FIELD_PARAMETER_SET_GET_LOCAL_DOF(dependentField,variable_type, &
                            & FIELD_BOUNDARY_CONDITIONS_SET_TYPE,condition_local_dof,FULL_LOAD,err,error,*999)
                          CALL DistributedVector_ValuesSet(BOUNDARY_CONDITIONS_VARIABLE%neumannBoundaryConditions% &
                            & pointValues,neumann_point_dof,FULL_LOAD*(REAL(ITERATION_NUMBER)/REAL(MAXIMUM_NUMBER_OF_ITERATIONS)), &
                            & err,error,*999)
                        END IF
                      END DO
                    ELSE
                      localError="Neumann boundary conditions for variable type "// &
                        & TRIM(NumberToVString(variable_type,"*",err,error))//" are not associated even though"// &
                        & TRIM(NumberToVString(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS( &
                        & BOUNDARY_CONDITION_NEUMANN_POINT_INCREMENTED), &
                        & '*',err,error))//" conditions of this type has been counted."
                      CALL FlagError(localError,err,error,*999)
                    END IF
                  END IF

                  !There might also be pressure incremented conditions
                  IF (BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED)>0) THEN
                    ! handle pressure incremented boundary conditions
                    IF(ASSOCIATED(BOUNDARY_CONDITIONS_VARIABLE%PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS)) THEN
                      PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS=>BOUNDARY_CONDITIONS_VARIABLE% &
                        & PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS
                      !Due to a variety of reasons, the pressure incremented type is setup differently to dirichlet conditions.
                      !We store two sets of vectors, the current and previous values
                      !   current: FIELD_PRESSURE_VALUES_SET_TYPE - always holds the current increment, even if not incremented
                      !   previous: FIELD_PREVIOUS_PRESSURE_SET_TYPE - holds the previously applied increment
                      !Grab the pointers for both
                      CALL Field_ParameterSetDataGet(dependentField,variable_type,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                        & PREV_LOADS,err,error,*999)
                      CALL Field_ParameterSetDataGet(dependentField,variable_type,FIELD_PRESSURE_VALUES_SET_TYPE, &
                        & CURRENT_LOADS,err,error,*999)
                      !Calculate the new load, update the old load
                      IF(ITERATION_NUMBER==1) THEN
                        !On the first iteration, FIELD_PRESSURE_VALUES_SET_TYPE actually contains the full load
                        DO condition_idx=1,BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS( &
                            & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          !Global dof index
                          condition_global_dof=PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS%PRESSURE_INCREMENTED_DOF_INDICES &
                            & (condition_idx)
                          !Must convert into local dof index
                          IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)%DOMAIN_NUMBER(1)== &
                            & myComputationalNodeNumber) THEN
                            condition_local_dof=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)% &
                              & LOCAL_NUMBER(1)
                            IF(0<condition_local_dof.AND.condition_local_dof<DOMAIN_MAPPING%GHOST_START) THEN
                              NEW_LOAD=CURRENT_LOADS(condition_local_dof)
                              NEW_LOAD=NEW_LOAD/MAXIMUM_NUMBER_OF_ITERATIONS
!if (condition_idx==1) write(*,*) "new load=",new_load
                              !Update current and previous loads
                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                & FIELD_PRESSURE_VALUES_SET_TYPE,condition_local_dof,NEW_LOAD,err,error,*999)
                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                & FIELD_PREVIOUS_PRESSURE_SET_TYPE,condition_local_dof,0.0_dp,err,error,*999)
                              IF(DIAGNOSTICS1) THEN
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx", &
                                    & condition_local_dof,err,error,*999)
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load", &
                                    & CURRENT_LOADS(condition_local_dof),err,error,*999)
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",NEW_LOAD,err,error,*999)
                              ENDIF
                            ENDIF !Non-ghost dof
                          ENDIF !Current domain
                        ENDDO !condition_idx
                      ELSE
                        !Calculate the new load, keep the current load
                        DO condition_idx=1,BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS( &
                            & BOUNDARY_CONDITION_PRESSURE_INCREMENTED)
                          !This is global dof idx
                          condition_global_dof=PRESSURE_INCREMENTED_BOUNDARY_CONDITIONS%PRESSURE_INCREMENTED_DOF_INDICES &
                            & (condition_idx)
                          !Must convert into local dof index
                          IF(DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)%DOMAIN_NUMBER(1)== &
                            & myComputationalNodeNumber) THEN
                            condition_local_dof=DOMAIN_MAPPING%GLOBAL_TO_LOCAL_MAP(condition_global_dof)% &
                              & LOCAL_NUMBER(1)
                            IF(0<condition_local_dof.AND.condition_local_dof<DOMAIN_MAPPING%GHOST_START) THEN
                              PREV_LOAD=PREV_LOADS(condition_local_dof)
                              CURRENT_LOAD=CURRENT_LOADS(condition_local_dof)
                              NEW_LOAD=CURRENT_LOAD+(CURRENT_LOAD-PREV_LOAD)  !This may be subject to numerical errors...
!if (condition_idx==1) write(*,*) "new load=",new_load
                              !Update current and previous loads
                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                & FIELD_PRESSURE_VALUES_SET_TYPE,condition_local_dof,NEW_LOAD,err,error,*999)
                              CALL Field_ParameterSetUpdateLocalDOF(dependentField,variable_type, &
                                & FIELD_PREVIOUS_PRESSURE_SET_TYPE,condition_local_dof,CURRENT_LOAD,err,error,*999)
                              IF(DIAGNOSTICS1) THEN
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"  dof idx", &
                                    & condition_local_dof,err,error,*999)
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    current load", &
                                    & CURRENT_LOADS(condition_local_dof),err,error,*999)
                                CALL WriteStringValue(DIAGNOSTIC_OUTPUT_TYPE,"    new load",NEW_LOAD,err,error,*999)
                              ENDIF
                            ENDIF !Non-ghost dof
                          ENDIF !Current domain
                        ENDDO !condition_idx
                      ENDIF
                      !Start transfer of dofs to neighbouring domains
                      CALL Field_ParameterSetUpdateStart(dependentField,variable_type,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                        & err,error,*999)
                      CALL Field_ParameterSetUpdateStart(dependentField,variable_type,FIELD_PRESSURE_VALUES_SET_TYPE, &
                        & err,error,*999)
                      !Restore the vector handles
                      CALL Field_ParameterSetDataRestore(dependentField,variable_type,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                        & PREV_LOADS,err,error,*999)
                      CALL Field_ParameterSetDataRestore(dependentField,variable_type,FIELD_PRESSURE_VALUES_SET_TYPE, &
                        & CURRENT_LOADS,err,error,*999)
                      !Finish transfer of dofs to neighbouring domains
                      CALL Field_ParameterSetUpdateFinish(dependentField,variable_type,FIELD_PREVIOUS_PRESSURE_SET_TYPE, &
                        & err,error,*999)
                      CALL Field_ParameterSetUpdateFinish(dependentField,variable_type,FIELD_PRESSURE_VALUES_SET_TYPE, &
                        & err,error,*999)
                    ELSE
                      localError="Pressure incremented boundary condition for variable type "// &
                        & TRIM(NumberToVString(variable_type,"*",err,error))//" is not associated even though"// &
                        & TRIM(NumberToVString(BOUNDARY_CONDITIONS_VARIABLE%DOF_COUNTS(BOUNDARY_CONDITION_PRESSURE_INCREMENTED), &
                        & '*',err,error))//" conditions of this type has been counted."
                      CALL FlagError(localError,err,error,*999)
                    ENDIF
                  ENDIF !Pressure incremented bc block
                ELSE
                  localError="Domain mapping is not associated for variable "// &
                    & TRIM(NumberToVString(variable_type,"*",err,error))//" of dependent field"
                  CALL FlagError(localError,err,error,*999)
                ENDIF !Domain mapping test
              ELSE
                ! do nothing - no boundary conditions variable type associated?
              ENDIF
            ENDDO !variable_idx
          ELSE
            CALL FlagError("Dependent field variables are not allocated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Dependent field is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Boundary conditions are not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF

    EXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT",err,error)
    RETURN 1

  END SUBROUTINE EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT

  !
  !================================================================================================================================
  !

  !> Apply load increments for equations sets
  SUBROUTINE EQUATIONS_SET_LOAD_INCREMENT_APPLY(EQUATIONS_SET,BOUNDARY_CONDITIONS,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS, &
    & err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS !<The boundary conditions to apply the increment to
    INTEGER(INTG), INTENT(IN) :: ITERATION_NUMBER !<The current load increment iteration index
    INTEGER(INTG), INTENT(IN) :: MAXIMUM_NUMBER_OF_ITERATIONS !<Final index for load increment loop
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string

    ENTERS("EQUATIONS_SET_LOAD_INCREMENT_APPLY",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      !Increment boundary conditions
      CALL EQUATIONS_SET_BOUNDARY_CONDITIONS_INCREMENT(EQUATIONS_SET,BOUNDARY_CONDITIONS,ITERATION_NUMBER, &
        & MAXIMUM_NUMBER_OF_ITERATIONS,err,error,*999)

      !Apply any other equation set specific increments
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<1) THEN
        CALL FlagError("Equations set specification must have at least one entry.",err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(1))
      CASE(EQUATIONS_SET_ELASTICITY_CLASS)
        CALL ELASTICITY_LOAD_INCREMENT_APPLY(EQUATIONS_SET,ITERATION_NUMBER,MAXIMUM_NUMBER_OF_ITERATIONS,err,error,*999)
      CASE DEFAULT
        !Do nothing
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF


    EXITS("EQUATIONS_SET_LOAD_INCREMENT_APPLY")
    RETURN
999 ERRORSEXITS("EQUATIONS_SET_LOAD_INCREMENT_APPLY",err,error)
    RETURN 1

  END SUBROUTINE EQUATIONS_SET_LOAD_INCREMENT_APPLY

 !
  !================================================================================================================================
  !

  !>Assembles the equations stiffness matrix, residuals and rhs for a nonlinear static equations set using a nodal method
  SUBROUTINE EquationsSet_AssembleStaticNonlinearNodal(equationsSet,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to assemble the equations for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: nodalMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
    
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

    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatrices_NodalInitialise(vectorMatrices,err,error,*999)
    nodalMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%nodes
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
    DO nodeIdx=nodalMapping%INTERNAL_START,nodalMapping%INTERNAL_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_NodeAdd(vectorMatrices,err,error,*999)
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
    DO nodeIdx=nodalMapping%BOUNDARY_START,nodalMapping%GHOST_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_NodeAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_NodalFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: nodeNumber !<The node number to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code 
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: matrixIdx
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(NodalMatrixType), POINTER :: nodalMatrix
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("EquationsSet_NodalJacobianEvaluate",err,error,*999)

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
    
    DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
      SELECT CASE(nonlinearMatrices%jacobians(matrixIdx)%ptr%jacobianCalculationType)
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
          & TRIM(NumberToVString(nonlinearMatrices%jacobians(matrixIdx)%ptr%jacobianCalculationType,"*",err,error))// &
          & " is not valid for matrix index number "//TRIM(NumberToVString(matrixIdx,"*",err,error))//"."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    END DO !matrixIdx
    IF(equations%outputType>=EQUATIONS_NODAL_MATRIX_OUTPUT) THEN
      CALL WriteString(GENERAL_OUTPUT_TYPE,"",err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal Jacobian matrix:",err,error,*999)
      CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"Node number = ",nodeNumber,err,error,*999)
      CALL WriteString(GENERAL_OUTPUT_TYPE,"Nodal Jacobian:",err,error,*999)
      DO matrixIdx=1,nonlinearMatrices%numberOfJacobians
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Jacobian number = ",matrixIdx,err,error,*999)
        CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Update Jacobian = ",nonlinearMatrices%jacobians(matrixIdx)%ptr% &
          & updateJacobian,err,error,*999)
        IF(nonlinearMatrices%jacobians(matrixIdx)%ptr%updateJacobian) THEN
          nodalMatrix=>nonlinearMatrices%jacobians(matrixIdx)%ptr%nodalJacobian
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of rows = ",nodalMatrix%numberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Number of columns = ",nodalMatrix%numberOfColumns,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of rows = ",nodalMatrix%maxNumberOfRows,err,error,*999)
          CALL WriteStringValue(GENERAL_OUTPUT_TYPE,"  Maximum number of columns = ",nodalMatrix%maxNumberOfColumns,err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,8,8,nodalMatrix%rowDofs, &
            & '("  Row dofs     :",8(X,I13))','(16X,8(X,I13))',err,error,*999)
          CALL WriteStringVector(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfColumns,8,8,nodalMatrix% &
            & columnDofs,'("  Column dofs  :",8(X,I13))','(16X,8(X,I13))',err,error,*999)
          CALL WriteStringMatrix(GENERAL_OUTPUT_TYPE,1,1,nodalMatrix%numberOfRows,1,1,nodalMatrix% &
            & numberOfColumns,8,8,nodalMatrix%matrix(1:nodalMatrix%numberOfRows,1:nodalMatrix% &
            & numberOfColumns),WRITE_STRING_MATRIX_NAME_AND_INDICES,'("  Matrix','(",I2,",:)',' :",8(X,E13.6))', &
            & '(16X,8(X,E13.6))',err,error,*999)
        END IF
      END DO !matrixIdx
    END IF

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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the Jacobian for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: nodalMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField
  
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_JACOBIAN_ONLY,0.0_DP,err,error,*999)
    !Assemble the nodes
    !Allocate the nodal matrices 
    CALL EquationsMatrices_NodalInitialise(vectorMatrices,err,error,*999)
    nodalMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%nodes
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
    DO nodeIdx=nodalMapping%INTERNAL_START,nodalMapping%INTERNAL_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_JacobianNodeAdd(vectorMatrices,err,error,*999)
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
    DO nodeIdx=nodalMapping%BOUNDARY_START,nodalMapping%GHOST_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalJacobianEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_JacobianNodeAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_NodalFinalise(vectorMatrices,err,error,*999)
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
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to evaluate the residual for
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: numberOfTimes
    INTEGER(INTG) :: nodeIdx,nodeNumber
    REAL(SP) :: nodeUserElapsed,nodeSystemElapsed,userElapsed,userTime1(1),userTime2(1),userTime3(1),userTime4(1), &
      & userTime5(1),userTime6(1),systemElapsed,systemTime1(1),systemTime2(1),systemTime3(1),systemTime4(1), &
      & systemTime5(1),systemTime6(1)
    TYPE(DOMAIN_MAPPING_TYPE), POINTER :: nodalMapping
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField
 
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
    
    IF(equations%outputType>=EQUATIONS_TIMING_OUTPUT) THEN
      CALL CPUTimer(USER_CPU,userTime1,err,error,*999)
      CALL CPUTimer(SYSTEM_CPU,systemTime1,err,error,*999)
    ENDIF
    !Initialise the matrices and rhs vector
    CALL EquationsMatrices_VectorValuesInitialise(vectorMatrices,EQUATIONS_MATRICES_NONLINEAR_ONLY,0.0_DP,err,error,*999)
    !Allocate the nodal matrices 
    CALL EquationsMatrices_NodalInitialise(vectorMatrices,err,error,*999)
    nodalMapping=>dependentField%decomposition%domain(dependentField%decomposition%MESH_COMPONENT_NUMBER)%ptr%mappings%nodes
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
    DO nodeIdx=nodalMapping%INTERNAL_START,nodalMapping%INTERNAL_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_NodeAdd(vectorMatrices,err,error,*999)
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
    DO nodeIdx=nodalMapping%BOUNDARY_START,nodalMapping%GHOST_FINISH
      nodeNumber=nodalMapping%DOMAIN_LIST(nodeIdx)
      numberOfTimes=numberOfTimes+1
      CALL EquationsMatrices_NodalCalculate(vectorMatrices,nodeNumber,err,error,*999)
      CALL EquationsSet_NodalResidualEvaluate(equationsSet,nodeNumber,err,error,*999)
      CALL EquationsMatrices_NodeAdd(vectorMatrices,err,error,*999)
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
    CALL EquationsMatrices_NodalFinalise(vectorMatrices,err,error,*999)
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
      CALL EquationsMatrices_VectorOutput(GENERAL_OUTPUT_TYPE,vectorMatrices,err,error,*999)
    ENDIF
       
    EXITS("EquationsSet_ResidualEvaluateStaticNodal")
    RETURN
999 ERRORSEXITS("EquationsSet_ResidualEvaluateStaticNodal",err,error)
    RETURN 1
    
  END SUBROUTINE EquationsSet_ResidualEvaluateStaticNodal

  !
  !================================================================================================================================
  !

END MODULE EQUATIONS_SET_ROUTINES
