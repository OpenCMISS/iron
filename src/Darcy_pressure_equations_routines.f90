!> \file
!> \author Adam Reeve
!> \brief This module handles all Darcy pressure equations routines to be coupled with finite elasticity.
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

!>This module handles all Darcy pressure equations routines.
MODULE DARCY_PRESSURE_EQUATIONS_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE CONTROL_LOOP_ROUTINES
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
  USE FINITE_ELASTICITY_ROUTINES
  USE INPUT_OUTPUT
  USE ISO_VARYING_STRING
  USE Kinds
  USE Maths
  USE MatrixVector
  USE NODE_ROUTINES
  USE PROBLEM_CONSTANTS
  USE Strings
  USE SOLVER_ROUTINES
  USE Timer
  USE Types

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC DarcyPressure_FiniteElementResidualEvaluate
  
  PUBLIC DarcyPressure_EquationsSetSetup
  
  PUBLIC DarcyPressure_EquationsSetSolutionMethodSet
  
  PUBLIC DarcyPressure_EquationsSetSpecificationSet

CONTAINS

  !
  !================================================================================================================================
  !

  !>Calculates the element residual vector and RHS for a Darcy pressure equation finite element equations set.
  SUBROUTINE DarcyPressure_FiniteElementResidualEvaluate(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,PGMSI(3),PGNSI(3)
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesNonlinearType), POINTER :: nonlinearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: dependentField,geometricField
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(FIELD_INTERPOLATED_POINT_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT,SOLID_DEPENDENT_INTERPOLATED_POINT, &
      & MATERIALS_INTERPOLATED_POINT,FIBRE_INTERPOLATED_POINT
    TYPE(FIELD_INTERPOLATED_POINT_METRICS_TYPE), POINTER :: GEOMETRIC_INTERPOLATED_POINT_METRICS, &
      & SOLID_DEPENDENT_INTERPOLATED_POINT_METRICS
    TYPE(FIELD_INTERPOLATION_PARAMETERS_TYPE), POINTER :: GEOMETRIC_INTERPOLATION_PARAMETERS,DEPENDENT_INTERPOLATION_PARAMETERS, &
      & FIBRE_INTERPOLATION_PARAMETERS,MATERIALS_INTERPOLATION_PARAMETERS,SOLID_DEPENDENT_INTERPOLATION_PARAMETERS
    REAL(DP), POINTER :: elementResidualVector(:)
    REAL(DP) :: K(3,3),density
    REAL(DP) :: DZDNU(3,3),SIGMA(3,3),DNUDZ(3,3),DNUDZT(3,3),TEMP_MATRIX(3,3)
    REAL(DP) :: Jznu
    INTEGER(INTG) :: SOLID_COMPONENT_NUMBER,SOLID_NUMBER_OF_XI,NUMBER_OF_DIMENSIONS
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_FiniteElementResidualEvaluate",err,error,*999)

    NULLIFY(GEOMETRIC_INTERPOLATED_POINT)
    NULLIFY(GEOMETRIC_INTERPOLATED_POINT_METRICS)
    NULLIFY(SOLID_DEPENDENT_INTERPOLATED_POINT)
    NULLIFY(SOLID_DEPENDENT_INTERPOLATED_POINT_METRICS)
    NULLIFY(MATERIALS_INTERPOLATED_POINT)
    NULLIFY(FIBRE_INTERPOLATED_POINT)
    NULLIFY(DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(SOLID_DEPENDENT_INTERPOLATION_PARAMETERS)
    NULLIFY(GEOMETRIC_INTERPOLATION_PARAMETERS)
    NULLIFY(FIBRE_INTERPOLATION_PARAMETERS)
    NULLIFY(MATERIALS_INTERPOLATION_PARAMETERS)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) &
        & CALL FlagError("Equations set specification must have three entries for a Darcy pressure type equations set.", &
        & err,error,*999)
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
            & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
          dependentField=>equations%interpolation%dependentField
          geometricField=>equations%interpolation%geometricField
          vectorMatrices=>vectorEquations%vectorMatrices
          nonlinearMatrices=>vectorMatrices%nonlinearMatrices
          elementResidualVector=>nonlinearMatrices%elementResidual%VECTOR
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping
          FIELD_VARIABLE=>dependentField%VARIABLE_TYPE_MAP(FIELD_V_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>dependentField%DECOMPOSITION%DOMAIN(FIELD_VARIABLE%COMPONENTS(1)%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>geometricField%DECOMPOSITION%DOMAIN(geometricField%DECOMPOSITION%MESH_COMPONENT_NUMBER)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr
          SOLID_DEPENDENT_INTERPOLATION_PARAMETERS=>equations%interpolation%dependentInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
          GEOMETRIC_INTERPOLATION_PARAMETERS=>equations%interpolation%geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
          FIBRE_INTERPOLATION_PARAMETERS=>equations%interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr
          MATERIALS_INTERPOLATION_PARAMETERS=>equations%interpolation%materialsInterpParameters(FIELD_U1_VARIABLE_TYPE)%ptr

          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & GEOMETRIC_INTERPOLATION_PARAMETERS,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & FIBRE_INTERPOLATION_PARAMETERS,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & MATERIALS_INTERPOLATION_PARAMETERS,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER, &
            & SOLID_DEPENDENT_INTERPOLATION_PARAMETERS,err,error,*999)

          GEOMETRIC_INTERPOLATED_POINT=>equations%interpolation%geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          GEOMETRIC_INTERPOLATED_POINT_METRICS=>equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr
          FIBRE_INTERPOLATED_POINT=>equations%interpolation%fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          MATERIALS_INTERPOLATED_POINT=>equations%interpolation%materialsInterpPoint(FIELD_U1_VARIABLE_TYPE)%ptr
          SOLID_DEPENDENT_INTERPOLATED_POINT=>equations%interpolation%dependentInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr
          SOLID_DEPENDENT_INTERPOLATED_POINT_METRICS=>equations%interpolation%dependentInterpPointMetrics( &
            & FIELD_U_VARIABLE_TYPE)%ptr

          SOLID_COMPONENT_NUMBER=dependentField%VARIABLE_TYPE_MAP(FIELD_U_VARIABLE_TYPE)%ptr%COMPONENTS(1)%MESH_COMPONENT_NUMBER
          SOLID_NUMBER_OF_XI=dependentField%DECOMPOSITION%DOMAIN(SOLID_COMPONENT_NUMBER)%ptr%TOPOLOGY%ELEMENTS% &
            & ELEMENTS(ELEMENT_NUMBER)%BASIS%NUMBER_OF_XI
          NUMBER_OF_DIMENSIONS=GEOMETRIC_BASIS%NUMBER_OF_XI

          IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
          ENDIF

          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
            !Calculate RWG.
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)

            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & SOLID_DEPENDENT_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & GEOMETRIC_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & FIBRE_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATE_GAUSS(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng, &
              & MATERIALS_INTERPOLATED_POINT,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(DEPENDENT_BASIS%NUMBER_OF_XI, &
              & SOLID_DEPENDENT_INTERPOLATED_POINT_METRICS,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,err,error,*999)

            !Get the permeability tensor
            K(1,1)=MATERIALS_INTERPOLATED_POINT%VALUES(1,1)
            K(1,2)=MATERIALS_INTERPOLATED_POINT%VALUES(2,1)
            K(1,3)=MATERIALS_INTERPOLATED_POINT%VALUES(3,1)
            K(2,2)=MATERIALS_INTERPOLATED_POINT%VALUES(4,1)
            K(2,3)=MATERIALS_INTERPOLATED_POINT%VALUES(5,1)
            K(3,3)=MATERIALS_INTERPOLATED_POINT%VALUES(6,1)
            K(2,1)=K(1,2)
            K(3,1)=K(1,3)
            K(3,2)=K(2,3)
            density=MATERIALS_INTERPOLATED_POINT%VALUES(7,1)

            !Calculate F=dZ/dNU, the deformation gradient tensor at the gauss point
            CALL FiniteElasticity_GaussDeformationGradientTensor(SOLID_DEPENDENT_INTERPOLATED_POINT_METRICS, &
              & GEOMETRIC_INTERPOLATED_POINT_METRICS,FIBRE_INTERPOLATED_POINT,DZDNU,ERR,ERROR,*999)

            CALL Invert(DZDNU,DNUDZ,Jznu,err,error,*999)
            CALL MatrixTranspose(DNUDZ,DNUDZT,err,error,*999)
            CALL MatrixProduct(K,DNUDZT,TEMP_MATRIX,err,error,*999)
            CALL MatrixProduct(DNUDZ,TEMP_MATRIX,SIGMA,err,error,*999)
            SIGMA=density * Jznu * SIGMA

            IF(DIAGNOSTICS1) THEN
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,DZDNU,WRITE_STRING_MATRIX_NAME_AND_INDICES,'(" DZDNU','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',err,error,*999)
              CALL WRITE_STRING_VALUE(DIAGNOSTIC_OUTPUT_TYPE,"Jznu",Jznu,err,error,*999)
              CALL WRITE_STRING_MATRIX(DIAGNOSTIC_OUTPUT_TYPE,1,1,3,1,1,3, &
                & 3,3,SIGMA,WRITE_STRING_MATRIX_NAME_AND_INDICES,'(" SIGMA','(",I1,",:)',' :",3(X,E13.6))', &
                & '(17X,3(X,E13.6))',err,error,*999)
            ENDIF

            CALL MatrixProduct(SIGMA,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%GU,TEMP_MATRIX,err,error,*999)

            !Loop over field components
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(nonlinearMatrices%updateResidual) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      ENDDO !ni
                      DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
                            elementResidualVector(mhs)=elementResidualVector(mhs)+ &
                              & PGMSI(mi)*PGNSI(ni)*TEMP_MATRIX(mi,ni)*DEPENDENT_INTERPOLATION_PARAMETERS%PARAMETERS(ns,1)* &
                              & equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr% &
                              & SCALE_FACTORS(ms,mh)*equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)% &
                              & PTR%SCALE_FACTORS(ns,nh)*RWG
                          ELSE
                            elementResidualVector(mhs)=elementResidualVector(mhs)+ &
                              & PGMSI(mi)*PGNSI(ni)*TEMP_MATRIX(mi,ni)*RWG*DEPENDENT_INTERPOLATION_PARAMETERS%PARAMETERS(ns,1)
                          ENDIF
                        ENDDO !ni
                      ENDDO !mi
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
          ENDDO !ng
          !Scale factor adjustment
          IF(dependentField%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr,err,error,*999)
            mhs=0
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_V_VARIABLE_TYPE)%ptr%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE DEFAULT
          localError="The third equations set specification of "// &
            & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a coupled Darcy fluid pressure type of a fluid mechanicsequations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DarcyPressure_FiniteElementResidualEvaluate")
    RETURN
999 ERRORSEXITS("DarcyPressure_FiniteElementResidualEvaluate",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_FiniteElementResidualEvaluate

  !
  !================================================================================================================================
  !

  !>Sets up the Darcy pressure equation type of a fluid mechanics equations set class.
  SUBROUTINE DarcyPressure_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Darcy pressure equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS,component_idx
    INTEGER(INTG) :: NUMBER_OF_DARCY_COMPONENTS,NUMBER_OF_COMPONENTS,NUMBER_OF_SOLID_COMPONENTS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_EquationsSetSetup",err,error,*999)

    NULLIFY(GEOMETRIC_DECOMPOSITION)
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(EQUATIONS_MATERIALS)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy pressure type equations set.", &
          & err,error,*999)
      END IF
      NUMBER_OF_DIMENSIONS = EQUATIONS_SET%REGION%COORDINATE_SYSTEM%NUMBER_OF_DIMENSIONS
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL DarcyPressure_EquationsSetSolutionMethodSet(EQUATIONS_SET, &
                & EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy pressure equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
              & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
            NUMBER_OF_COMPONENTS=NUMBER_OF_DIMENSIONS !Solid components
            NUMBER_OF_DARCY_COMPONENTS=1 !Only solving for the fluid pressure at the moment
          END SELECT
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              !Create the auto created dependent field
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & DEPENDENT_FIELD,err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL FIELD_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,"Dependent Field",err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,4,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,(/FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/),err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                & NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

              !Set labels
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"U",err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del U/del n", &
                & err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE,"V",err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE,"del V/del n", &
                & err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"x1",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,2,"x2",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,3,"x3",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & "del x1/del n",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,2, &
                & "del x2/del n",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,3, &
                & "del x3/del n",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"p",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & "del p/del n",err,error,*999)

              !Elasticity: Default to the geometric interpolation setup
              DO component_idx=1,NUMBER_OF_DIMENSIONS
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ENDDO !component_idx
              !Darcy: Default pressure and mass increase to the first geometric component
              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ENDDO !component_idx

              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Elasticity: Set the displacement components to node based interpolation
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !component_idx
                !Darcy: Set the pressure and mass increase components to node based interpolation
                DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_V_VARIABLE_TYPE, &
                    & component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD, &
                    & FIELD_DELVDELN_VARIABLE_TYPE,component_idx,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(localError,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GEOMETRIC_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,(/FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                    & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE/) &
                  & ,err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)

              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & NUMBER_OF_DIMENSIONS,err,error,*999)

              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & NUMBER_OF_DARCY_COMPONENTS,err,error,*999)

              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                !Elasticity:
                DO component_idx=1,NUMBER_OF_DIMENSIONS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                ENDDO !component_idx
                !Darcy:
                DO component_idx=1,NUMBER_OF_DARCY_COMPONENTS
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,component_idx, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
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
              & " is invalid for a finite elasticity equation"
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
          NUMBER_OF_COMPONENTS=8 !permeability tensor + density + original porosity
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE)
            NUMBER_OF_SOLID_COMPONENTS=6
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE)
            NUMBER_OF_SOLID_COMPONENTS=4
          CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
            NUMBER_OF_SOLID_COMPONENTS=6
          CASE DEFAULT
            localError="The third equations set specification of "// &
              & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Darcy pressure type of a fluid mechanics equation set."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
              !Create the auto created materials field, which is shared with the solid equations
              CALL FIELD_CREATE_START(EQUATIONS_SET_SETUP%FIELD_USER_NUMBER,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                & MATERIALS_FIELD,err,error,*999)
              CALL FIELD_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL FIELD_MESH_DECOMPOSITION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL FIELD_GEOMETRIC_FIELD_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,EQUATIONS_SET%GEOMETRY% &
                & GEOMETRIC_FIELD,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,3,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE],err,error,*999)

              !U variable, elasticity parameters
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                 & NUMBER_OF_SOLID_COMPONENTS,err,error,*999)

              !V variable, solid density
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                 & 1,err,error,*999)

              !U1 variable, fluid parameters
              CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                 & NUMBER_OF_COMPONENTS,err,error,*999)

              CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              DO component_idx=1,NUMBER_OF_SOLID_COMPONENTS
                !Default the materials components to the geometric interpolation setup with constant interpolation
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ENDDO
              !Solid density field variable
              CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                & 1,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_V_VARIABLE_TYPE, &
                & 1,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              DO component_idx=1,NUMBER_OF_COMPONENTS
                !Default the materials components to the geometric interpolation setup with constant interpolation
                CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & component_idx,GEOMETRIC_MESH_COMPONENT,err,error,*999)
              ENDDO

              !Default the field scaling to that of the geometric field
              CALL FIELD_SCALING_TYPE_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
              CALL FIELD_SCALING_TYPE_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,GEOMETRIC_SCALING_TYPE,err,error,*999)
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,3,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                & FIELD_U1_VARIABLE_TYPE],err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_SOLID_COMPONENTS,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE, &
                  & 1,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & NUMBER_OF_COMPONENTS,err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED) THEN
                !Finish creating the materials field
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                !Default to K=I
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,1,1.0_DP,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,2,0.0_DP,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,3,0.0_DP,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,4,1.0_DP,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,5,0.0_DP,err,error,*999)
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,6,1.0_DP,err,error,*999)
                !Density
                CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_VALUES_SET_TYPE,7,1.0_DP,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy pressure equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a finite elasticity fluid pressure equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_NONLINEAR,err,error,*999)
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_STATIC,err,error,*999)
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,0,err,error,*999)
              CALL EquationsMapping_ResidualVariablesNumberSet(vectorMapping,2,err,error,*999)
              CALL EquationsMapping_ResidualVariableTypesSet(vectorMapping, &
                  & [FIELD_V_VARIABLE_TYPE,FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(equations%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_BLOCK_STORAGE_TYPE,err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                CALL EquationsMatrices_NonlinearStorageTypeSet(vectorMatrices,MATRIX_COMPRESSED_ROW_STORAGE_TYPE, &
                  & err,error,*999)
                CALL EquationsMatrices_NonlinearStructureTypeSet(vectorMatrices,EQUATIONS_MATRIX_FEM_STRUCTURE, &
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
                localError="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(localError,err,error,*999)
            END SELECT
          CASE DEFAULT
            localError="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Darcy pressure equation."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Darcy pressure equation."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Darcy pressure type of a fluid mechanics equation set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DarcyPressure_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("DarcyPressure_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Darcy pressure equation type of an fluid mechanics equations set class.
  SUBROUTINE DarcyPressure_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    
    ENTERS("DarcyPressure_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy pressure type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
        & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
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
        localError="The third equations set specification of "// &
          & TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Darcy pressure type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("DarcyPressure_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORSEXITS("DarcyPressure_EquationsSetSolutionMethodSet",err,error)
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets/changes the equation specification for a Darcy pressure type of a fluid mechanics equations set.
  SUBROUTINE DarcyPressure_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("DarcyPressure_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Darcy pressure type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_STATIC_INRIA_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRESSURE_HOLMES_MOW_SUBTYPE, &
          & EQUATIONS_SET_ELASTICITY_FLUID_PRES_HOLMES_MOW_ACTIVE_SUBTYPE)
         !ok
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(specification(3),"*",err,error))// &
          & " is not valid for a Darcy pressure type of a fluid mechanics equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_FLUID_MECHANICS_CLASS,EQUATIONS_SET_DARCY_PRESSURE_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("DarcyPressure_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("DarcyPressure_EquationsSetSpecificationSet",err,error)
    EXITS("DarcyPressure_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE DarcyPressure_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Darcy pressure problem post solve.
  SUBROUTINE DARCY_PRESSURE_POST_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DARCY_PRESSURE_POST_SOLVE",err,error,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            !nothing
          CASE DEFAULT
            localError="The third problem specification of "// &
              & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Darcy pressure type of a fluid mechanics problem."
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

    EXITS("DARCY_PRESSURE_POST_SOLVE")
    RETURN
999 ERRORSEXITS("DARCY_PRESSURE_POST_SOLVE",err,error)
    RETURN 1
    
  END SUBROUTINE DARCY_PRESSURE_POST_SOLVE

  !
  !================================================================================================================================
  !

  !>Darcy pressure problem pre solve.
  SUBROUTINE DARCY_PRESSURE_PRE_SOLVE(CONTROL_LOOP,SOLVER,err,error,*)

    !Argument variables
    TYPE(CONTROL_LOOP_TYPE), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER!<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError

    ENTERS("DARCY_PRESSURE_PRE_SOLVE",err,error,*999)
    
    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(ASSOCIATED(SOLVER)) THEN
        IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN 
          IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a Darcy pressure problem.",err,error,*999)
          END IF
          SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
          CASE(PROBLEM_STANDARD_ELASTICITY_FLUID_PRESSURE_SUBTYPE)
            !nothing
          CASE DEFAULT
            localError="The third problem specification of "// &
              & TRIM(NumberToVString(CONTROL_LOOP%PROBLEM%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a Darcy pressure type of a fluid mechanics problem."
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

    EXITS("DARCY_PRESSURE_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("DARCY_PRESSURE_PRE_SOLVE",err,error)
    RETURN 1
  END SUBROUTINE DARCY_PRESSURE_PRE_SOLVE

  !
  !================================================================================================================================
  !

END MODULE DARCY_PRESSURE_EQUATIONS_ROUTINES
!
