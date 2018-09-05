!> \file
!> \author Chris Bradley
!> \brief This module handles all Helmholtz equations routines.
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

!>This module handles all Helmholtz equations routines.
MODULE HELMHOLTZ_EQUATIONS_ROUTINES

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

!!MERGE: move

  PUBLIC Helmholtz_BoundaryConditionsAnalyticCalculate
  
  PUBLIC Helmholtz_EquationsSetSolutionMethodSet
  
  PUBLIC HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP

  PUBLIC Helmholtz_EquationsSetSpecificationSet
  
  PUBLIC HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE

  PUBLIC HELMHOLTZ_EQUATION_PROBLEM_SETUP

  PUBLIC Helmholtz_ProblemSpecificationSet

CONTAINS

  !
  !================================================================================================================================
  !


  !>Calculates the analytic solution and sets the boundary conditions for an analytic problem.
  SUBROUTINE Helmholtz_BoundaryConditionsAnalyticCalculate(EQUATIONS_SET,BOUNDARY_CONDITIONS,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET
    TYPE(BOUNDARY_CONDITIONS_TYPE), POINTER :: BOUNDARY_CONDITIONS
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,node_idx,NUMBER_OF_DIMENSIONS,variable_idx,variable_type,local_ny,deriv_idx,dim_idx
    REAL(DP), POINTER :: GEOMETRIC_PARAMETERS(:), MATERIALS_PARAMETERS(:)
    TYPE(DOMAIN_TYPE), POINTER :: DOMAIN
    TYPE(DOMAIN_NODES_TYPE), POINTER :: DOMAIN_NODES
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE,MATERIALS_VARIABLE
    REAL(DP) :: VALUE,X(3),k,mu
    TYPE(VARYING_STRING) :: LOCAL_ERROR    
    
    ENTERS("Helmholtz_BoundaryConditionsAnalyticCalculate",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
        DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
        IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
          GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
          IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN            
            MATERIALS_FIELD=>EQUATIONS_SET%MATERIALS%MATERIALS_FIELD
            IF(ASSOCIATED(MATERIALS_FIELD)) THEN            
              CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
              NULLIFY(GEOMETRIC_VARIABLE)
              NULLIFY(GEOMETRIC_PARAMETERS)
              CALL Field_VariableGet(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,GEOMETRIC_VARIABLE,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,GEOMETRIC_PARAMETERS, &
                & err,error,*999)
              NULLIFY(MATERIALS_VARIABLE)
              NULLIFY(MATERIALS_PARAMETERS)
              CALL Field_VariableGet(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,MATERIALS_VARIABLE,err,error,*999)
              CALL FIELD_PARAMETER_SET_DATA_GET(MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,MATERIALS_PARAMETERS, &
                & err,error,*999)
              IF(ASSOCIATED(BOUNDARY_CONDITIONS)) THEN
                DO variable_idx=1,DEPENDENT_FIELD%NUMBER_OF_VARIABLES
                  variable_type=DEPENDENT_FIELD%VARIABLES(variable_idx)%VARIABLE_TYPE
                  FIELD_VARIABLE=>DEPENDENT_FIELD%VARIABLE_TYPE_MAP(variable_type)%PTR
                  IF(ASSOCIATED(FIELD_VARIABLE)) THEN
                    CALL FIELD_PARAMETER_SET_CREATE(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE,err,error,*999)
                    DO component_idx=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                      IF(FIELD_VARIABLE%COMPONENTS(component_idx)%INTERPOLATION_TYPE==FIELD_NODE_BASED_INTERPOLATION) THEN
                        IF(MATERIALS_VARIABLE%COMPONENTS(1)%INTERPOLATION_TYPE==FIELD_CONSTANT_INTERPOLATION) THEN
                          local_ny=MATERIALS_VARIABLE%COMPONENTS(1)%PARAM_TO_DOF_MAP%CONSTANT_PARAM2DOF_MAP
                          k=MATERIALS_PARAMETERS(local_ny)
                          mu=k/SQRT(2.0_DP)
                          DOMAIN=>FIELD_VARIABLE%COMPONENTS(component_idx)%DOMAIN
                          IF(ASSOCIATED(DOMAIN)) THEN
                            IF(ASSOCIATED(DOMAIN%TOPOLOGY)) THEN
                              DOMAIN_NODES=>DOMAIN%TOPOLOGY%NODES
                              IF(ASSOCIATED(DOMAIN_NODES)) THEN
                                !Loop over the local nodes excluding the ghosts.
                                DO node_idx=1,DOMAIN_NODES%NUMBER_OF_NODES
                                  !!TODO \todo We should interpolate the geometric field here and the node position.
                                  DO dim_idx=1,NUMBER_OF_DIMENSIONS
                                    local_ny=GEOMETRIC_VARIABLE%COMPONENTS(dim_idx)%PARAM_TO_DOF_MAP%NODE_PARAM2DOF_MAP% &
                                      & NODES(node_idx)%DERIVATIVES(1)%VERSIONS(1)
                                    X(dim_idx)=GEOMETRIC_PARAMETERS(local_ny)
                                  ENDDO !dim_idx
                                  !Loop over the derivatives
                                  DO deriv_idx=1,DOMAIN_NODES%NODES(node_idx)%NUMBER_OF_DERIVATIVES
                                    SELECT CASE(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE)
                                    CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1)
                                      !u=cos(k*x/sqrt(2))*sin(k*y/sqrt(2))
                                      SELECT CASE(variable_type)
                                      CASE(FIELD_U_VARIABLE_TYPE)
                                        SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                        CASE(NO_GLOBAL_DERIV)
                                          VALUE=COS(mu*X(1))*SIN(mu*X(2))
                                        CASE(GLOBAL_DERIV_S1)
                                          VALUE=-mu*SIN(mu*X(1))*SIN(mu*X(2))
                                        CASE(GLOBAL_DERIV_S2)
                                          VALUE=mu*COS(mu*X(1))*COS(mu*X(2))
                                        CASE(GLOBAL_DERIV_S1_S2)
                                          VALUE=-mu*mu*SIN(mu*X(1))*COS(mu*X(2))
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                            & err,error))//" is invalid."
                                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                                        END SELECT
                                      CASE(FIELD_DELUDELN_VARIABLE_TYPE)
                                       SELECT CASE(DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX)
                                        CASE(NO_GLOBAL_DERIV)
                                          VALUE=0.0_DP !!TODO
                                        CASE(GLOBAL_DERIV_S1)
                                          CALL FlagError("Not implemented.",err,error,*999)
                                        CASE(GLOBAL_DERIV_S2)
                                          CALL FlagError("Not implemented.",err,error,*999)
                                        CASE(GLOBAL_DERIV_S1_S2)
                                          CALL FlagError("Not implemented.",err,error,*999)
                                        CASE DEFAULT
                                          LOCAL_ERROR="The global derivative index of "//TRIM(NUMBER_TO_VSTRING( &
                                            & DOMAIN_NODES%NODES(node_idx)%DERIVATIVES(deriv_idx)%GLOBAL_DERIVATIVE_INDEX,"*", &
                                            & err,error))//" is invalid."
                                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                                        END SELECT
                                      CASE DEFAULT
                                        LOCAL_ERROR="The variable type of "//TRIM(NUMBER_TO_VSTRING(variable_type,"*",err,error))//&
                                          & " is invalid."
                                        CALL FlagError(LOCAL_ERROR,err,error,*999)
                                      END SELECT
                                    CASE DEFAULT
                                      LOCAL_ERROR="The analytic function type of "// &
                                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                                        & " is invalid."
                                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                                    END SELECT
                                    local_ny=FIELD_VARIABLE%COMPONENTS(component_idx)%PARAM_TO_DOF_MAP% &
                                      & NODE_PARAM2DOF_MAP%NODES(node_idx)%DERIVATIVES(deriv_idx)%VERSIONS(1)
                                    CALL FIELD_PARAMETER_SET_UPDATE_LOCAL_DOF(DEPENDENT_FIELD,variable_type, &
                                      & FIELD_ANALYTIC_VALUES_SET_TYPE,local_ny,VALUE,err,error,*999)
                                    IF(variable_type==FIELD_U_VARIABLE_TYPE) THEN
                                      IF(DOMAIN_NODES%NODES(node_idx)%BOUNDARY_NODE) THEN
                                        !If we are a boundary node then set the analytic value on the boundary
                                        CALL BOUNDARY_CONDITIONS_SET_LOCAL_DOF(BOUNDARY_CONDITIONS,DEPENDENT_FIELD,variable_type, &
                                          & local_ny,BOUNDARY_CONDITION_FIXED,VALUE,err,error,*999)
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
                          CALL FlagError("Only constant interpolation of the material field is implemented.",err,error,*999)
                        ENDIF
                      ELSE
                        CALL FlagError("Only node based interpolation is implemented.",err,error,*999)
                      ENDIF
                    ENDDO !component_idx
                    CALL FIELD_PARAMETER_SET_UPDATE_START(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                      & err,error,*999)
                    CALL FIELD_PARAMETER_SET_UPDATE_FINISH(DEPENDENT_FIELD,variable_type,FIELD_ANALYTIC_VALUES_SET_TYPE, &
                      & err,error,*999)
                  ELSE
                    CALL FlagError("Field variable is not associated.",err,error,*999)
                  ENDIF

                ENDDO !variable_idx
                CALL FIELD_PARAMETER_SET_DATA_RESTORE(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
                  & GEOMETRIC_PARAMETERS,err,error,*999)
              ELSE
                CALL FlagError("Boundary conditions is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials field is not associated.",err,error,*999)
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
    
    EXITS("Helmholtz_BoundaryConditionsAnalyticCalculate")
    RETURN
999 ERRORS("Helmholtz_BoundaryConditionsAnalyticCalculate",err,error)
    EXITS("Helmholtz_BoundaryConditionsAnalyticCalculate")
    RETURN 1
    
  END SUBROUTINE Helmholtz_BoundaryConditionsAnalyticCalculate
  
  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Helmholtz equation finite element equations set.
  SUBROUTINE HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,mi,ms,nh,nhs,ni,ns
    REAL(DP) :: RWG,SUM,PGMSI(3),PGNSI(3),PGM,PGN,k
    TYPE(BASIS_TYPE), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingLinearType), POINTER :: linearMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesLinearType), POINTER :: linearMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: equationsMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FIELD_TYPE), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,MATERIALS_FIELD
    TYPE(FIELD_VARIABLE_TYPE), POINTER :: FIELD_VARIABLE
    TYPE(QUADRATURE_SCHEME_TYPE), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
#ifdef TAUPROF
    CHARACTER(26) :: CVAR
    INTEGER :: GAUSS_POINT_LOOP_PHASE(2) = [ 0, 0 ]
    SAVE GAUSS_POINT_LOOP_PHASE
#endif

    ENTERS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Helmholtz type equations set.", &
          & err,error,*999)
      END IF
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE)
!!TODO: move these and scale factor adjustment out once generalised Helmholtz is put in.
          !Store all these in equations matrices/somewhere else?????
          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          MATERIALS_FIELD=>equations%interpolation%materialsField
          vectorMatrices=>vectorEquations%vectorMatrices
          linearMatrices=>vectorMatrices%linearMatrices
          equationsMatrix=>linearMatrices%matrices(1)%PTR
          rhsVector=>vectorMatrices%rhsVector
          vectorMapping=>vectorEquations%vectorMapping
          linearMapping=>vectorMapping%linearMapping
          FIELD_VARIABLE=>linearMapping%equationsMatrixToVarMaps(1)%variable
          FIELD_VAR_TYPE=FIELD_VARIABLE%VARIABLE_TYPE
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%DECOMPOSITION%MESH_COMPONENT_NUMBER)%PTR% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%QUADRATURE_SCHEME_MAP(BASIS_DEFAULT_QUADRATURE_SCHEME)%PTR
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
          CALL FIELD_INTERPOLATION_PARAMETERS_ELEMENT_GET(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
          !Loop over gauss points
          DO ng=1,QUADRATURE_SCHEME%NUMBER_OF_GAUSS
#ifdef TAUPROF
              WRITE (CVAR,'(a17,i2)') 'Gauss Point Loop ',ng
              CALL TAU_PHASE_CREATE_DYNAMIC(GAUSS_POINT_LOOP_PHASE,CVAR)
              CALL TAU_PHASE_START(GAUSS_POINT_LOOP_PHASE)
#endif
            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
            CALL FIELD_INTERPOLATED_POINT_METRICS_CALCULATE(GEOMETRIC_BASIS%NUMBER_OF_XI,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)

            CALL FIELD_INTERPOLATE_GAUSS(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR,err,error,*999)
            k = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%PTR%VALUES(1,NO_PART_DERIV) ! wave number

            !Calculate RWG.
!!TODO: Think about symmetric problems. 
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR%JACOBIAN* &
              & QUADRATURE_SCHEME%GAUSS_WEIGHTS(ng)
            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
!!TODO: CHANGE ELEMENT CALCULATE TO WORK OF ns ???
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1
                nhs=0
                IF(equationsMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      PGM=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,NO_PART_DERIV,ng)
                      PGN=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,NO_PART_DERIV,ng)
                      DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        PGMSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                        PGNSI(ni)=QUADRATURE_SCHEME%GAUSS_BASIS_FNS(ns,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)
                      ENDDO !ni
                      SUM=0.0_DP
                      SUM=SUM-k*k*PGM*PGN
                      DO mi=1,DEPENDENT_BASIS%NUMBER_OF_XI
                        DO ni=1,DEPENDENT_BASIS%NUMBER_OF_XI
                          SUM=SUM+PGMSI(mi)*PGNSI(ni)*equations%interpolation% &
                            & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%PTR%GU(mi,ni)
                        ENDDO !ni
                      ENDDO !mi
                      equationsMatrix%elementMatrix%matrix(mhs,nhs)=equationsMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
              ENDDO !ms
            ENDDO !mh
#ifdef TAUPROF
            CALL TAU_PHASE_STOP(GAUSS_POINT_LOOP_PHASE)
#endif
          ENDDO !ng
          
          !Scale factor adjustment
          IF(DEPENDENT_FIELD%SCALINGS%SCALING_TYPE/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%PTR,err,error,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                mhs=mhs+1                    
                nhs=0
                IF(equationsMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%NUMBER_OF_COMPONENTS
                    DO ns=1,DEPENDENT_BASIS%NUMBER_OF_ELEMENT_PARAMETERS
                      nhs=nhs+1
                      equationsMatrix%elementMatrix%matrix(mhs,nhs)=equationsMatrix%elementMatrix%matrix(mhs,nhs)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)* &
                        & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ns,nh)
                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%PTR%SCALE_FACTORS(ms,mh)
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE(EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
            & " is not valid for a Helmholtz equation type of a classical field equations set class."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
        
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !

  !>Sets up the Helmholtz equation type of a classical field equations set class.
  SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Helmholtz equation on.
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Helmholtz type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE)
        CALL Helmholtz_EquationsSetStandardSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*999)
      CASE(EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Helmholtz equation type of a classical field equation set class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP",err,error)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Helmholtz equation type of an classical field equations set class.
  SUBROUTINE Helmholtz_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Helmholtz_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Helmholtz type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE)        
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
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)        
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
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",err,error))//" is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Helmholtz equation type of an classical field equations set class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Helmholtz_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Helmholtz_EquationsSetSolutionMethodSet",err,error)
    EXITS("Helmholtz_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a Helmholtz equation type of a classical field equations set class.
  SUBROUTINE Helmholtz_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: subtype

    ENTERS("Helmholtz_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Helmholtz type equations set.", &
          & err,error,*999)
      END IF
      subtype=specification(3)
      SELECT CASE(subtype)
      CASE(EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE)
        !ok
      CASE(EQUATIONS_SET_GENERALISED_HELMHOLTZ_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        localError="The third equations set specification of "//TRIM(NumberToVstring(subtype,"*",err,error))// &
          & " is not valid for a Helmholtz type of a classical field equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_CLASSICAL_FIELD_CLASS,EQUATIONS_SET_HELMHOLTZ_EQUATION_TYPE,subtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Helmholtz_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Helmholtz_EquationsSetSpecificationSet",err,error)
    EXITS("Helmholtz_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Helmholtz equation.
  SUBROUTINE Helmholtz_EquationsSetStandardSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EQUATIONS_SET_TYPE), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EQUATIONS_SET_SETUP_TYPE), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,NUMBER_OF_DIMENSIONS
    TYPE(DECOMPOSITION_TYPE), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(FIELD_TYPE), POINTER :: ANALYTIC_FIELD,DEPENDENT_FIELD,GEOMETRIC_FIELD
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EQUATIONS_SET_MATERIALS_TYPE), POINTER :: EQUATIONS_MATERIALS
    INTEGER(INTG) :: MATERIAL_FIELD_NUMBER_OF_VARIABLES, MATERIAL_FIELD_NUMBER_OF_COMPONENTS
    INTEGER(INTG) :: i,GEOMETRIC_COMPONENT_NUMBER
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Helmholtz_EquationsSetStandardSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)
   
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a Helmholtz type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_STANDARD_HELMHOLTZ_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%SETUP_TYPE)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Helmholtz_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD,err,error, &
              & *999)               
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
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
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,"Phi",err,error,*999)
              CALL FIELD_VARIABLE_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,"del Phi/del n", &
                & err,error,*999)
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
              CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,1,"Phi",err,error,*999)
              CALL FIELD_COMPONENT_LABEL_SET(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & "del Phi/del n",err,error,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL FIELD_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL FIELD_DEPENDENT_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_VARIABLES_CHECK(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
              CALL FIELD_VARIABLE_TYPES_CHECK(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                & err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL FIELD_DIMENSION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_DATA_TYPE_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
              SELECT CASE(EQUATIONS_SET%SOLUTION_METHOD)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL FIELD_COMPONENT_INTERPOLATION_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD_AUTO_CREATED) THEN
              CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          MATERIAL_FIELD_NUMBER_OF_VARIABLES=1
          MATERIAL_FIELD_NUMBER_OF_COMPONENTS=1 ! wave number k
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
                CALL FIELD_NUMBER_OF_VARIABLES_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD, &
                  & MATERIAL_FIELD_NUMBER_OF_VARIABLES,err,error,*999)
                CALL FIELD_VARIABLE_TYPES_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,[FIELD_U_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL FIELD_DIMENSION_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
                CALL FIELD_DATA_TYPE_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL FIELD_NUMBER_OF_COMPONENTS_SET_AND_LOCK(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & MATERIAL_FIELD_NUMBER_OF_COMPONENTS,err,error,*999)
                CALL FIELD_COMPONENT_MESH_COMPONENT_GET(EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE, &
                  & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)                

                !Default to constant interpolation of material field
                DO i=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_MESH_COMPONENT_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & i,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                  CALL FIELD_COMPONENT_INTERPOLATION_SET(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & i,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                END DO

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
                CALL FIELD_NUMBER_OF_COMPONENTS_CHECK(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF( ASSOCIATED(EQUATIONS_MATERIALS) ) THEN
              IF( EQUATIONS_MATERIALS%MATERIALS_FIELD_AUTO_CREATED ) THEN
                CALL FIELD_CREATE_FINISH(EQUATIONS_MATERIALS%MATERIALS_FIELD,err,error,*999)
                !Set the default values for the materials field
                DO i=1,MATERIAL_FIELD_NUMBER_OF_COMPONENTS
                  CALL FIELD_COMPONENT_VALUES_INITIALISE(EQUATIONS_MATERIALS%MATERIALS_FIELD,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE, i, 1.0_DP, ERR, ERROR, *999)
                ENDDO
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              DEPENDENT_FIELD=>EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD
              IF(ASSOCIATED(DEPENDENT_FIELD)) THEN
                GEOMETRIC_FIELD=>EQUATIONS_SET%GEOMETRY%GEOMETRIC_FIELD
                IF(ASSOCIATED(GEOMETRIC_FIELD)) THEN
                  CALL FIELD_NUMBER_OF_COMPONENTS_GET(GEOMETRIC_FIELD,FIELD_U_VARIABLE_TYPE,NUMBER_OF_DIMENSIONS,err,error,*999)
                  SELECT CASE(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE)
                  CASE(EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1)
                    !Check that we are in 2D
                    IF(NUMBER_OF_DIMENSIONS/=2) THEN
                      LOCAL_ERROR="The number of geometric dimensions of "// &
                        & TRIM(NUMBER_TO_VSTRING(NUMBER_OF_DIMENSIONS,"*",err,error))// &
                        & " is invalid. The analytic function type of "// &
                        & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                        & " requires that there be 2 geometric dimensions."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
                    ENDIF
                    !Create analytic field if required
                    !Set analtyic function type
                    EQUATIONS_SET%ANALYTIC%ANALYTIC_FUNCTION_TYPE=EQUATIONS_SET_HELMHOLTZ_EQUATION_TWO_DIM_1
                  CASE DEFAULT
                    LOCAL_ERROR="The specified analytic function type of "// &
                      & TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ANALYTIC_FUNCTION_TYPE,"*",err,error))// &
                      & " is invalid for a standard Helmholtz equation."
                    CALL FlagError(LOCAL_ERROR,err,error,*999)
                  END SELECT
                ELSE
                  CALL FlagError("Equations set geometric field is not associated.",err,error,*999)
                ENDIF
             ELSE
                CALL FlagError("Equations set dependent field is not associated.",err,error,*999)
              ENDIF
            ELSE
              CALL FlagError("Equations set dependent field has not been finished.",err,error,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(ASSOCIATED(EQUATIONS_SET%ANALYTIC)) THEN
              ANALYTIC_FIELD=>EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD
              IF(ASSOCIATED(ANALYTIC_FIELD)) THEN
                IF(EQUATIONS_SET%ANALYTIC%ANALYTIC_FIELD_AUTO_CREATED) THEN
                  CALL FIELD_CREATE_FINISH(EQUATIONS_SET%DEPENDENT%DEPENDENT_FIELD,err,error,*999)
                ENDIF
              ENDIF
            ELSE
              CALL FlagError("Equations set analytic is not associated.",err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%ACTION_TYPE)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%DEPENDENT_FINISHED) THEN
              CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
              CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
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
              CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,err,error,*999)
              CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE],err,error,*999)
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE],err,error,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE],err,error,*999)
                CALL EquationsMatrices_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%sparsityType,"*",err,error))//" is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SOLUTION_METHOD,"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Helmholtz equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a standard Helmholtz equation subtype."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Helmholtz_EquationsSetStandardSetup")
    RETURN
999 ERRORSEXITS("Helmholtz_EquationsSetStandardSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_EquationsSetStandardSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the Helmholtz problem.
  SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: PROBLEM !<A pointer to the problem set to setup a Helmholtz equation on.
    TYPE(PROBLEM_SETUP_TYPE), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Helmholtz equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE)
        CALL HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*999)
      CASE(PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE)
        CALL FlagError("Not implemented.",err,error,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " is not valid for a Helmholtz equation type of a classical field problem class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a Helmholtz equation type.
  SUBROUTINE Helmholtz_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(PROBLEM_TYPE), POINTER :: problem !<A pointer to the problem to set the problem specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemSubtype

    ENTERS("Helmholtz_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemSubtype)
        CASE(PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE)
          !ok
        CASE(PROBLEM_GENERALISED_HELMHOLTZ_SUBTYPE)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
            & " is not valid for a Helmholtz equation type of a classical field problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_CLASSICAL_FIELD_CLASS,PROBLEM_HELMHOLTZ_EQUATION_TYPE,problemSubtype]
      ELSE
        CALL FlagError("Helmholtz equation problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("Helmholtz_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Helmholtz_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Helmholtz_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets up the standard Helmholtz equations problem.
  SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

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
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP",err,error,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a Helmholtz equation problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_STANDARD_HELMHOLTZ_SUBTYPE) THEN
        SELECT CASE(PROBLEM_SETUP%SETUP_TYPE)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%ACTION_TYPE)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
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
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
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
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
            !Set solver defaults
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
                & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
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
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%ACTION_TYPE,"*",err,error))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
              & " is invalid for a standard Helmholtz equation."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
       CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%SETUP_TYPE,"*",err,error))// &
            & " is invalid for a standard Helmholtz equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",err,error))// &
          & " does not equal a standard Helmholtz equation subtype."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
       
    EXITS("HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP")
    RETURN
999 ERRORSEXITS("HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP",err,error)
    RETURN 1
  END SUBROUTINE HELMHOLTZ_EQUATION_PROBLEM_STANDARD_SETUP

  !
  !================================================================================================================================
  !
 
END MODULE HELMHOLTZ_EQUATIONS_ROUTINES
