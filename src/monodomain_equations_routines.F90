!> \file
!> \author Ishani Roy & Sander Land
!> \brief This module handles all Monodomain equations routines which use Strang Splitting with hard-coded cell models.
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

!>This module handles all Monodomain equations routines.
MODULE MONODOMAIN_EQUATIONS_ROUTINES

  USE BaseRoutines
  USE BasisRoutines
  USE BasisAccessRoutines
  USE BOUNDARY_CONDITIONS_ROUTINES
  USE Constants
  USE ControlLoopRoutines
  USE ControlLoopAccessRoutines
  USE DistributedMatrixVector
  USE DistributedMatrixVectorAccessRoutines
  USE DomainMappings
  USE ELECTROPHYSIOLOGY_CELL_ROUTINES
  USE EquationsRoutines
  USE EquationsAccessRoutines
  USE EquationsMappingRoutines
  USE EquationsMatricesRoutines
  USE EquationsSetConstants
  USE EquationsSetAccessRoutines
  USE FIELD_IO_ROUTINES
  USE FieldRoutines
  USE FieldAccessRoutines
  USE InputOutput
  USE ISO_VARYING_STRING
  USE Kinds
  USE MatrixVector
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

  PUBLIC MONODOMAIN_CONTROL_LOOP_POST_LOOP
 
  PUBLIC MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP
  
  PUBLIC Monodomain_FiniteElementCalculate

  PUBLIC Monodomain_EquationsSetSolutionMethodSet
  
  PUBLIC Monodomain_EquationsSetSpecificationSet

  PUBLIC Monodomain_ProblemSpecificationSet
  
  PUBLIC MONODOMAIN_EQUATION_PROBLEM_SETUP
  
  PUBLIC MONODOMAIN_PRE_SOLVE,MONODOMAIN_POST_SOLVE
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE MONODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: equations_set_idx
    TYPE(ControlLoopTimeType), POINTER :: TIME_LOOP,TIME_LOOP_PARENT
    TYPE(ControlLoopType), POINTER :: PARENT_LOOP
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD
    TYPE(ProblemType), POINTER :: PROBLEM
    TYPE(RegionType), POINTER :: DEPENDENT_REGION   
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER,CURRENT_LOOP_ITERATION

    ENTERS("MONODOMAIN_CONTROL_LOOP_POST_LOOP",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP)) THEN
      IF(CONTROL_LOOP%outputType>=CONTROL_LOOP_PROGRESS_OUTPUT) THEN
        SELECT CASE(CONTROL_LOOP%loopType)
        CASE(CONTROL_SIMPLE_TYPE)
          !do nothing
        CASE(CONTROL_FIXED_LOOP_TYPE)
          !do nothing
        CASE(CONTROL_TIME_LOOP_TYPE)
          !Export the dependent field for this time step
          TIME_LOOP=>CONTROL_LOOP%timeLoop
          IF(ASSOCIATED(TIME_LOOP)) THEN
            PROBLEM=>CONTROL_LOOP%PROBLEM
            IF(ASSOCIATED(PROBLEM)) THEN
              NULLIFY(SOLVERS)
              NULLIFY(SOLVER)
              !Get the solver. 
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)            
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
              !Loop over the equations sets associated with the solver
              SOLVER_EQUATIONS=>SOLVER%SOLVER_EQUATIONS
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%dependent%dependentField
                      NULLIFY(DEPENDENT_REGION)
                      CALL Field_RegionGet(DEPENDENT_FIELD,DEPENDENT_REGION,ERR,ERROR,*999)
                      NULLIFY(PARENT_LOOP)
                      PARENT_LOOP=>CONTROL_LOOP%parentLoop
                      IF(ASSOCIATED(PARENT_LOOP)) THEN
                        !add the iteration number of the parent loop to the filename
                        NULLIFY(TIME_LOOP_PARENT)
                        TIME_LOOP_PARENT=>PARENT_LOOP%timeLoop
                        IF(ASSOCIATED(TIME_LOOP_PARENT)) THEN
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP_PARENT%outputNumber
                          CURRENT_LOOP_ITERATION=TIME_LOOP_PARENT%globalIterationNumber
                          FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%userNumber,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP_PARENT%globalIterationNumber,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%iterationNumber,"*",ERR,ERROR))
                        ELSE
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP%outputNumber
                          CURRENT_LOOP_ITERATION=TIME_LOOP%globalIterationNumber
                          FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%userNumber,"*",ERR,ERROR))// &
                            & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%globalIterationNumber,"*",ERR,ERROR))
                        ENDIF
                      ELSE
                        OUTPUT_ITERATION_NUMBER=TIME_LOOP%outputNumber
                        CURRENT_LOOP_ITERATION=TIME_LOOP%globalIterationNumber
                        FILENAME="Time_"//TRIM(NUMBER_TO_VSTRING(DEPENDENT_REGION%userNumber,"*",ERR,ERROR))// &
                          & "_"//TRIM(NUMBER_TO_VSTRING(TIME_LOOP%globalIterationNumber,"*",ERR,ERROR))
                      ENDIF
                      METHOD="FORTRAN"
                      IF(OUTPUT_ITERATION_NUMBER/=0.AND.MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0) THEN
                        CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,ERR,ERROR,*999)
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NUMBER_TO_VSTRING(equations_set_idx,"*",ERR,ERROR))// &
                        & " in the solver mapping."
                      CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
                    ENDIF
                  ENDDO !equations_set_idx
                ELSE
                  CALL FlagError("Solver equations solver mapping is not associated.",ERR,ERROR,*999)
                ENDIF
              ELSE
                CALL FlagError("Solver solver equations are not associated.",ERR,ERROR,*999)
              ENDIF
            ELSE
              CALL FlagError("Control loop problem is not associated.",ERR,ERROR,*999)
            ENDIF
          ELSE
            CALL FlagError("Time loop is not associated.",ERR,ERROR,*999)
          ENDIF
        CASE(CONTROL_WHILE_LOOP_TYPE)
          !do nothing
        CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE DEFAULT
          LOCAL_ERROR="The control loop type of "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%loopType,"*",ERR,ERROR))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",ERR,ERROR,*999)
    ENDIF

    EXITS("MONODOMAIN_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_CONTROL_LOOP_POST_LOOP",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE MONODOMAIN_CONTROL_LOOP_POST_LOOP
  !
  !================================================================================================================================
  !
  !

  !>Sets the problem specification for a monodomain equation set class.
  SUBROUTINE Monodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    INTEGER(INTG) :: subtype
    TYPE(VARYING_STRING) :: localError

    ENTERS("Monodomain_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)<3) THEN
        CALL FlagError("Equations set specification must have at least three entries for a monodomain class equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(specification(2))
      CASE(EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        subtype=specification(3)
        SELECT CASE(subtype)
        CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="Equations set subtype "//TRIM(NumberToVstring(subtype,"*",err,error))// &
            & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        !Set full specification
        IF(ALLOCATED(equationsSet%specification)) THEN
          CALL FlagError("Equations set specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(equationsSet%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
        END IF
        equationsSet%specification(1:3)=[EQUATIONS_SET_BIOELECTRICS_CLASS, &
          & EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE,subtype]
      CASE DEFAULT
        localError="Equations set equation type "//TRIM(NumberToVstring(specification(2),"*",err,error))// &
          & " is not valid for a monodomain equations set class."
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated",err,error,*999)
    END IF

    EXITS("Monodomain_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Monodomain_EquationsSetSpecificationSet",err,error)
    EXITS("Monodomain_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a monodomain problem class.
  SUBROUTINE Monodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType,problemSubtype

    ENTERS("Monodomain_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)>=3) THEN
        problemType=problemSpecification(2)
        SELECT CASE(problemType)
        CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
          problemSubtype=problemSpecification(3)
          SELECT CASE(problemSubtype)
          CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE, &
            & PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
            !ok
          CASE DEFAULT
            localError="Problem subtype "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
              & " is not valid for a Monodomain equation type of a Strang splitting  problem class."
            CALL FlagError(localError,err,error,*999)
          END SELECT
          IF(ALLOCATED(problem%specification)) THEN
            CALL FlagError("Problem specification is already allocated.",err,error,*999)
          ELSE
            ALLOCATE(problem%specification(3),stat=err)
            IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
          END IF
          problem%specification(1:3)=[PROBLEM_BIOELECTRICS_CLASS,PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE,problemSubtype]
        CASE DEFAULT
          localError="Problem equation type "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a monodomain problem class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Monodomain problem specification must have a type.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated",err,error,*999)
    END IF

    EXITS("Monodomain_ProblemSpecificationSet")
    RETURN
999 ERRORSEXITS("Monodomain_ProblemSpecificationSet",err,error)
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a Monodomain equation finite element equations set.
  SUBROUTINE Monodomain_FiniteElementCalculate(EQUATIONS_SET,ELEMENT_NUMBER,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,ng,mh,mhs,ms,nh,nhs,ni,ns,nj
    REAL(DP) :: RWG,SUM,Df, Dt, D(3,3), f(3), fnorm
    REAL(DP) :: DPHIDX(3,8) ! assumes <= 8 basis functions / DEPENDENT_BASIS%numberOfElementParameters <= 8
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,materialsField
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QuadratureSchemeType), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
     
    ENTERS("Monodomain_FiniteElementCalculate",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
          CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
            & err,error,*999)
        END IF
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
        CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)

          DEPENDENT_FIELD=>equations%interpolation%dependentField
          GEOMETRIC_FIELD=>equations%interpolation%geometricField
          materialsField=>equations%interpolation%materialsField
          vectorEquations=>equations%vectorEquations
          vectorMatrices=>vectorEquations%vectorMatrices
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          rhsVector=>vectorMatrices%rhsVector

          IF(.NOT.(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix.OR.rhsVector%updateVector)) RETURN ! no updates -> return immediately
 
          ! set up opencmiss interpolation and basis
          vectorMapping=>vectorEquations%vectorMapping
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%VARIABLE
          FIELD_VAR_TYPE=FIELD_VARIABLE%variableType
          GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr
          DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%decomposition%meshComponentNumber)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%decomposition%meshComponentNumber)%ptr% &
            & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
          QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr

          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)
          CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
            & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)


          DO ng=1,QUADRATURE_SCHEME%numberOfGauss
            ! get interpolated geometric and material interpolated point
            CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)
            CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%numberOfXi,equations%interpolation% &
              & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)
            CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
              & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,ERR,ERROR,*999)

            ! calculate weight = det J * gauss pt weight
            RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
              & QUADRATURE_SCHEME%gaussWeights(ng)

            ! basis function chain rule taken out of the inner loop.
            DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
            DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
              DPHIDX(nj,ms)=0.0_DP
              DO ni=1,DEPENDENT_BASIS%numberOfXi                          
                DPHIDX(nj,ms)=DPHIDX(nj,ms) + &
                             & QUADRATURE_SCHEME%gaussBasisFunctions(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                             & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
              ENDDO !ni
            ENDDO !nj
            ENDDO !ms

            ! Diffusion tensor  D  =  Dt I + (Df - Dt) f f^T  where Dt and Df are diffusivity/conductivity in fiber/transverse directions
            Df = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,NO_PART_DERIV) ! 2 = Df
            Dt = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3,NO_PART_DERIV) ! 3 = Dt
            fnorm = 0.0
            DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
              f(nj) = equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(3+nj,NO_PART_DERIV) ! 4,5[,6] = f
              fnorm = fnorm + f(nj)*f(nj)
            ENDDO
            ! normalize f, and fill in default for 0,0,0 -> 1,0,0
            fnorm = SQRT(fnorm)
            IF(fnorm < 1e-6) THEN
              f = [ 1.0, 0.0, 0.0 ] ! default
            ELSE
              f = f / fnorm
            ENDIF
            DO ni=1,GEOMETRIC_VARIABLE%numberOfComponents
              D(ni,:)  = 0.0
              D(ni,ni) = Dt
              DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                D(ni,nj) = D(ni,nj) + (Df - Dt) * f(ni) * f(nj)
              ENDDO
            ENDDO

            !Loop over field components
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns. TODO: use symmetry?
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        SUM=0.0_DP
                        DO ni=1,GEOMETRIC_VARIABLE%numberOfComponents 
                        DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                          SUM=SUM + D(ni,nj) * DPHIDX(ni,ms) * DPHIDX(nj,ns)
                        ENDDO !nj
                        ENDDO !ni
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG ! Aij = int D_ij * dphi_m/dx_i * dphi_n/dx_j 
                      ENDIF
 
                      IF(dampingMatrix%updateMatrix) THEN ! non mass lumped version
                         dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                            & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                            & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG ! int phi_m phi_n
                      ENDIF 

                    ENDDO !ns
                  ENDDO !nh
                ENDIF
                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
!                IF(dampingMatrix%updateMatrix) THEN  ! mass lumnped version
!                  dampingMatrix%elementMatrix%matrix(mhs,mhs)=dampingMatrix%elementMatrix%matrix(mhs,mhs)+ &
!                  & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)*RWG   !  // int phi_m
!                ENDIF
              ENDDO !ms
            ENDDO !mh
          IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP 
          ENDDO !ng
          IF(DEPENDENT_FIELD%SCALINGS%scalingType/=FIELD_NO_SCALING) THEN
            CALL Field_InterpolationParametersScaleFactorsElementGet(ELEMENT_NUMBER,equations%interpolation% &
              & dependentInterpParameters(FIELD_VAR_TYPE)%ptr,ERR,ERROR,*999)
            mhs=0          
            DO mh=1,FIELD_VARIABLE%numberOfComponents
              !Loop over element rows
              DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                mhs=mhs+1                    
                nhs=0
                IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix) THEN
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      IF(stiffnessMatrix%updateMatrix) THEN
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                          & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN ! non mass lumped version
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
                           & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)* &
                           & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ns,nh)
                      ENDIF 
                    ENDDO !ns
                  ENDDO !nh
                ENDIF

                IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=rhsVector%elementVector%vector(mhs)* &
                  & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh) 
 !               IF(dampingMatrix%updateMatrix) THEN ! mass lumped version
 !                  dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)* &
 !                    & equations%interpolation%dependentInterpParameters(FIELD_VAR_TYPE)%ptr%scaleFactors(ms,mh)
 !               ENDIF
              ENDDO !ms
            ENDDO !mh
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
            & " is not valid for a Monodomain equation type of a Strang splitting equations set class."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
 
    EXITS("Monodomain_FiniteElementCalculate")
    RETURN
999 ERRORSEXITS("Monodomain_FiniteElementCalculate",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Monodomain_FiniteElementCalculate

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation type of a Strang splitting equations set class.
  SUBROUTINE MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a Monodomain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        CALL Monodomain_EquationsSetSubtypeSetUP(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*999)        
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of a Strang splitting equation set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_EQUATIONS_SET_SETUP

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a Monodomain equation type of an Strang splitting equations set class.
  SUBROUTINE Monodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Monodomain_EquationsSetSolutionMethodSet",ERR,ERROR,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
      CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE,EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)        
        SELECT CASE(SOLUTION_METHOD)
        CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
          EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
        CASE(EQUATIONS_SET_BEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FD_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_FV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFEM_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE(EQUATIONS_SET_GFV_SOLUTION_METHOD)
          CALL FlagError("Not implemented.",ERR,ERROR,*999)
        CASE DEFAULT
          LOCAL_ERROR="The specified solution method of "//TRIM(NUMBER_TO_VSTRING(SOLUTION_METHOD,"*",ERR,ERROR))//" is invalid."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of monodomain Strang splitting equations set class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Monodomain_EquationsSetSolutionMethodSet",ERR,ERROR)
    EXITS("Monodomain_EquationsSetSolutionMethodSet")
    RETURN 1

  END SUBROUTINE Monodomain_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets up the Monodomain equation.
  SUBROUTINE Monodomain_EquationsSetSubtypeSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,GEOMETRIC_MESH_COMPONENT,GEOMETRIC_SCALING_TYPE,numberOfDimensions, &
      & NUMBER_OF_MATERIALS_COMPONENTS, NUM_COMP
    TYPE(DecompositionType), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Monodomain_EquationsSetSubtypeSetup",ERR,ERROR,*999)
 
    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)


    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a monodomain type equations set.", &
          & err,error,*999)
      END IF
      IF(EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE .OR. &
          & EQUATIONS_SET%SPECIFICATION(3)==EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE) THEN
        SELECT CASE(EQUATIONS_SET_SETUP%setupType)
        CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL Monodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
              & ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & dependentField,ERR,ERROR,*999)
              CALL Field_TypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DEPENDENT_TYPE,ERR,ERROR,*999)
              CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL Field_DecompositionSetAndLock(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%dependent%dependentField,EQUATIONS_SET%GEOMETRY% &
                & geometricField,ERR,ERROR,*999)
              CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%dependent%dependentField,2,ERR,ERROR,*999)
              CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                & FIELD_DELUDELN_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999) 

          

              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1,&
                & ERR,ERROR,*999)
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,&
                & 1 ,ERR,ERROR,*999)

              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%solutionMethod,"*",ERR,ERROR))// &
                  & " is invalid or not implemented"
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
            ELSE
              ! user specified field
              CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
            ENDIF
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
              CALL Field_CreateFinish(EQUATIONS_SET%dependent%dependentField,ERR,ERROR,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation"
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%independentFieldAutoCreated) THEN
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
                 NUM_COMP = 4
              CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
                 NUM_COMP = 19
              CASE DEFAULT
                CALL FlagError("Invalid cell model equations set subtype",ERR,ERROR,*999)
              END SELECT    

              CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                & independentField,ERR,ERROR,*999)
              CALL Field_TypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_GENERAL_TYPE,ERR,ERROR,*999)
              CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,&
                   & FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
              CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
              CALL Field_DecompositionSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_DECOMPOSITION, &
                & ERR,ERROR,*999)
              CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,EQUATIONS_SET%GEOMETRY% &
                & geometricField,ERR,ERROR,*999)
              CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,1,ERR,ERROR,*999)
              CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,&
                   &[FIELD_U_VARIABLE_TYPE],ERR,ERROR,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,ERR,ERROR,*999)
         
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,&
                   &NUM_COMP,ERR,ERROR,*999)
              CALL Field_DOFOrderTypeSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,&
              &  FIELD_CONTIGUOUS_COMPONENT_DOF_ORDER,ERR,ERROR,*999) ! dofs continuous, so first + (x-1) is x'th component index

              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)

              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              END SELECT
            ELSE
              ! user specified field
              CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
            ENDIF

          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            IF(EQUATIONS_SET%INDEPENDENT%independentFieldAutoCreated) THEN
              CALL Field_CreateFinish(EQUATIONS_SET%INDEPENDENT%independentField,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
              CASE(EQUATIONS_SET_MONODOMAIN_BUENOOROVIO_SUBTYPE)
                CALL BUENO_OROVIO_INITIALIZE(EQUATIONS_SET%INDEPENDENT%independentField,ERR,ERROR,*999) ! initialize to y0
              CASE(EQUATIONS_SET_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
                CALL TENTUSSCHER06_INITIALIZE(EQUATIONS_SET%INDEPENDENT%independentField,ERR,ERROR,*999) ! initialize to y0 
              CASE DEFAULT
                CALL FlagError("Invalid cell model equations set subtype",ERR,ERROR,*999)
              END SELECT              
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation"
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT

        CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
           CASE(EQUATIONS_SET_SETUP_START_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
             IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
               IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
                !Create the auto created materials field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
                  & materialsField,ERR,ERROR,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_MATERIAL_TYPE,ERR,ERROR,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_INDEPENDENT_TYPE,ERR,ERROR,*999)
                CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,ERR,ERROR,*999)
                CALL Field_DecompositionSetAndLock(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_DECOMPOSITION, &
                  & ERR,ERROR,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%materialsField,EQUATIONS_SET%GEOMETRY% &
                  & geometricField,ERR,ERROR,*999)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_MATERIALS%materialsField,1,ERR,ERROR,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%materialsField,[FIELD_U_VARIABLE_TYPE], &
                  & ERR,ERROR,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_VECTOR_DIMENSION_TYPE,ERR,ERROR,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,ERR,ERROR,*999)
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,ERR,ERROR,*999)
                  !Materials field components are
                  ! 1. activation  factor (usually 0.0 or 1.0)
                  ! 2,3 for fiber/transverse conductivity   . defaults to constant interpolation 
                  ! 4,5[,6] : fiber unit vector in dimension
                  ! 7: out - activation times
                NUMBER_OF_MATERIALS_COMPONENTS= 7 !numberOfDimensions + 3
                 !Set the number of materials components
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & NUMBER_OF_MATERIALS_COMPONENTS,ERR,ERROR,*999)

                ! 1st = activation = node based, 2 3 diffusion constants
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 1,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 2,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 3,FIELD_CONSTANT_INTERPOLATION,ERR,ERROR,*999)
                ! 4 5 (6) fiber unit vector
                DO component_idx=1,3 !numberOfDimensions
                  CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                    & component_idx,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & component_idx+3,GEOMETRIC_MESH_COMPONENT,ERR,ERROR,*999)
                  CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & component_idx+3,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)
                ENDDO !component_idx

                  CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & NUMBER_OF_MATERIALS_COMPONENTS,FIELD_NODE_BASED_INTERPOLATION,ERR,ERROR,*999)

                  !Default the field scaling to that of the geometric field
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_SCALING_TYPE,ERR,ERROR,*999)
              ELSE
                ! user specified field
                CALL FlagError("No user specified field supported!",ERR,ERROR,*999)
              ENDIF
             ELSE

               CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
             ENDIF
           CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
            IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
              IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
                !Finish creating the materials field
                CALL Field_CreateFinish(EQUATIONS_MATERIALS%materialsField,ERR,ERROR,*999)
                !Set the default values for the materials field
                CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                  & numberOfDimensions,ERR,ERROR,*999)
                !Materials field components are 1 for each dimension 
                !i.e., k in div(k.grad(u(x)))
                NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions                             
                !First set the k values to 1.0
                DO component_idx=1,numberOfDimensions
                  CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,component_idx,1.0_DP,ERR,ERROR,*999)
                ENDDO !component_idx
              ENDIF
            ELSE
              CALL FlagError("Equations set materials is not associated.",ERR,ERROR,*999)
            ENDIF
! ! ! Upto here
           CASE DEFAULT
             LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
               & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
               & " is invalid for a Monodomain equation."
             CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
           END SELECT
        CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            !Do nothing
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            !Do nothing
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
          CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
          SELECT CASE(EQUATIONS_SET_SETUP%actionType)
          CASE(EQUATIONS_SET_SETUP_START_ACTION)
            CALL EquationsSet_AssertDependentIsFinished(EQUATIONS_SET,err,error,*999)
            CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
            CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
          CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
            SELECT CASE(EQUATIONS_SET%solutionMethod)
            CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
              !Finish the equations creation
              CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,ERR,ERROR,*999)
              CALL Equations_CreateFinish(equations,ERR,ERROR,*999)
              NULLIFY(vectorEquations)
              CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
              !Create the equations mapping.
              CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,ERR,ERROR,*999)
              !!! Check this 
              !CALL EquationsMapping_LinearMatricesNumberSet(vectorMapping,1,ERR,ERROR,*999)
              CALL EquationsMapping_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,ERR,ERROR,*999)
              !CALL EquationsMapping_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_U_VARIABLE_TYPE], &
               ! & ERR,ERROR,*999)
              CALL EquationsMapping_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,ERR,ERROR,*999)
              !!!! Check the above two lines
              CALL EquationsMapping_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,ERR,ERROR,*999)
              CALL EquationsMapping_VectorCreateFinish(vectorMapping,ERR,ERROR,*999)
              !Create the equations matrices
              CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,ERR,ERROR,*999)
              SELECT CASE(EQUATIONS%sparsityType)
              CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                CALL EquationsMatrices_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                  & ERR,ERROR,*999)
              CASE(EQUATIONS_MATRICES_SPARSE_MATRICES) 
                !CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                 ! & ERR,ERROR,*999)
                !CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                 ! & ERR,ERROR,*999)
                 CALL EquationsMatrices_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & ERR,ERROR,*999)
                  CALL EquationsMatrices_DynamicStructureTypeSet(vectorMatrices, &
                    [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],ERR,ERROR,*999)     
              CASE DEFAULT
                LOCAL_ERROR="The equations matrices sparsity type of "// &
                  & TRIM(NUMBER_TO_VSTRING(EQUATIONS%sparsityType,"*",ERR,ERROR))//" is invalid."
                CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
              END SELECT
              CALL EquationsMatrices_VectorCreateFinish(vectorMatrices,ERR,ERROR,*999)
            CASE DEFAULT
                LOCAL_ERROR="The solution method of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%solutionMethod,"*",ERR,ERROR))// &
                & " is invalid or not implemented."
              CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a Monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET_SETUP%setupType,"*",ERR,ERROR))// &
            & " is invalid for a Monodomain equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The equations set subtype of "//TRIM(NUMBER_TO_VSTRING(EQUATIONS_SET%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Monodomain equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_EquationsSetSubtypeSetup")
    RETURN
999 ERRORSEXITS("Monodomain_EquationsSetSubtypeSetup",ERR,ERROR)
    RETURN 1
    
  END SUBROUTINE Monodomain_EquationsSetSubtypeSetup

  !
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must have at least two entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        CALL Monodomain_ProblemStrangSplittingSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation Strang splitting problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_PROBLEM_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SETUP
  
!
  !================================================================================================================================
  !
 
  !>Sets up the Monodomain solution.
  SUBROUTINE Monodomain_ProblemStrangSplittingSetup(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)
    
    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the solutions set to setup a Monodomain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("Monodomain_ProblemStrangSplittingSetup",ERR,ERROR,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain Strang splitting problem.",err,error,*999)
      END IF
      SELECT CASE(PROBLEM%SPECIFICATION(3))
      CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE,PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
        CALL MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*999)
      CASE DEFAULT
        LOCAL_ERROR="Problem subtype "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " is not valid for a Monodomain equation type of a Strang splitting problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN
999 ERRORS("Monodomain_ProblemStrangSplittingSetup",ERR,ERROR)
    EXITS("Monodomain_ProblemStrangSplittingSetup")
    RETURN 1
    
  END SUBROUTINE Monodomain_ProblemStrangSplittingSetup
  
  !
  !================================================================================================================================
  !

  !>Sets up the linear Monodomain equations solution.
  SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP(PROBLEM,PROBLEM_SETUP,ERR,ERROR,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem to setup
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SOLVER_TYPE), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVERS_TYPE), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP",ERR,ERROR,*999)

    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVER_EQUATIONS)
    NULLIFY(SOLVERS)
    
    IF(ASSOCIATED(PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain equation problem.",err,error,*999)
      END IF
      IF(PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE.OR. &
          & PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)THEN
        SELECT CASE(PROBLEM_SETUP%setupType)
        CASE(PROBLEM_SETUP_INITIAL_TYPE)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Do nothing????
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Do nothing???
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_CONTROL_TYPE)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Set up a simple control loop
            CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,CONTROL_TIME_LOOP_TYPE,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Finish the control loops
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,ERR,ERROR,*999)            
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVERS_TYPE)
           !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
          SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Start the solvers creation
            CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_NUMBER_SET(SOLVERS,1,ERR,ERROR,*999)
            !Set the solver to be a first order dynamic solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,ERR,ERROR,*999)
            !Set solver defaults
            CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,ERR,ERROR,*999)
            CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,ERR,ERROR,*999)
            CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the solvers
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            !Finish the solvers creation
            CALL SOLVERS_CREATE_FINISH(SOLVERS,ERR,ERROR,*999)
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
                & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
           SELECT CASE(PROBLEM_SETUP%actionType)
          CASE(PROBLEM_SETUP_START_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            !Create the solver equations
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,ERR,ERROR,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,ERR,ERROR,*999)
          CASE(PROBLEM_SETUP_FINISH_ACTION)
            !Get the control loop
            CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
            CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,ERR,ERROR,*999)
            !Get the solver equations
            CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,ERR,ERROR,*999)
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,ERR,ERROR,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,ERR,ERROR,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,ERR,ERROR,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The action type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%actionType,"*",ERR,ERROR))// &
              & " for a setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
              & " is invalid for a monodomain equation."
            CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The setup type of "//TRIM(NUMBER_TO_VSTRING(PROBLEM_SETUP%setupType,"*",ERR,ERROR))// &
            & " is invalid for a Monodomain equation."
          CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
        END SELECT
      ELSE
        LOCAL_ERROR="The problem subtype of "//TRIM(NUMBER_TO_VSTRING(PROBLEM%SPECIFICATION(3),"*",ERR,ERROR))// &
          & " does not equal a Monodomain equation subtype."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      ENDIF
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
    
    EXITS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP")
    RETURN
999 ERRORSEXITS("MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_EQUATION_PROBLEM_SUBTYPE_SETUP

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE MONODOMAIN_PRE_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD, independentField

    ENTERS("MONODOMAIN_PRE_SOLVE",ERR,ERROR,*999)

    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<2) THEN
        CALL FlagError("Problem specification must at least two entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
        DEPENDENT_FIELD=>SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
          & equationsSets(1)%ptr%dependent%dependentField
        independentField=>SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
          & equationsSets(1)%ptr%INDEPENDENT%independentField

        CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 1,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,ERR,ERROR,*999)
        CALL Field_ParametersToFieldParametersCopy(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
          & 1,DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_PREVIOUS_VALUES_SET_TYPE,1,ERR,ERROR,*999) ! also to prev.

      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a monodomain problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("MONODOMAIN_PRE_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_PRE_SOLVE

  !
  !================================================================================================================================
  !

  !>Sets up the output type for a monodomain problem class.
  SUBROUTINE MONODOMAIN_POST_SOLVE(CONTROL_LOOP,SOLVER,ERR,ERROR,*)

    !Argument variables
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP !<A pointer to the control loop to solve.
    TYPE(SOLVER_TYPE), POINTER :: SOLVER !<A pointer to the solver
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    REAL(DP) :: TMPV, TMPA
    INTEGER(INTG) :: I
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD, MATERIAL_FIELD,  independentField, GEOMETRIC_FIELD

    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
 
    ENTERS("MONODOMAIN_POST_SOLVE",ERR,ERROR,*999)
  
    IF(ASSOCIATED(CONTROL_LOOP%PROBLEM)) THEN
      IF(.NOT.ALLOCATED(CONTROL_LOOP%PROBLEM%SPECIFICATION)) THEN
        CALL FlagError("Problem specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(CONTROL_LOOP%PROBLEM%SPECIFICATION,1)<3) THEN
        CALL FlagError("Problem specification must have three entries for a monodomain equation problem.",err,error,*999)
      END IF
      SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(2))
      CASE(PROBLEM_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE) ! CAN NOT GET EQN TYPE?! 


       ! Don't update in next step
        EQUATIONS=>SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping%equationsSets(1)%ptr%EQUATIONS
        vectorEquations=>equations%vectorEquations
        vectorMatrices=>vectorEquations%vectorMatrices
        dynamicMatrices=>vectorMatrices%dynamicMatrices
        stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
        dampingMatrix=>dynamicMatrices%matrices(2)%ptr
        rhsVector=>vectorMatrices%rhsVector
        stiffnessMatrix%updateMatrix = .FALSE.
        dampingMatrix%updateMatrix   = .FALSE.
        rhsVector%updateVector       = .FALSE.

       ! integrate   cell models
        GEOMETRIC_FIELD => SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
                           & equationsSets(1)%ptr%GEOMETRY%geometricField
        DEPENDENT_FIELD => SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
                           & equationsSets(1)%ptr%dependent%dependentField
        MATERIAL_FIELD  => SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
                           & equationsSets(1)%ptr%MATERIALS%materialsField
        independentField => SOLVER%SOLVERS%SOLVERS(1)%ptr%SOLVER_EQUATIONS%solverMapping% &
                           & equationsSets(1)%ptr%INDEPENDENT%independentField

        CALL Field_ParametersToFieldParametersCopy(DEPENDENT_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, &
         & 1, independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE, 1,ERR,ERROR,*999) ! dependent -> independent

        SELECT CASE(CONTROL_LOOP%PROBLEM%SPECIFICATION(3))
        CASE(PROBLEM_MONODOMAIN_BUENOOROVIO_SUBTYPE)
          CALL BUENO_OROVIO_INTEGRATE(independentField,MATERIAL_FIELD,&
          &  CONTROL_LOOP%timeLoop%currentTime-CONTROL_LOOP%timeLoop%timeIncrement, CONTROL_LOOP%timeLoop%currentTime&  ! from t-dt to t
          & ,ERR,ERROR,*999)
        CASE(PROBLEM_MONODOMAIN_TENTUSSCHER06_SUBTYPE)
           CALL TENTUSSCHER06_INTEGRATE(independentField,MATERIAL_FIELD,&
          &  CONTROL_LOOP%timeLoop%currentTime-CONTROL_LOOP%timeLoop%timeIncrement, CONTROL_LOOP%timeLoop%currentTime&  ! from t-dt to t
          & ,ERR,ERROR,*999)
        CASE DEFAULT
          CALL FlagError("Invalid cell model subtype",ERR,ERROR,*999)
        END SELECT

        DO I=1,independentField%DECOMPOSITION%DOMAIN(1)%ptr%TOPOLOGY%NODES%numberOfNodes
          !Default to version 1 of each derivative
          CALL Field_ParameterSetGetNode(independentField,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,1,TMPV,&  ! get local node?
               & ERR,ERROR,*999)
          IF (TMPV > 0) THEN
           !Default to version 1 of each derivative
           CALL Field_ParameterSetGetNode(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,7,TMPA,&
                & ERR,ERROR,*999) 
            IF(ABS(TMPA)<ZERO_TOLERANCE) THEN
              !Default to version 1 of each derivative
              CALL Field_ParameterSetUpdateNode(MATERIAL_FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VALUES_SET_TYPE,1,1,I,7,&
                   &CONTROL_LOOP%timeLoop%currentTime, ERR,ERROR,*999)  
            ENDIF      
          ENDIF
        ENDDO

        IF(MOD(CONTROL_LOOP%timeLoop%currentTime+1e-6,5.0)<1e-3) THEN
          WRITE(*,*)  'T=',CONTROL_LOOP%timeLoop%currentTime
        ENDIF

      CASE DEFAULT
        LOCAL_ERROR="Problem type "//TRIM(NUMBER_TO_VSTRING(CONTROL_LOOP%PROBLEM%SPECIFICATION(2),"*",ERR,ERROR))// &
          & " is not valid for a monodomain problem class."
        CALL FlagError(LOCAL_ERROR,ERR,ERROR,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",ERR,ERROR,*999)
    ENDIF
       
    EXITS("MONODOMAIN_POST_SOLVE")
    RETURN
999 ERRORSEXITS("MONODOMAIN_POST_SOLVE",ERR,ERROR)
    RETURN 1
  END SUBROUTINE MONODOMAIN_POST_SOLVE

  !
  !================================================================================================================================
  !

END MODULE MONODOMAIN_EQUATIONS_ROUTINES
