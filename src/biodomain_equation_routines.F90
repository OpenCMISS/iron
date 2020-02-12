!> \file
!> \author Chris Bradley
!> \brief This module handles all bioelectric domain equation routines.
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

!>This module handles all bioelectric domain equation routines.
MODULE BIODOMAIN_EQUATION_ROUTINES

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

#include "macros.h"

  IMPLICIT NONE

  PRIVATE

  !Module parameters

  !Module types

  !Module variables

  !Interfaces

  PUBLIC BIODOMAIN_CONTROL_LOOP_POST_LOOP
  
  PUBLIC Biodomain_EquationsSetSetup

  PUBLIC Biodomain_EquationsSetSolutionMethodSet
  
  PUBLIC Biodomain_EquationsSetSpecificationSet

  PUBLIC BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE
  
  PUBLIC BIODOMAIN_PRE_SOLVE
  
  PUBLIC BIODOMAIN_EQUATION_PROBLEM_SETUP

  PUBLIC Biodomain_ProblemSpecificationSet
  
CONTAINS

  !
  !================================================================================================================================
  !

  !>Runs after each control loop iteration
  SUBROUTINE BIODOMAIN_CONTROL_LOOP_POST_LOOP(CONTROL_LOOP,err,error,*)

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
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SOLVER_MAPPING_TYPE), POINTER :: SOLVER_MAPPING
    TYPE(SolversType), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: FILENAME,LOCAL_ERROR,METHOD
    INTEGER(INTG) :: OUTPUT_ITERATION_NUMBER,CURRENT_LOOP_ITERATION

    ENTERS("BIODOMAIN_CONTROL_LOOP_POST_LOOP",err,error,*999)

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
              !Get the solver. For Biodomain problems of any split the 2nd solver will contain the dependent field equation set
              CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)            
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              !Loop over the equations sets associated with the solver
              SOLVER_EQUATIONS=>SOLVER%solverEquations
              IF(ASSOCIATED(SOLVER_EQUATIONS)) THEN
                SOLVER_MAPPING=>SOLVER_EQUATIONS%solverMapping
                IF(ASSOCIATED(SOLVER_MAPPING)) THEN
                  DO equations_set_idx=1,SOLVER_MAPPING%numberOfEquationsSets
                    EQUATIONS_SET=>SOLVER_MAPPING%equationsSets(equations_set_idx)%ptr
                    IF(ASSOCIATED(EQUATIONS_SET)) THEN
                      DEPENDENT_FIELD=>EQUATIONS_SET%dependent%dependentField
                      NULLIFY(DEPENDENT_REGION)
                      CALL Field_RegionGet(DEPENDENT_FIELD,DEPENDENT_REGION,err,error,*999)
                      NULLIFY(PARENT_LOOP)
                      PARENT_LOOP=>CONTROL_LOOP%parentLoop
                      IF(ASSOCIATED(PARENT_LOOP)) THEN
                        !add the iteration number of the parent loop to the filename
                        NULLIFY(TIME_LOOP_PARENT)
                        TIME_LOOP_PARENT=>PARENT_LOOP%timeLoop
                        IF(ASSOCIATED(TIME_LOOP_PARENT)) THEN
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP_PARENT%outputNumber
                          CURRENT_LOOP_ITERATION=TIME_LOOP_PARENT%globalIterationNumber
                          FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%userNumber,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP_PARENT%globalIterationNumber,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP%iterationNumber,"*",err,error))
                        ELSE
                          OUTPUT_ITERATION_NUMBER=TIME_LOOP%outputNumber
                          CURRENT_LOOP_ITERATION=TIME_LOOP%globalIterationNumber
                          FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%userNumber,"*",err,error))// &
                            & "_"//TRIM(NumberToVString(TIME_LOOP%globalIterationNumber,"*",err,error))
                        ENDIF
                      ELSE
                        OUTPUT_ITERATION_NUMBER=TIME_LOOP%outputNumber
                        CURRENT_LOOP_ITERATION=TIME_LOOP%globalIterationNumber
                        FILENAME="Time_"//TRIM(NumberToVString(DEPENDENT_REGION%userNumber,"*",err,error))// &
                          & "_"//TRIM(NumberToVString(TIME_LOOP%globalIterationNumber,"*",err,error))
                      ENDIF
                      METHOD="FORTRAN"
                      IF(OUTPUT_ITERATION_NUMBER>0) THEN
                        IF(MOD(CURRENT_LOOP_ITERATION,OUTPUT_ITERATION_NUMBER)==0) THEN
                          CALL FIELD_IO_NODES_EXPORT(DEPENDENT_REGION%FIELDS,FILENAME,METHOD,err,error,*999)
                        ENDIF
                      ENDIF
                    ELSE
                      LOCAL_ERROR="Equations set is not associated for equations set index "// &
                        & TRIM(NumberToVString(equations_set_idx,"*",err,error))// &
                        & " in the solver mapping."
                      CALL FlagError(LOCAL_ERROR,err,error,*999)
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
        CASE(CONTROL_WHILE_LOOP_TYPE)
          !do nothing
        CASE(CONTROL_LOAD_INCREMENT_LOOP_TYPE)
          !do nothing
        CASE DEFAULT
          LOCAL_ERROR="The control loop type of "//TRIM(NumberToVString(CONTROL_LOOP%loopType,"*",err,error))// &
            & " is invalid."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ENDIF
    ELSE
      CALL FlagError("Control loop is not associated.",err,error,*999)
    ENDIF

    EXITS("BIODOMAIN_CONTROL_LOOP_POST_LOOP")
    RETURN
999 ERRORSEXITS("BIODOMAIN_CONTROL_LOOP_POST_LOOP",err,error)
    RETURN 1
    
  END SUBROUTINE BIODOMAIN_CONTROL_LOOP_POST_LOOP

  !
  !================================================================================================================================
  !

  !>Sets up the bioelectric domain equation type of a bioelectric equations set class.
  SUBROUTINE Biodomain_EquationsSetSetup(EQUATIONS_SET,EQUATIONS_SET_SETUP,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to setup a bioelectric domain equation on.
    TYPE(EquationsSetSetupType), INTENT(INOUT) :: EQUATIONS_SET_SETUP !<The equations set setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) :: component_idx,dimension_idx,DIMENSION_MULTIPLIER,GEOMETRIC_COMPONENT_NUMBER,GEOMETRIC_SCALING_TYPE, &
      & numberOfDimensions,NUMBER_OF_MATERIALS_COMPONENTS,GEOMETRIC_MESH_COMPONENT
    TYPE(DecompositionType), POINTER :: GEOMETRIC_DECOMPOSITION
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(EquationsSetMaterialsType), POINTER :: EQUATIONS_MATERIALS
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: EQUATIONS_SET_SPEC_TYPE,EQUATIONS_SET_SPEC_SUBTYPE
    
    ENTERS("Biodomain_EquationsSetSetup",err,error,*999)

    NULLIFY(EQUATIONS)
    NULLIFY(vectorMapping)
    NULLIFY(vectorMatrices)
    NULLIFY(GEOMETRIC_DECOMPOSITION)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a biodomain equation class.",err,error,*999)
      END IF
      EQUATIONS_SET_SPEC_TYPE=EQUATIONS_SET%SPECIFICATION(2)
      EQUATIONS_SET_SPEC_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      SELECT CASE(EQUATIONS_SET_SETUP%setupType)
      CASE(EQUATIONS_SET_SETUP_INITIAL_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL Biodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,EQUATIONS_SET_FEM_SOLUTION_METHOD, &
            & err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
!!Todo: CHECK VALID SETUP
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_GEOMETRY_TYPE)
        !\todo Check geometric dimension
      CASE(EQUATIONS_SET_SETUP_DEPENDENT_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
              !Create the auto created dependent field
              CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                & dependentField,err,error,*999)
              CALL Field_LabelSet(EQUATIONS_SET%dependent%dependentField,"Dependent Field",err,error,*999)
              CALL Field_TypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
              CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
              CALL Field_DecompositionSetAndLock(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_DECOMPOSITION, &
                & err,error,*999)
              CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%dependent%dependentField,EQUATIONS_SET%GEOMETRY% &
                & geometricField,err,error,*999)
              SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
              CASE(EQUATIONS_SET_NO_SUBTYPE)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%dependent%dependentField,2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE],err,error,*999)
              CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
                & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%dependent%dependentField,3,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE,"GeometryM3D",ERR, &
                  & ERROR,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE,3, &
                  & err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third equations set specification of "// &
                  & TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                  & " is not valid for a monodomain equation set."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Vm",err,error,*999)
              CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn", &
                & err,error,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                & FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & err,error,*999)
              CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & err,error,*999)
              !Default to the geometric interpolation setup
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                & GEOMETRIC_MESH_COMPONENT,err,error,*999)
              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                !Default the scaling to the geometric field scaling
                CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
                CALL Field_ScalingTypeSet(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ELSE
              !Check the user specified field
              CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
              CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
              SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
              CASE(EQUATIONS_SET_NO_SUBTYPE)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE], &
                  & err,error,*999)
              CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
                & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,3,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,3,err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="The third equations set specification of "// &
                  & TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                  & " is not valid for a monodomain equation set."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
              CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                & err,error,*999)
              CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
              SELECT CASE(EQUATIONS_SET%solutionMethod)
              CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                  & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
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
                LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                  & " is invalid."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            ENDIF
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
            CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
              IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
                !Create the auto created dependent field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%DEPENDENT% &
                  & dependentField,err,error,*999)
                CALL Field_LabelSet(EQUATIONS_SET%dependent%dependentField,"Dependent Field",err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DEPENDENT_TYPE,err,error,*999)
                CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_DecompositionSetAndLock(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%dependent%dependentField,EQUATIONS_SET%GEOMETRY% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%dependent%dependentField,4,err,error,*999)
!! \todo allow for no rhs variable and so eliminate one of the flux variables
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%dependent%dependentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,"Vm",err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,"dVm/dn", &
                  & err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE,"Phie", &
                  & err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,"dPhie/dn", &
                  & err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DimensionSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                  & err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                  & FIELD_DELVDELN_VARIABLE_TYPE,1,err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_V_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%dependent%dependentField,FIELD_DELVDELN_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                    & FIELD_DELUDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                    & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%dependent%dependentField, &
                    & FIELD_DELVDELN_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  !Default the scaling to the geometric field scaling
                  CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL Field_ScalingTypeSet(EQUATIONS_SET%dependent%dependentField,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
!! \todo allow for no rhs variable and so eliminate one of the flux variables
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
              IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
                CALL FlagError("The dependent field for the second bidomain equation cannot be auto-created. "// &
                  & "You must pass in the field from the first bidomain equation.",err,error,*999)
              ELSE
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DEPENDENT_TYPE,err,error,*999)
!! \todo allow for no rhs variable and so eliminate one of the flux variables
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_DELUDELN_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE,FIELD_DELVDELN_VARIABLE_TYPE],err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                  & err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELUDELN_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_DELVDELN_VARIABLE_TYPE,1, &
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                & " is invalid for a bidomain equations set type."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The equation set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
              & " is invalid for a biodomain equations set class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(EQUATIONS_SET%DEPENDENT%dependentFieldAutoCreated) THEN
            CALL Field_CreateFinish(EQUATIONS_SET%dependent%dependentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation"
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT

      CASE(EQUATIONS_SET_SETUP_INDEPENDENT_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
            CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
              & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
              IF(EQUATIONS_SET%INDEPENDENT%independentFieldAutoCreated) THEN
                !Create the auto created independent field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                  & independentField,err,error,*999)
                CALL Field_LabelSet(EQUATIONS_SET%INDEPENDENT%independentField,"Independent Field",err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_DecompositionSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,EQUATIONS_SET%GEOMETRY% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,4,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE,FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & "XB_state_variables",err,error,*999)
                ELSE                
                  CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & "Active_stress",err,error,*999)
                ENDIF
                CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE, &
                  & "sarcomere half length",err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE, &
                  & "contraction velocity",err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                  CALL Field_DimensionSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_SCALAR_DIMENSION_TYPE,err,error,*999)
                ENDIF
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & 1,err,error,*999)
                ELSEIF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & 6,err,error,*999)
                ELSEIF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                    & 4,err,error,*999)
                ENDIF
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,5, &
                  & err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE, &
                    & 4,err,error,*999)
                ELSE
                  CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE, &
                    & 3,err,error,*999)
                ENDIF
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE, &
                  & 6,err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,2, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,3, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,4, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,5, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,6, &
                    & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDIF
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,2, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,3, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,4, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,5, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE,2, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE,3, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U1_VARIABLE_TYPE, &
                    & 4,GEOMETRIC_MESH_COMPONENT,err,error,*999)
                ENDIF
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,2, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,3, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,4, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,5, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U2_VARIABLE_TYPE,6, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U_VARIABLE_TYPE,6,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U1_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U1_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U1_VARIABLE_TYPE,3,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                      & FIELD_U1_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,2,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U2_VARIABLE_TYPE,6,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  !Default the scaling to the geometric field scaling
                  CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL Field_ScalingTypeSet(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,4,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_U1_VARIABLE_TYPE,FIELD_U2_VARIABLE_TYPE],err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                  CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_SCALAR_DIMENSION_TYPE, &
                    & err,error,*999)
                ENDIF
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                ELSE IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,6,err,error,*999)
                ELSEIF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4,err,error,*999)
                ENDIF
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,5,err,error,*999)
                IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,4,err,error,*999)
                ELSE
                  CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,3,err,error,*999)
                ENDIF
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,6,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,3, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,5, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,6, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF  
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,2, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,3, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,4, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,5, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,2, &
                    & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,3, &
                    & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  IF(EQUATIONS_SET_SPEC_SUBTYPE==EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE) THEN
                    CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U1_VARIABLE_TYPE,4, &
                      & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  ENDIF
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,2, &
                    & FIELD_CONSTANT_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,3, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,4, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,5, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U2_VARIABLE_TYPE,6, &
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE)
              IF(EQUATIONS_SET%INDEPENDENT%independentFieldAutoCreated) THEN
                !Create the auto created independent field
                CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_SET%INDEPENDENT% &
                  & independentField,err,error,*999)
                CALL Field_LabelSet(EQUATIONS_SET%INDEPENDENT%independentField,"Independent Field",err,error,*999)
                CALL Field_TypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_INDEPENDENT_TYPE, &
                  & err,error,*999)
                CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
                CALL Field_DecompositionSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_DECOMPOSITION, &
                  & err,error,*999)
                CALL Field_GeometricFieldSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,EQUATIONS_SET%GEOMETRY% &
                  & geometricField,err,error,*999)
                CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,2,err,error,*999)
                CALL Field_VariableTypesSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,[FIELD_U_VARIABLE_TYPE, &
                  & FIELD_V_VARIABLE_TYPE],err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                  & "Active_stress",err,error,*999)
                CALL Field_VariableLabelSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
                  & "fibre_info",err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                  & FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE, &
                  & FIELD_INTG_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE, &
                  & 1,err,error,*999)
                CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,5, &
                  & err,error,*999)
                !Default to the geometric interpolation setup
                CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_U_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,1, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,2, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,3, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,4, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                CALL Field_ComponentMeshComponentSet(EQUATIONS_SET%INDEPENDENT%independentField,FIELD_V_VARIABLE_TYPE,5, &
                  & GEOMETRIC_MESH_COMPONENT,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_U_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,1,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,2,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,3,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,4,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationSetAndLock(EQUATIONS_SET%INDEPENDENT%independentField, &
                    & FIELD_V_VARIABLE_TYPE,5,FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  !Default the scaling to the geometric field scaling
                  CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
                  CALL Field_ScalingTypeSet(EQUATIONS_SET%INDEPENDENT%independentField,GEOMETRIC_SCALING_TYPE,err,error,*999)
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ELSE
                !Check the user specified field
                CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_GENERAL_TYPE,err,error,*999)
                CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
                CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,2,err,error,*999)
                CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE,FIELD_V_VARIABLE_TYPE],ERR, &
                  & ERROR,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
                CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,FIELD_INTG_TYPE,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1,err,error,*999)
                CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,5,err,error,*999)
                SELECT CASE(EQUATIONS_SET%solutionMethod)
                CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,1, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,2, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,3, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,4, &
                    & FIELD_NODE_BASED_INTERPOLATION,err,error,*999)
                  CALL Field_ComponentInterpolationCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_V_VARIABLE_TYPE,5, &
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
                  LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
                    & " is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE DEFAULT
              LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                " is not implemented for an equations set setup independent type."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            LOCAL_ERROR="Equations set setup independent type is not implemented for an equations set bidomain equation type"
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The equation set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
              & " is invalid for a biodomain equations set class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          IF(EQUATIONS_SET%INDEPENDENT%independentFieldAutoCreated) THEN
            CALL Field_CreateFinish(EQUATIONS_SET%INDEPENDENT%independentField,err,error,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation"
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT

      CASE(EQUATIONS_SET_SETUP_MATERIALS_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(EQUATIONS_SET,err,error,*999)
          CALL EquationsSet_AssertMaterialsIsCreated(EQUATIONS_SET,err,error,*999)
          EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
          IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
            !Create the auto created materials field
            CALL Field_CreateStart(EQUATIONS_SET_SETUP%fieldUserNumber,EQUATIONS_SET%REGION,EQUATIONS_MATERIALS% &
              & materialsField,err,error,*999)
            CALL Field_LabelSet(EQUATIONS_MATERIALS%materialsField,"Materials Field",err,error,*999)
            CALL Field_TypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_DecompositionGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_DECOMPOSITION,err,error,*999)
            CALL Field_DecompositionSetAndLock(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_DECOMPOSITION, &
              & err,error,*999)
            CALL Field_GeometricFieldSetAndLock(EQUATIONS_MATERIALS%materialsField,EQUATIONS_SET%GEOMETRY% &
              & geometricField,err,error,*999)
            CALL Field_NumberOfVariablesSetAndLock(EQUATIONS_MATERIALS%materialsField,1,err,error,*999)
            CALL Field_VariableTypesSetAndLock(EQUATIONS_MATERIALS%materialsField,[FIELD_U_VARIABLE_TYPE], &
              & err,error,*999)
            CALL Field_VariableLabelSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE,"Materials",err,error,*999)
            CALL Field_DimensionSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_VECTOR_DIMENSION_TYPE,err,error,*999)
            CALL Field_DataTypeSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
              & FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
              & numberOfDimensions,err,error,*999)
            IF(EQUATIONS_SET_SPEC_TYPE==EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE) THEN
              !Monodomain. Materials field components are 2 plus one for each dimension i.e., Am, Cm and \sigma
              NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions+2
              DIMENSION_MULTIPLIER=1
            ELSE
              !Bidomain. Materials field components are 2 plus two for each dimension i.e., Am, C, \sigma_i and \sigma_e
              NUMBER_OF_MATERIALS_COMPONENTS=2*numberOfDimensions+2
              DIMENSION_MULTIPLIER=2
            ENDIF
            !Set the number of materials components
            CALL Field_NumberOfComponentsSetAndLock(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
              & NUMBER_OF_MATERIALS_COMPONENTS,err,error,*999)
            !Default the Am and Cm materials components to the first component geometric interpolation with const interpolation
            CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
              & 1,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
            DO component_idx=1,2
              CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & component_idx,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
            ENDDO !components_idx
            CALL Field_ComponentLabelSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE,1,"Am", &
              & err,error,*999)
            CALL Field_ComponentLabelSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE,2,"Cm", &
              & err,error,*999)
            !Default the \sigma materials components to the geometric interpolation setup with constant interpolation
            DO component_idx=1,numberOfDimensions
              CALL Field_ComponentMeshComponentGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & component_idx,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
              DO dimension_idx=1,DIMENSION_MULTIPLIER
                CALL Field_ComponentMeshComponentSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 2+component_idx+(dimension_idx-1)*numberOfDimensions,GEOMETRIC_COMPONENT_NUMBER,err,error,*999)
                CALL Field_ComponentInterpolationSet(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                  & 2+component_idx+(dimension_idx-1)*numberOfDimensions,FIELD_CONSTANT_INTERPOLATION,err,error,*999)
              ENDDO !dimension_idx
            ENDDO !component_idx
            !Default the field scaling to that of the geometric field
            CALL Field_ScalingTypeGet(EQUATIONS_SET%GEOMETRY%geometricField,GEOMETRIC_SCALING_TYPE,err,error,*999)
            CALL Field_ScalingTypeSet(EQUATIONS_MATERIALS%materialsField,GEOMETRIC_SCALING_TYPE,err,error,*999)
          ELSE
            !Check the user specified field
            CALL Field_TypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_MATERIAL_TYPE,err,error,*999)
            CALL Field_DependentTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_INDEPENDENT_TYPE,err,error,*999)
            CALL Field_NumberOfVariablesCheck(EQUATIONS_SET_SETUP%FIELD,1,err,error,*999)
            CALL Field_VariableTypesCheck(EQUATIONS_SET_SETUP%FIELD,[FIELD_U_VARIABLE_TYPE],err,error,*999)
            CALL Field_DimensionCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_VECTOR_DIMENSION_TYPE, &
              & err,error,*999)
            CALL Field_DataTypeCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,FIELD_DP_TYPE,err,error,*999)
            CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
              & numberOfDimensions,err,error,*999)
            SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
            CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE,EQUATIONS_SET_MONODOMAIN_STRANG_SPLITTING_EQUATION_TYPE)
              !Monodomain. Materials field components are 2 plus one for each dimension i.e., Am, Cm and \sigma
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,numberOfDimensions+2, &
                & err,error,*999)
            CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
              !Bidomain. Materials field components are 2 plus two for each dimension i.e., Am, C, \sigma_i and \sigma_e
              CALL Field_NumberOfComponentsCheck(EQUATIONS_SET_SETUP%FIELD,FIELD_U_VARIABLE_TYPE,2*numberOfDimensions+2, &
                & err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The equations set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
                & " is invalid for a bioelectrics class."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ENDIF
       CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          EQUATIONS_MATERIALS=>EQUATIONS_SET%MATERIALS
          IF(ASSOCIATED(EQUATIONS_MATERIALS)) THEN
            IF(EQUATIONS_MATERIALS%materialsFieldAutoCreated) THEN
              !Finish creating the materials field
              CALL Field_CreateFinish(EQUATIONS_MATERIALS%materialsField,err,error,*999)
              !Set the default values for the materials field
              CALL Field_NumberOfComponentsGet(EQUATIONS_SET%GEOMETRY%geometricField,FIELD_U_VARIABLE_TYPE, &
                & numberOfDimensions,err,error,*999)
              IF(EQUATIONS_SET_SPEC_TYPE==EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE) THEN
                !Monodomain. Materials field components are 2 plus one for each dimension i.e., Am, Cm and \sigma
                NUMBER_OF_MATERIALS_COMPONENTS=numberOfDimensions+2
                DIMENSION_MULTIPLIER=1
              ELSE
                !Bidomain. Materials field components are 2 plus two for each dimension i.e., Am, C, \sigma_i and \sigma_e
                NUMBER_OF_MATERIALS_COMPONENTS=2*numberOfDimensions+2
                DIMENSION_MULTIPLIER=2
              ENDIF
              !First set Am
              CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,1,200.0_DP,err,error,*999)
              !Now set Cm
              CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                & FIELD_VALUES_SET_TYPE,2,0.0025_DP,err,error,*999)
              !Now set the sigmas to be 1.0
              DO component_idx=1,numberOfDimensions
                DO dimension_idx=1,DIMENSION_MULTIPLIER
                  CALL Field_ComponentValuesInitialise(EQUATIONS_MATERIALS%materialsField,FIELD_U_VARIABLE_TYPE, &
                    & FIELD_VALUES_SET_TYPE,2+component_idx+(dimension_idx-1)*numberOfDimensions,1.0_DP,err,error,*999)
                ENDDO !dimension_idx
              ENDDO !component_idx
            ENDIF
          ELSE
            CALL FlagError("Equations set materials is not associated.",err,error,*999)
          ENDIF
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_SOURCE_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          !Do nothing
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          !Do nothing
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_ANALYTIC_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          CALL FlagError("Not implemented.",err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_SETUP_EQUATIONS_TYPE)
        SELECT CASE(EQUATIONS_SET_SETUP%actionType)
        CASE(EQUATIONS_SET_SETUP_START_ACTION)
          CALL EquationsSet_AssertDependentIsFinished(EQUATIONS_SET,err,error,*999)
          CALL EquationsSet_AssertMaterialsIsFInished(EQUATIONS_SET,err,error,*999)
          !Create the equations
          CALL Equations_CreateStart(EQUATIONS_SET,EQUATIONS,err,error,*999)
          CALL Equations_LinearityTypeSet(EQUATIONS,EQUATIONS_LINEAR,err,error,*999)
          SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
          CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
            CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
          CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
            CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
              CALL Equations_TimeDependenceTypeSet(EQUATIONS,EQUATIONS_STATIC,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The third equations set specification of "// &
                & TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                & " is invalid for a bidomain equation."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The second equations set specification of "// &
              & TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
              & " is invalid for a bioelectrics class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SETUP_FINISH_ACTION)
          SELECT CASE(EQUATIONS_SET%solutionMethod)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            !Finish the creation of the equations
            CALL EquationsSet_EquationsGet(EQUATIONS_SET,EQUATIONS,err,error,*999)
            CALL Equations_CreateFinish(EQUATIONS,err,error,*999)
            NULLIFY(vectorEquations)
            CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
            !Create the equations mapping.
            CALL EquationsMapping_VectorCreateStart(vectorEquations,FIELD_DELUDELN_VARIABLE_TYPE,vectorMapping,err,error,*999)
            SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
            CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
              CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
              CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
              CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
            CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
              SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
              CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
                CALL EquationsMappingVector_DynamicMatricesSet(vectorMapping,.TRUE.,.TRUE.,err,error,*999)
                CALL EquationsMappingVector_DynamicVariableTypeSet(vectorMapping,FIELD_U_VARIABLE_TYPE,err,error,*999)
                CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELUDELN_VARIABLE_TYPE,err,error,*999)
              CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
                CALL EquationsMappingVector_NumberOfLinearMatricesSet(vectorMapping,1,err,error,*999)
                CALL EquationsMappingVector_LinearMatricesVariableTypesSet(vectorMapping,[FIELD_V_VARIABLE_TYPE], &
                  & err,error,*999)
                CALL EquationsMappingVector_RHSVariableTypeSet(vectorMapping,FIELD_DELVDELN_VARIABLE_TYPE,err,error,*999)
              CASE DEFAULT
                LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                  & " is invalid for a bidomain equation type."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The equations set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
                & " is invalid for a bioelectrics class."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            CALL EquationsMapping_VectorCreateFinish(vectorMapping,err,error,*999)
            !Create the equations matrices
            CALL EquationsMatrices_VectorCreateStart(vectorEquations,vectorMatrices,err,error,*999)
            SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
            CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
              !Set up matrix storage and structure
              IF(equations%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
              ELSE
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
            CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
              SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
              CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
              !Set up matrix storage and structure
              IF(equations%lumpingType==EQUATIONS_LUMPED_MATRICES) THEN
                !Set up lumping
                CALL EquationsMatricesVector_DynamicLumpingTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_UNLUMPED,EQUATIONS_MATRIX_LUMPED],err,error,*999)
                CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                  & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_DIAGONAL_STORAGE_TYPE], &
                  & err,error,*999)
                CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                  & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_DIAGONAL_STRUCTURE],err,error,*999)
              ELSE
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES) 
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE,DISTRIBUTED_MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatricesVector_DynamicStorageTypeSet(vectorMatrices, &
                    & [DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE,DISTRIBUTED_MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_DynamicStructureTypeSet(vectorMatrices, &
                    & [EQUATIONS_MATRIX_FEM_STRUCTURE,EQUATIONS_MATRIX_FEM_STRUCTURE],err,error,*999)                  
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              ENDIF
              CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
                SELECT CASE(equations%sparsityType)
                CASE(EQUATIONS_MATRICES_FULL_MATRICES)
                  CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_BLOCK_STORAGE_TYPE], &
                    & err,error,*999)
                CASE(EQUATIONS_MATRICES_SPARSE_MATRICES)
                  CALL EquationsMatricesVector_LinearStorageTypeSet(vectorMatrices,[MATRIX_COMPRESSED_ROW_STORAGE_TYPE], &
                    & err,error,*999)
                  CALL EquationsMatricesVector_LinearStructureTypeSet(vectorMatrices,[EQUATIONS_MATRIX_FEM_STRUCTURE], &
                    & err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The equations matrices sparsity type of "// &
                    & TRIM(NumberToVString(equations%sparsityType,"*",err,error))//" is invalid."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
                  & " is invalid for a bidomain equation type."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The equations set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
                & " is invalid for a bioelectrics class."
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
            LOCAL_ERROR="The solution method of "//TRIM(NumberToVString(EQUATIONS_SET%solutionMethod,"*",err,error))// &
              & " is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(EQUATIONS_SET_SETUP%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
    
    EXITS("Biodomain_EquationsSetSetup")
    RETURN
999 ERRORSEXITS("Biodomain_EquationsSetSetup",err,error)
    RETURN 1
    
  END SUBROUTINE Biodomain_EquationsSetSetup

  !
  !================================================================================================================================
  !

  !>Sets/changes the solution method for a bioelectric domain equation type of an bioelectrics equations set class.
  SUBROUTINE Biodomain_EquationsSetSolutionMethodSet(EQUATIONS_SET,SOLUTION_METHOD,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to set the solution method for
    INTEGER(INTG), INTENT(IN) :: SOLUTION_METHOD !<The solution method to set
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    INTEGER(INTG) :: EQUATIONS_SET_SPEC_TYPE,EQUATIONS_SET_SPEC_SUBTYPE
    
    ENTERS("Biodomain_EquationsSetSolutionMethodSet",err,error,*999)
    
    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
        CALL FlagError("Equations set specification is not allocated.",err,error,*999)
      ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
        CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
      END IF
      EQUATIONS_SET_SPEC_TYPE=EQUATIONS_SET%SPECIFICATION(2)
      EQUATIONS_SET_SPEC_SUBTYPE=EQUATIONS_SET%SPECIFICATION(3)
      SELECT CASE(EQUATIONS_SET_SPEC_TYPE)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)        
        SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
        CASE(EQUATIONS_SET_NO_SUBTYPE,EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
          SELECT CASE(SOLUTION_METHOD)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
            LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
            & " is not valid for a bioelectric monodomain equation type of an bioelectrics equations set class."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)        
        SELECT CASE(EQUATIONS_SET_SPEC_SUBTYPE)
        CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)        
          SELECT CASE(SOLUTION_METHOD)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
            LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
          SELECT CASE(SOLUTION_METHOD)
          CASE(EQUATIONS_SET_FEM_SOLUTION_METHOD)
            EQUATIONS_SET%solutionMethod=EQUATIONS_SET_FEM_SOLUTION_METHOD
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
            LOCAL_ERROR="The specified solution method of "//TRIM(NumberToVString(SOLUTION_METHOD,"*",err,error))//" is invalid."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="Equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_SUBTYPE,"*",err,error))// &
            & " is not valid for a bioelectric bidomain equation type of an bioelectrics equations set class."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="Equations set type of "//TRIM(NumberToVString(EQUATIONS_SET_SPEC_TYPE,"*",err,error))// &
          & " is not valid for a bioelectrics equations set class."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("Biodomain_EquationsSetSolutionMethodSet")
    RETURN
999 ERRORS("Biodomain_EquationsSetSolutionMethodSet",err,error)
    EXITS("Biodomain_EquationsSetSolutionMethodSet")
    RETURN 1
    
  END SUBROUTINE Biodomain_EquationsSetSolutionMethodSet

  !
  !================================================================================================================================
  !

  !>Sets the equation specification for a bioelectric domain equation type of a bioelectric equations set class.
  SUBROUTINE Biodomain_EquationsSetSpecificationSet(equationsSet,specification,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: equationsSet !<A pointer to the equations set to set the specification for
    INTEGER(INTG), INTENT(IN) :: specification(:) !<The equations set specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: equationsSetType,equationsSetSubtype

    ENTERS("Biodomain_EquationsSetSpecificationSet",err,error,*999)

    IF(ASSOCIATED(equationsSet)) THEN
      IF(SIZE(specification,1)/=3) THEN
        CALL FlagError("Equations set specification must have three entries for a biodomain equation type equations set.", &
          & err,error,*999)
      END IF
      equationsSetType=specification(2)
      equationsSetSubtype=specification(3)
      SELECT CASE(equationsSetType)
      CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_NO_SUBTYPE)
          !ok
        CASE(EQUATIONS_SET_1D3D_MONODOMAIN_ELASTICITY_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE, &
          & EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
          & EQUATIONS_SET_1D3D_MONODOMAIN_ACTIVE_STRAIN_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third equations set specification of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a monodomain type of a bioelectric equations set."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
        SELECT CASE(equationsSetSubtype)
        CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
          !ok
        CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
          !ok
        CASE DEFAULT
          localError="The third equations set specification of "//TRIM(NumberToVstring(equationsSetSubtype,"*",err,error))// &
            & " is not valid for a bidomain equation type of a bioelectric equations set class."
          CALL FlagError(localError,err,error,*999)
        END SELECT
      CASE DEFAULT
        localError="The second equations set specification of "//TRIM(NumberToVstring(equationsSetType,"*",err,error))// &
          & " is not valid for a bioelectric equations set."
        CALL FlagError(localError,err,error,*999)
      END SELECT
      !Set full specification
      IF(ALLOCATED(equationsSet%specification)) THEN
        CALL FlagError("Equations set specification is already allocated.",err,error,*999)
      ELSE
        ALLOCATE(equationsSet%specification(3),stat=err)
        IF(err/=0) CALL FlagError("Could not allocate equations set specification.",err,error,*999)
      END IF
      equationsSet%specification(1:3)=[EQUATIONS_SET_BIOELECTRICS_CLASS,equationsSetType,equationsSetSubtype]
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    END IF

    EXITS("Biodomain_EquationsSetSpecificationSet")
    RETURN
999 ERRORS("Biodomain_EquationsSetSpecificationSet",err,error)
    EXITS("Biodomain_EquationsSetSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Biodomain_EquationsSetSpecificationSet

  !
  !================================================================================================================================
  !

  !>Performs pre-solve actions for mono- and bi-domain problems.
  SUBROUTINE BIODOMAIN_PRE_SOLVE(SOLVER,err,error,*)

    !Argument variables
    TYPE(SolverType), POINTER :: SOLVER !<A pointer to the solver to perform the pre-solve actions for.
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    REAL(DP) :: CURRENT_TIME,TIME_INCREMENT
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP
    TYPE(ProblemType), POINTER :: PROBLEM
    TYPE(SolversType), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    ENTERS("BIODOMAIN_PRE_SOLVE",err,error,*999)

    IF(ASSOCIATED(SOLVER)) THEN
      SOLVERS=>SOLVER%solvers
      IF(ASSOCIATED(SOLVERS)) THEN
        CONTROL_LOOP=>SOLVERS%controlLoop
        IF(ASSOCIATED(CONTROL_LOOP)) THEN
          CALL CONTROL_LOOP_CURRENT_TIMES_GET(CONTROL_LOOP,CURRENT_TIME,TIME_INCREMENT,err,error,*999)
          PROBLEM=>CONTROL_LOOP%PROBLEM
          IF(ASSOCIATED(PROBLEM)) THEN
            IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
              CALL FlagError("Problem specification is not allocated.",err,error,*999)
            ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
              CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
            END IF
            SELECT CASE(PROBLEM%SPECIFICATION(2))
            CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE)
                SELECT CASE(SOLVER%globalNumber)
                CASE(1)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT,err,error,*999)
                CASE(2)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
                    & " is invalid for a Gudunov split monodomain problem."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE(PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
                SELECT CASE(SOLVER%globalNumber)
                CASE(1)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT/2.0_DP,err,error,*999)
                CASE(2)
                  !Do nothing
                CASE(3)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME+TIME_INCREMENT/2.0_DP,CURRENT_TIME+TIME_INCREMENT, &
                    & err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
                    & " is invalid for a Strang split monodomain problem."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is invalid for a monodomain problem type."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
                SELECT CASE(SOLVER%globalNumber)
                CASE(1)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT,err,error,*999)
                CASE(2)
                  !Do nothing
                CASE(3)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
                    & " is invalid for a Gudunov split bidomain problem."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
                SELECT CASE(SOLVER%globalNumber)
                CASE(1)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT/2.0_DP,err,error,*999)
                CASE(2)
                  !Do nothing
                CASE(3)
                  !Do nothing
                CASE(4)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME+TIME_INCREMENT/2.0_DP,CURRENT_TIME+TIME_INCREMENT, &
                    & err,error,*999)
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
                    & " is invalid for a Gudunov split bidomain problem."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is invalid for a bidomain problem type."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE(PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE)
              SELECT CASE(PROBLEM%SPECIFICATION(3))
              CASE(PROBLEM_GUDUNOV_MONODOMAIN_SIMPLE_ELASTICITY_SUBTYPE,PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE, &
                & PROBLEM_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE,EQUATIONS_SET_MONODOMAIN_ELASTICITY_VELOCITY_SUBTYPE, &
                & PROBLEM_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE)
                SELECT CASE(SOLVER%globalNumber)
                CASE(1)
                  CALL SOLVER_DAE_TIMES_SET(SOLVER,CURRENT_TIME,CURRENT_TIME+TIME_INCREMENT,err,error,*999)
                CASE(2)
                  !Do nothing
                CASE DEFAULT
                  LOCAL_ERROR="The solver global number of "//TRIM(NumberToVString(SOLVER%globalNumber,"*",err,error))// &
                    & " is invalid for a bioelectrics finite elasticity problem."
                  CALL FlagError(LOCAL_ERROR,err,error,*999)
                END SELECT
              CASE DEFAULT
                LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
                  & " is invalid for a monodomain problem type."
                CALL FlagError(LOCAL_ERROR,err,error,*999)
              END SELECT
            CASE DEFAULT
              LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))// &
                & " is invalid."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          ELSE
            CALL FlagError("Control loop problem is not associated.",err,error,*999)
          ENDIF
        ELSE
          CALL FlagError("Solvers control loop is not associated.",err,error,*999)
        ENDIF
      ELSE
        CALL FlagError("Solver solvers is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Solver is not associated.",err,error,*999)
    ENDIF
       
    EXITS("BIODOMAIN_PRE_SOLVE")
    RETURN
999 ERRORSEXITS("BIODOMAIN_PRE_SOLVE",err,error)
    RETURN 1
    
  END SUBROUTINE BIODOMAIN_PRE_SOLVE

  !
  !================================================================================================================================
  !
 
  !>Sets up the bioelectric domain problem.
  SUBROUTINE BIODOMAIN_EQUATION_PROBLEM_SETUP(PROBLEM,PROBLEM_SETUP,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER :: PROBLEM !<A pointer to the problem set to setup a bioelectric domain equation on.
    TYPE(ProblemSetupType), INTENT(INOUT) :: PROBLEM_SETUP !<The problem setup information
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    TYPE(CellMLEquationsType), POINTER :: CELLML_EQUATIONS
    TYPE(ControlLoopType), POINTER :: CONTROL_LOOP,CONTROL_LOOP_ROOT
    TYPE(SolverType), POINTER :: SOLVER
    TYPE(SOLVER_EQUATIONS_TYPE), POINTER :: SOLVER_EQUATIONS
    TYPE(SolversType), POINTER :: SOLVERS
    TYPE(VARYING_STRING) :: LOCAL_ERROR

    NULLIFY(CELLML_EQUATIONS)
    NULLIFY(CONTROL_LOOP)
    NULLIFY(SOLVER)
    NULLIFY(SOLVERS)
    NULLIFY(SOLVER_EQUATIONS)
    
    ENTERS("BIODOMAIN_EQUATION_PROBLEM_SETUP",err,error,*999)

    IF(ASSOCIATED(PROBLEM)) THEN
      SELECT CASE(PROBLEM_SETUP%setupType)
      CASE(PROBLEM_SETUP_INITIAL_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Do nothing????
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Do nothing????
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CONTROL_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Set up a time control loop
          CALL CONTROL_LOOP_CREATE_START(PROBLEM,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_TYPE_SET(CONTROL_LOOP,CONTROL_TIME_LOOP_TYPE,err,error,*999)
          CALL CONTROL_LOOP_LABEL_SET(CONTROL_LOOP,"Time Loop",err,error,*999)
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Finish the control loops
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_CREATE_FINISH(CONTROL_LOOP,err,error,*999)            
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric domain equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVERS_TYPE)
        !Get the control loop
        CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
        CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Start the solvers creation
          CALL SOLVERS_CREATE_START(CONTROL_LOOP,SOLVERS,err,error,*999)
          IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE)
              CALL Solvers_NumberOfSolversSet(SOLVERS,2,err,error,*999)
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE Solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE(PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
              CALL Solvers_NumberOfSolversSet(SOLVERS,3,err,error,*999)
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"First ODE solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the third solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Second ODE solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is invalid for a monodomain problem type of a bioelectric problem class."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
              CALL Solvers_NumberOfSolversSet(SOLVERS,3,err,error,*999)
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"ODE solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the third solver to be a linear solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Elliptic solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
              CALL Solvers_NumberOfSolversSet(SOLVERS,4,err,error,*999)
              !Set the first solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"First ODE solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the second solver to be a dynamic solver 
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DYNAMIC_TYPE,err,error,*999)
              CALL SOLVER_DYNAMIC_ORDER_SET(SOLVER,SOLVER_DYNAMIC_FIRST_ORDER,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Parabolic solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_DYNAMIC_DEGREE_SET(SOLVER,SOLVER_DYNAMIC_FIRST_DEGREE,err,error,*999)
              CALL SOLVER_DYNAMIC_SCHEME_SET(SOLVER,SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,err,error,*999)
              CALL SOLVER_DYNAMIC_RESTART_SET(SOLVER,.TRUE.,err,error,*999)
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the third solver to be a differential-algebraic equations solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_DAE_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Second ODE solver",err,error,*999)
              !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_CMISS_LIBRARY,err,error,*999)
              !Set the fourth solver to be a linear solver
              NULLIFY(SOLVER)
              CALL SOLVERS_SOLVER_GET(SOLVERS,4,SOLVER,err,error,*999)
              CALL SOLVER_TYPE_SET(SOLVER,SOLVER_LINEAR_TYPE,err,error,*999)
              CALL SOLVER_LABEL_SET(SOLVER,"Elliptic solver",err,error,*999)
             !Set solver defaults
              CALL SOLVER_LIBRARY_TYPE_SET(SOLVER,SOLVER_PETSC_LIBRARY,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))// &
                & " is invalid for a monodomain problem type of a bioelectric problem class."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
          CASE DEFAULT
            LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))// &
              & " is invalid for a bioelectric problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the solvers
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          !Finish the solvers creation
          CALL SOLVERS_CREATE_FINISH(SOLVERS,err,error,*999)
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_SOLVER_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          !Get the solver
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
            !Create the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
            !Create the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_FIRST_ORDER_DYNAMIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
            !Create the solver equations for the elliptic solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,4,SOLVER,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))//  &
                & " is invalid for a bidomain problem type."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            CALL SOLVER_EQUATIONS_CREATE_START(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            CALL SOLVER_EQUATIONS_LINEARITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_LINEAR,err,error,*999)
            CALL SOLVER_EQUATIONS_TIME_DEPENDENCE_TYPE_SET(SOLVER_EQUATIONS,SOLVER_EQUATIONS_STATIC,err,error,*999)
            CALL SOLVER_EQUATIONS_SPARSITY_TYPE_SET(SOLVER_EQUATIONS,SOLVER_SPARSE_MATRICES,err,error,*999)
          CASE DEFAULT
            LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))//  &
              & " is invalid for a bioelectric problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
            !Get the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
            !Get the solver equations for the second (parabolic) solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,2,SOLVER,err,error,*999)
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
            !Get the solver equations for the elliptic solver
            NULLIFY(SOLVER)
            NULLIFY(SOLVER_EQUATIONS)
            SELECT CASE(PROBLEM%SPECIFICATION(3))
            CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
            CASE(PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
              CALL SOLVERS_SOLVER_GET(SOLVERS,4,SOLVER,err,error,*999)
            CASE DEFAULT
              LOCAL_ERROR="The problem subtype of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(3),"*",err,error))//  &
                & " is invalid for a bidomain problem type."
              CALL FlagError(LOCAL_ERROR,err,error,*999)
            END SELECT
            CALL SOLVER_SOLVER_EQUATIONS_GET(SOLVER,SOLVER_EQUATIONS,err,error,*999)
            !Finish the solver equations creation
            CALL SOLVER_EQUATIONS_CREATE_FINISH(SOLVER_EQUATIONS,err,error,*999)             
          CASE DEFAULT
            LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))//  &
              & " is invalid for a bioelectric problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE(PROBLEM_SETUP_CELLML_EQUATIONS_TYPE)
        SELECT CASE(PROBLEM_SETUP%actionType)
        CASE(PROBLEM_SETUP_START_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          !Get the solver
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL CellMLEquations_CreateStart(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            IF(PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE) THEN
              !Create the CellML equations for the second DAE solver
              NULLIFY(SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL CellMLEquations_CreateStart(SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            ENDIF
          CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
            !Create the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL CellMLEquations_CreateStart(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Set the time dependence
            CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
            !Set the linearity
            CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            IF(PROBLEM%SPECIFICATION(3)==PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE) THEN
              !Create the CellML equations for the second DAE solver
              NULLIFY(SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,4,SOLVER,err,error,*999)
              CALL CellMLEquations_CreateStart(SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Set the time dependence
              CALL CellMLEquations_TimeDependenceTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_DYNAMIC,err,error,*999)
              !Set the linearity
              CALL CellMLEquations_LinearityTypeSet(CELLML_EQUATIONS,CELLML_EQUATIONS_NONLINEAR,err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))//  &
              & " is invalid for a bioelectric problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE(PROBLEM_SETUP_FINISH_ACTION)
          !Get the control loop
          CONTROL_LOOP_ROOT=>PROBLEM%controlLoop
          CALL CONTROL_LOOP_GET(CONTROL_LOOP_ROOT,CONTROL_LOOP_NODE,CONTROL_LOOP,err,error,*999)
          CALL CONTROL_LOOP_SOLVERS_GET(CONTROL_LOOP,SOLVERS,err,error,*999)
          IF(.NOT.ALLOCATED(PROBLEM%SPECIFICATION)) THEN
            CALL FlagError("Problem specification is not allocated.",err,error,*999)
          ELSE IF(SIZE(PROBLEM%SPECIFICATION,1)<3) THEN
            CALL FlagError("Problem specification must have three entries for a biodomain problem.",err,error,*999)
          END IF
          SELECT CASE(PROBLEM%SPECIFICATION(2))
          CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(CELLML_EQUATIONS,err,error,*999)
            IF(PROBLEM%SPECIFICATION(3)==PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE) THEN
              !Get the CellML equations for the second DAE solver
              NULLIFY(SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,3,SOLVER,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CellMLEquations_CreateFinish(CELLML_EQUATIONS,err,error,*999)
            ENDIF
          CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
            !Get the CellML equations for the first DAE solver
            CALL SOLVERS_SOLVER_GET(SOLVERS,1,SOLVER,err,error,*999)
            CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
            !Finish the CellML equations creation
            CALL CellMLEquations_CreateFinish(CELLML_EQUATIONS,err,error,*999)
            IF(PROBLEM%SPECIFICATION(3)==PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE) THEN
              !Get the CellML equations for the second DAE solver
              NULLIFY(SOLVER)
              NULLIFY(CELLML_EQUATIONS)
              CALL SOLVERS_SOLVER_GET(SOLVERS,4,SOLVER,err,error,*999)
              CALL SOLVER_CELLML_EQUATIONS_GET(SOLVER,CELLML_EQUATIONS,err,error,*999)
              !Finish the CellML equations creation
              CALL CellMLEquations_CreateFinish(CELLML_EQUATIONS,err,error,*999)
            ENDIF
          CASE DEFAULT
            LOCAL_ERROR="The problem type of "//TRIM(NumberToVString(PROBLEM%SPECIFICATION(2),"*",err,error))//  &
              & " is invalid for a bioelectric problem class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT            
        CASE DEFAULT
          LOCAL_ERROR="The action type of "//TRIM(NumberToVString(PROBLEM_SETUP%actionType,"*",err,error))// &
            & " for a setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
            & " is invalid for a bioelectric equation."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      CASE DEFAULT
        LOCAL_ERROR="The setup type of "//TRIM(NumberToVString(PROBLEM_SETUP%setupType,"*",err,error))// &
          & " is invalid for a bioelectric domain equation."
        CALL FlagError(LOCAL_ERROR,err,error,*999)
      END SELECT
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    ENDIF
      
    EXITS("BIODOMAIN_EQUATION_PROBLEM_SETUP")
    RETURN
999 ERRORSEXITS("BIODOMAIN_EQUATION_PROBLEM_SETUP",err,error)
    RETURN 1
  END SUBROUTINE BIODOMAIN_EQUATION_PROBLEM_SETUP
  
  !
  !================================================================================================================================
  !

  !>Sets the problem specification for a bioelectric domain problem class.
  SUBROUTINE Biodomain_ProblemSpecificationSet(problem,problemSpecification,err,error,*)

    !Argument variables
    TYPE(ProblemType), POINTER, INTENT(IN) :: problem !<A pointer to the problem to set the specification for
    INTEGER(INTG), INTENT(IN) :: problemSpecification(:) !<The problem specification to set
    INTEGER(INTG), INTENT(OUT) :: err !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: error !<The error string
    !Local Variables
    TYPE(VARYING_STRING) :: localError
    INTEGER(INTG) :: problemType
    INTEGER(INTG) :: problemSubtype

    ENTERS("Biodomain_ProblemSpecificationSet",err,error,*999)

    IF(ASSOCIATED(problem)) THEN
      IF(SIZE(problemSpecification,1)==3) THEN
        problemType=problemSpecification(2)
        problemSubtype=problemSpecification(3)
        SELECT CASE(problemType)
        CASE(PROBLEM_MONODOMAIN_EQUATION_TYPE)
          SELECT CASE(problemSubtype)
          CASE(PROBLEM_MONODOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
              & PROBLEM_MONODOMAIN_STRANG_SPLIT_SUBTYPE)
            !ok
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
              & " is not valid for a monodomain type of a bioelectric problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE(PROBLEM_BIDOMAIN_EQUATION_TYPE)
          SELECT CASE(problemSubtype)
          CASE(PROBLEM_BIDOMAIN_GUDUNOV_SPLIT_SUBTYPE, &
              & PROBLEM_BIDOMAIN_STRANG_SPLIT_SUBTYPE)
            !ok
          CASE DEFAULT
            localError="The third problem specification of "//TRIM(NumberToVstring(problemSubtype,"*",err,error))// &
              & " is not valid for a bidomain type of a bioelectric problem."
            CALL FlagError(localError,err,error,*999)
          END SELECT
        CASE DEFAULT
          localError="The second problem specification of "//TRIM(NumberToVstring(problemType,"*",err,error))// &
            & " is not valid for a bioelectric problem."
          CALL FlagError(localError,err,error,*999)
        END SELECT
        IF(ALLOCATED(problem%specification)) THEN
          CALL FlagError("Problem specification is already allocated.",err,error,*999)
        ELSE
          ALLOCATE(problem%specification(3),stat=err)
          IF(err/=0) CALL FlagError("Could not allocate problem specification.",err,error,*999)
        END IF
        problem%specification(1:3)=[PROBLEM_BIOELECTRICS_CLASS,problemType,problemSubtype]
      ELSE
        CALL FlagError("Biodomain problem specification must have three entries.",err,error,*999)
      END IF
    ELSE
      CALL FlagError("Problem is not associated.",err,error,*999)
    END IF

    EXITS("Biodomain_ProblemSpecificationSet")
    RETURN
999 ERRORS("Biodomain_ProblemSpecificationSet",err,error)
    EXITS("Biodomain_ProblemSpecificationSet")
    RETURN 1
    
  END SUBROUTINE Biodomain_ProblemSpecificationSet

  !
  !================================================================================================================================
  !

  !>Calculates the element stiffness matrices and RHS for a bioelectric domain equation finite element equations set.
  SUBROUTINE BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE(EQUATIONS_SET,ELEMENT_NUMBER,err,error,*)

    !Argument variables
    TYPE(EquationsSetType), POINTER :: EQUATIONS_SET !<A pointer to the equations set to perform the finite element calculations on
    INTEGER(INTG), INTENT(IN) :: ELEMENT_NUMBER !<The element number to calculate
    INTEGER(INTG), INTENT(OUT) :: ERR !<The error code
    TYPE(VARYING_STRING), INTENT(OUT) :: ERROR !<The error string
    !Local Variables
    INTEGER(INTG) FIELD_VAR_TYPE,mh,mhs,ms,ng,nh,nhs,ni,nj,ns
    LOGICAL :: USE_FIBRES
    REAL(DP) :: CONDUCTIVITY(3,3),DPHIDX(3,64),RWG,SUM
    TYPE(BasisType), POINTER :: DEPENDENT_BASIS,GEOMETRIC_BASIS,FIBRE_BASIS
    TYPE(EquationsType), POINTER :: equations
    TYPE(EquationsMappingVectorType), POINTER :: vectorMapping
    TYPE(EquationsMappingDynamicType), POINTER :: dynamicMapping
    TYPE(EquationsMatricesVectorType), POINTER :: vectorMatrices
    TYPE(EquationsMatricesDynamicType), POINTER :: dynamicMatrices
    TYPE(EquationsMatricesRHSType), POINTER :: rhsVector
    TYPE(EquationsMatrixType), POINTER :: dampingMatrix,stiffnessMatrix
    TYPE(EquationsVectorType), POINTER :: vectorEquations
    TYPE(FieldType), POINTER :: DEPENDENT_FIELD,GEOMETRIC_FIELD,FIBRE_FIELD,materialsField
    TYPE(FieldVariableType), POINTER :: FIELD_VARIABLE,GEOMETRIC_VARIABLE
    TYPE(QuadratureSchemeType), POINTER :: QUADRATURE_SCHEME
    TYPE(VARYING_STRING) :: LOCAL_ERROR
    
    ENTERS("BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE",err,error,*999)

    IF(ASSOCIATED(EQUATIONS_SET)) THEN
      EQUATIONS=>EQUATIONS_SET%EQUATIONS
      IF(ASSOCIATED(EQUATIONS)) THEN
        NULLIFY(vectorEquations)
        CALL Equations_VectorEquationsGet(equations,vectorEquations,err,error,*999)
        DEPENDENT_FIELD=>equations%interpolation%dependentField
        GEOMETRIC_FIELD=>equations%interpolation%geometricField
        materialsField=>equations%interpolation%materialsField
        FIBRE_FIELD=>equations%interpolation%fibreField
        USE_FIBRES=ASSOCIATED(FIBRE_FIELD)
        vectorMapping=>vectorEquations%vectorMapping
        vectorMatrices=>vectorEquations%vectorMatrices
        DEPENDENT_BASIS=>DEPENDENT_FIELD%DECOMPOSITION%DOMAIN(DEPENDENT_FIELD%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        GEOMETRIC_BASIS=>GEOMETRIC_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        GEOMETRIC_VARIABLE=>GEOMETRIC_FIELD%variableTypeMap(FIELD_U_VARIABLE_TYPE)%ptr
        IF(USE_FIBRES) FIBRE_BASIS=>FIBRE_FIELD%DECOMPOSITION%DOMAIN(GEOMETRIC_FIELD%decomposition%meshComponentNumber)%ptr% &
          & TOPOLOGY%ELEMENTS%ELEMENTS(ELEMENT_NUMBER)%BASIS
        QUADRATURE_SCHEME=>DEPENDENT_BASIS%QUADRATURE%quadratureSchemeMap(BASIS_DEFAULT_QUADRATURE_SCHEME)%ptr
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
          & geometricInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations%interpolation% &
          & materialsInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        IF(USE_FIBRES) CALL Field_InterpolationParametersElementGet(FIELD_VALUES_SET_TYPE,ELEMENT_NUMBER,equations% &
          & interpolation%fibreInterpParameters(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
        
        IF(.NOT.ALLOCATED(EQUATIONS_SET%SPECIFICATION)) THEN
          CALL FlagError("Equations set specification is not allocated.",err,error,*999)
        ELSE IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<2) THEN
          CALL FlagError("Equations set specification does not have a type set.",err,error,*999)
        END IF
        SELECT CASE(EQUATIONS_SET%SPECIFICATION(2))
        CASE(EQUATIONS_SET_MONODOMAIN_EQUATION_TYPE)
          
          dynamicMatrices=>vectorMatrices%dynamicMatrices
          stiffnessMatrix=>dynamicMatrices%matrices(1)%ptr
          dampingMatrix=>dynamicMatrices%matrices(2)%ptr
          rhsVector=>vectorMatrices%rhsVector
          dynamicMapping=>vectorMapping%dynamicMapping
          FIELD_VARIABLE=>dynamicMapping%equationsMatrixToVarMaps(1)%variable
          FIELD_VAR_TYPE=FIELD_VARIABLE%variableType

          IF(stiffnessMatrix%updateMatrix.OR.dampingMatrix%updateMatrix.OR.rhsVector%updateVector) THEN
            !Loop over gauss points
            DO ng=1,QUADRATURE_SCHEME%numberOfGauss
              CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & geometricInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolatedPointMetricsCalculate(GEOMETRIC_BASIS%numberOfXi,equations%interpolation% &
                & geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              CALL Field_InterpolateGauss(NO_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                & materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              IF(USE_FIBRES) THEN
                CALL Field_InterpolateGauss(FIRST_PART_DERIV,BASIS_DEFAULT_QUADRATURE_SCHEME,ng,equations%interpolation% &
                  & fibreInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
                CALL Field_InterpolatedPointMetricsCalculate(FIBRE_BASIS%numberOfXi,equations%interpolation% &
                  & fibreInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr,err,error,*999)
              ENDIF
              !Calculate RWG.
              RWG=equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%jacobian* &
                & QUADRATURE_SCHEME%gaussWeights(ng)
              !Calculate the conductivity tensor
              CONDUCTIVITY=0.0_DP
              IF(USE_FIBRES) THEN
                !Calculate the conductivity tensor in fibre coordinates
                CALL FlagError("Not implemented.",err,error,*999)
              ELSE
                !Use the conductivity tensor in geometric coordinates
                DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                  CONDUCTIVITY(nj,nj)=equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(nj+2,1)
                ENDDO !nj
              ENDIF
              !Compute basis dPhi/dx terms
              DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                  DPHIDX(nj,ms)=0.0_DP
                  DO ni=1,DEPENDENT_BASIS%numberOfXi
                    DPHIDX(nj,ms)=DPHIDX(nj,ms)+ &
                      & QUADRATURE_SCHEME%gaussBasisFunctions(ms,PARTIAL_DERIVATIVE_FIRST_DERIVATIVE_MAP(ni),ng)* &
                      & equations%interpolation%geometricInterpPointMetrics(FIELD_U_VARIABLE_TYPE)%ptr%dXidX(ni,nj)
                  ENDDO !ni
                ENDDO !ms
              ENDDO !nj            
              !Loop over field components
              mhs=0          
              DO mh=1,FIELD_VARIABLE%numberOfComponents
                !Loop over element rows
                DO ms=1,DEPENDENT_BASIS%numberOfElementParameters
                  mhs=mhs+1
                  nhs=0
                  !Loop over element columns
                  DO nh=1,FIELD_VARIABLE%numberOfComponents
                    DO ns=1,DEPENDENT_BASIS%numberOfElementParameters
                      nhs=nhs+1
                      SUM=0.0_DP
                      IF(stiffnessMatrix%updateMatrix) THEN
                        DO ni=1,GEOMETRIC_VARIABLE%numberOfComponents
                          DO nj=1,GEOMETRIC_VARIABLE%numberOfComponents
                            SUM=SUM+CONDUCTIVITY(ni,nj)*DPHIDX(ni,mhs)*DPHIDX(nj,nhs)
                          ENDDO !nj
                        ENDDO !ni
                        IF((equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)<ZERO_TOLERANCE)&
                          & .OR. (equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1) &
                          & <ZERO_TOLERANCE)) THEN
                          LOCAL_ERROR="The value of the surface area to volume ratio or the capacitance is below zero tolerance"
                          CALL FlagError(LOCAL_ERROR,err,error,*999)
                        ENDIF
                        stiffnessMatrix%elementMatrix%matrix(mhs,nhs)=stiffnessMatrix%elementMatrix%matrix(mhs,nhs)+SUM*RWG/ &
                          & equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(1,1)/ &
                          & equations%interpolation%materialsInterpPoint(FIELD_U_VARIABLE_TYPE)%ptr%VALUES(2,1)
                      ENDIF
                      IF(dampingMatrix%updateMatrix) THEN
                        dampingMatrix%elementMatrix%matrix(mhs,nhs)=dampingMatrix%elementMatrix%matrix(mhs,nhs)+ &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ms,NO_PART_DERIV,ng)* &
                          & QUADRATURE_SCHEME%gaussBasisFunctions(ns,NO_PART_DERIV,ng)*RWG
                      ENDIF
                    ENDDO !ns
                  ENDDO !nh
                  IF(rhsVector%updateVector) rhsVector%elementVector%vector(mhs)=0.0_DP
                ENDDO !ms
              ENDDO !mh
            ENDDO !ng
          ENDIF
        CASE(EQUATIONS_SET_BIDOMAIN_EQUATION_TYPE)
          IF(SIZE(EQUATIONS_SET%SPECIFICATION,1)<3) THEN
            CALL FlagError("Equations set specification does not have a subtype set.",err,error,*999)
          END IF
          SELECT CASE(EQUATIONS_SET%SPECIFICATION(3))
          CASE(EQUATIONS_SET_FIRST_BIDOMAIN_SUBTYPE)
          CASE(EQUATIONS_SET_SECOND_BIDOMAIN_SUBTYPE)
          CASE DEFAULT
            LOCAL_ERROR="The equations set subtype of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(3),"*",err,error))// &
              & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
            CALL FlagError(LOCAL_ERROR,err,error,*999)
          END SELECT
        CASE DEFAULT
          LOCAL_ERROR="The equations set type of "//TRIM(NumberToVString(EQUATIONS_SET%SPECIFICATION(2),"*",err,error))// &
            & " is not valid for a bioelectric domain type of a bioelectrics equations set class."
          CALL FlagError(LOCAL_ERROR,err,error,*999)
        END SELECT
      ELSE
        CALL FlagError("Equations set equations is not associated.",err,error,*999)
      ENDIF
    ELSE
      CALL FlagError("Equations set is not associated.",err,error,*999)
    ENDIF
       
    EXITS("BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE")
    RETURN
999 ERRORSEXITS("BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE",err,error)
    RETURN 1
  END SUBROUTINE BIODOMAIN_EQUATION_FINITE_ELEMENT_CALCULATE

  !
  !================================================================================================================================
  !
   
END MODULE BIODOMAIN_EQUATION_ROUTINES
